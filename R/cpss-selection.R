.cpss_feature_names <- function(X) {
  feature_names <- colnames(X)
  source <- "column_names"
  if (is.null(feature_names)) {
    feature_names <- paste0("x", seq_len(ncol(X)))
    source <- "generated_positionally"
  }
  feature_names <- as.character(feature_names)
  if (length(feature_names) != ncol(X) || anyNA(feature_names) ||
      any(!nzchar(feature_names)) || anyDuplicated(feature_names)) {
    stop("X column names must be unique, nonmissing, and nonempty.",
         call. = FALSE)
  }
  list(names = feature_names, source = source)
}

.cpss_resolve_index_spec <- function(
    index, p, feature_names, name, default = integer()) {
  if (is.null(index)) return(as.integer(default))
  if (is.character(index)) {
    if (!length(index) || anyNA(index) || any(!nzchar(index)) ||
        anyDuplicated(index)) {
      stop(name, " must contain unique nonmissing feature names.",
           call. = FALSE)
    }
    resolved <- match(index, feature_names)
    if (anyNA(resolved)) {
      stop(
        name, " contains feature names absent from X: ",
        paste(index[is.na(resolved)], collapse = ", "),
        call. = FALSE
      )
    }
    return(as.integer(resolved))
  }
  .validate_indices(index, p, name)
}

.cpss_partialize_mandatory <- function(
    X, y, mandatory_core, candidate_set, rank_tolerance) {
  y_centered <- y - mean(y)
  candidate_design <- X[, candidate_set, drop = FALSE]
  candidate_centered <- sweep(
    candidate_design, 2L, colMeans(candidate_design), "-"
  )
  candidate_norm <- sqrt(colSums(candidate_centered^2))
  if (!length(mandatory_core)) {
    relative_norm <- as.numeric(candidate_norm > 0)
    return(list(
      X = candidate_centered,
      y = y_centered,
      rank = 0L,
      condition_number = 1,
      relative_candidate_norm = relative_norm,
      eligible = is.finite(candidate_norm) &
        candidate_norm > sqrt(.Machine$double.eps)
    ))
  }

  mandatory_design <- X[, mandatory_core, drop = FALSE]
  mandatory_centered <- sweep(
    mandatory_design, 2L, colMeans(mandatory_design), "-"
  )
  decomposition <- svd(
    mandatory_centered,
    nu = min(nrow(mandatory_centered), ncol(mandatory_centered)),
    nv = 0L
  )
  if (!length(decomposition$d) || !is.finite(decomposition$d[1L]) ||
      decomposition$d[1L] <= 0) {
    stop("The centered mandatory core is rank deficient.", call. = FALSE)
  }
  cutoff <- rank_tolerance * max(dim(mandatory_centered)) *
    decomposition$d[1L]
  rank <- sum(decomposition$d > cutoff)
  if (rank != length(mandatory_core)) {
    stop(
      "The centered mandatory core is not full column rank at ",
      "rank_tolerance.", call. = FALSE
    )
  }
  basis <- decomposition$u[, seq_len(rank), drop = FALSE]
  residualized_X <- candidate_centered -
    basis %*% crossprod(basis, candidate_centered)
  residualized_y <- y_centered - basis %*% crossprod(basis, y_centered)
  residualized_norm <- sqrt(colSums(residualized_X^2))
  relative_norm <- residualized_norm / pmax(
    candidate_norm, sqrt(.Machine$double.eps)
  )
  list(
    X = residualized_X,
    y = as.numeric(residualized_y),
    rank = as.integer(rank),
    condition_number = decomposition$d[1L] /
      decomposition$d[rank],
    relative_candidate_norm = relative_norm,
    eligible = is.finite(relative_norm) &
      relative_norm > sqrt(.Machine$double.eps)
  )
}

.cpss_complementary_halves <- function(n, strata = NULL, seed = NULL) {
  if (is.null(strata)) strata <- rep("all", n)
  strata <- as.character(strata)
  if (length(strata) != n || anyNA(strata) || any(!nzchar(strata))) {
    stop("strata must contain one nonmissing label per observation.",
         call. = FALSE)
  }
  if (any(table(strata) < 2L)) {
    stop("Every CPSS stratum must contain at least two observations.",
         call. = FALSE)
  }
  .with_seed(seed, function() {
    first <- integer()
    second <- integer()
    offset <- 0L
    for (level in sort(unique(strata), method = "radix")) {
      index <- sample(which(strata == level), replace = FALSE)
      n_first <- floor(length(index) / 2L)
      if (length(index) %% 2L == 1L && offset %% 2L == 1L) {
        n_first <- n_first + 1L
      }
      n_first <- max(1L, min(length(index) - 1L, n_first))
      first <- c(first, index[seq_len(n_first)])
      second <- c(second, index[(n_first + 1L):length(index)])
      offset <- offset + 1L
    }
    list(first = sort(first), second = sort(second))
  })
}

.cpss_choose_budget_path <- function(beta, lambda, budget, tolerance) {
  support_size <- colSums(abs(beta) > tolerance)
  usable <- which(support_size > 0L & support_size <= budget)
  if (!length(usable)) {
    return(list(index = NA_integer_, support = integer(),
                support_size = support_size))
  }
  target <- max(support_size[usable])
  candidates <- usable[support_size[usable] == target]
  selected <- candidates[which.max(lambda[candidates])]
  list(
    index = as.integer(selected),
    support = which(abs(beta[, selected]) > tolerance),
    support_size = support_size
  )
}

.cpss_base_fit <- function(
    X, y, selector, budget, path_points, lambda_min_ratio,
    coefficient_tolerance, lasso_maxit, mcp_gamma, mcp_max_iter) {
  x_center <- colMeans(X)
  X_centered <- sweep(X, 2L, x_center, "-")
  x_scale <- sqrt(colMeans(X_centered^2))
  valid <- which(is.finite(x_scale) &
                   x_scale > sqrt(.Machine$double.eps))
  output_beta <- numeric(ncol(X))
  if (!length(valid)) {
    return(list(
      support = integer(), beta = output_beta, lambda = NA_real_,
      selected_path_index = NA_integer_,
      selected_path_locally_convex = NA, converged = TRUE,
      warnings = character(), iterations = 0L,
      path_support_size = integer()
    ))
  }
  y_centered <- y - mean(y)
  y_scale <- sqrt(mean(y_centered^2))
  if (!is.finite(y_scale) || y_scale <= sqrt(.Machine$double.eps)) {
    stop("A CPSS half-sample has a degenerate response.", call. = FALSE)
  }
  X_standardized <- sweep(
    X_centered[, valid, drop = FALSE], 2L, x_scale[valid], "/"
  )
  y_standardized <- y_centered / y_scale
  warnings_seen <- character()
  fit <- withCallingHandlers(
    if (identical(selector, "lasso")) {
      glmnet_arguments <- list(
        x = X_standardized, y = y_standardized, family = "gaussian",
        alpha = 1, intercept = FALSE, standardize = FALSE,
        nlambda = path_points, lambda.min.ratio = lambda_min_ratio
      )
      glmnet_limits <- list(
        dfmax = as.integer(max(2L * budget, budget + 5L)),
        maxit = lasso_maxit
      )
      if ("control" %in% names(formals(glmnet::glmnet))) {
        glmnet_arguments$control <- glmnet_limits
      } else {
        glmnet_arguments <- c(glmnet_arguments, glmnet_limits)
      }
      do.call(glmnet::glmnet, glmnet_arguments)
    } else {
      if (!requireNamespace("grpreg", quietly = TRUE)) {
        stop(
          "selector = \"mcp\" requires the optional package 'grpreg'.",
          call. = FALSE
        )
      }
      path_limit <- as.integer(max(2L * budget, budget + 5L))
      grpreg::grpreg(
        X = X_standardized, y = y_standardized, family = "gaussian",
        group = seq_len(ncol(X_standardized)), penalty = "grMCP",
        gamma = mcp_gamma, nlambda = path_points,
        lambda.min = lambda_min_ratio, log.lambda = TRUE, alpha = 1,
        eps = 1e-4, max.iter = mcp_max_iter,
        dfmax = path_limit, gmax = path_limit,
        group.multiplier = rep(1, ncol(X_standardized)),
        warn = FALSE, returnX = FALSE
      )
    },
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (identical(selector, "lasso")) {
    beta_path <- as.matrix(fit$beta)
    lambda <- as.numeric(fit$lambda)
    iterations <- as.integer(fit$npasses %||% NA_integer_)
    converged <- is.na(iterations) || iterations < lasso_maxit
    convex_min <- NA_integer_
  } else {
    beta_path <- as.matrix(fit$beta[-1L, , drop = FALSE])
    lambda <- as.numeric(fit$lambda)
    iterations <- as.integer(fit$iter %||% NA_integer_)
    if (!length(lambda) || nrow(beta_path) != length(valid) ||
        ncol(beta_path) != length(lambda) ||
        length(iterations) != length(lambda)) {
      stop("The grpreg MCP path is malformed.", call. = FALSE)
    }
    converged <- length(iterations) > 0L &&
      all(is.finite(iterations)) &&
      sum(iterations) < mcp_max_iter
    convex_min <- NA_integer_
  }
  choice <- .cpss_choose_budget_path(
    beta_path, lambda, budget, coefficient_tolerance
  )
  support <- integer()
  selected_lambda <- NA_real_
  locally_convex <- NA
  if (!is.na(choice$index)) {
    support <- valid[choice$support]
    output_beta[support] <- beta_path[choice$support, choice$index]
    selected_lambda <- lambda[choice$index]
    if (identical(selector, "mcp") && !is.na(convex_min)) {
      locally_convex <- choice$index <= convex_min
    }
  }
  list(
    support = support,
    beta = output_beta,
    lambda = selected_lambda,
    selected_path_index = choice$index,
    selected_path_locally_convex = locally_convex,
    converged = isTRUE(converged),
    warnings = unique(warnings_seen),
    iterations = iterations,
    path_support_size = choice$support_size
  )
}

#' Select a stable high-dimensional submodel core with CPSS
#'
#' `cpss_select_core()` implements complementary-pairs stability selection
#' (CPSS) with a size-budgeted LASSO or MCP base selector. Selection must be
#' performed on data independent of any sample later used for the
#' conditional-Gaussian exact-null test. CPSS itself is due to Shah and
#' Samworth (2013); this function supplies the core-selection stage used by
#' the package's CPSS--SM endpoint.
#'
#' @param X Numeric selection-sample design matrix.
#' @param y Numeric selection-sample response.
#' @param selector Base selector, either `"lasso"` or `"mcp"`.
#' @param complementary_pairs Number of complementary half-sample pairs.
#' @param base_selection_size Maximum support size of each base fit.
#' @param stability_threshold Inclusion-frequency threshold. The usual
#'   Meinshausen--Buehlmann PFER expression is only reported above one half.
#' @param seed Optional reproducible resampling seed.
#' @param path_points Number of penalty-path points per base fit.
#' @param lambda_min_ratio Path lower-bound control. It is passed as
#'   `lambda.min.ratio` to `glmnet` for LASSO and as the `lambda.min` fraction
#'   of `lambda.max` to `grpreg` for the singleton-group MCP sensitivity path.
#' @param coefficient_tolerance Numerical support threshold.
#' @param lasso_maxit Maximum `glmnet` iterations.
#' @param mcp_gamma MCP concavity parameter.
#' @param mcp_max_iter Maximum total `grpreg` iterations over the MCP path.
#' @param no_top_k_fallback Must be `TRUE`. An empty thresholded core is
#'   returned without a data-dependent top-k replacement. If a mandatory core
#'   is supplied, an empty CPSS extension leaves the mandatory core unchanged.
#' @param strata Optional resampling strata with one entry per row.
#' @param mandatory_core Optional one-based indices or feature names that are
#'   always retained. Within every half-sample they are centered and projected
#'   out by a checked economy SVD before the base selector is fitted. Mandatory
#'   variables are outside the CPSS frequency and PFER family.
#' @param candidate_set Optional one-based indices or feature names eligible
#'   for CPSS selection. The default is every predictor outside
#'   `mandatory_core`.
#' @param rank_tolerance Relative SVD tolerance for mandatory-core
#'   residualization.
#'
#' @return A `cpss_core_selection` object containing the selected indices,
#'   selection frequencies, feature identities, and base-fit diagnostics.
#' @references
#' Meinshausen, N. and Buehlmann, P. (2010). Stability selection.
#' *Journal of the Royal Statistical Society: Series B*, 72, 417--473.
#' \doi{10.1111/j.1467-9868.2010.00740.x}
#'
#' Shah, R. D. and Samworth, R. J. (2013). Variable selection with error
#' control: another look at stability selection. *Journal of the Royal
#' Statistical Society: Series B*, 75, 55--80.
#' \doi{10.1111/j.1467-9868.2011.01034.x}
#' @export
cpss_select_core <- function(
    X,
    y,
    selector = c("lasso", "mcp"),
    complementary_pairs = 50L,
    base_selection_size = 20L,
    stability_threshold = 0.60,
    seed = NULL,
    path_points = 100L,
    lambda_min_ratio = 0.01,
    coefficient_tolerance = 1e-8,
    lasso_maxit = 100000L,
    mcp_gamma = 3,
    mcp_max_iter = 100000L,
    no_top_k_fallback = TRUE,
    strata = NULL,
    mandatory_core = NULL,
    candidate_set = NULL,
    rank_tolerance = 1e-10) {
  call <- match.call()
  selector <- match.arg(selector)
  X <- .validate_numeric_matrix(X, "X")
  y <- .validate_response(y, nrow(X))
  feature_state <- .cpss_feature_names(X)
  n <- nrow(X)
  p <- ncol(X)
  mandatory_core <- .cpss_resolve_index_spec(
    mandatory_core, p, feature_state$names, "mandatory_core"
  )
  candidate_set <- .cpss_resolve_index_spec(
    candidate_set, p, feature_state$names, "candidate_set",
    default = setdiff(seq_len(p), mandatory_core)
  )
  if (!length(candidate_set)) {
    stop("candidate_set must contain at least one optional predictor.",
         call. = FALSE)
  }
  overlap <- intersect(mandatory_core, candidate_set)
  if (length(overlap)) {
    stop("mandatory_core and candidate_set must be disjoint.",
         call. = FALSE)
  }
  rank_tolerance <- .check_positive_scalar(
    rank_tolerance, "rank_tolerance"
  )
  whole_scale <- sqrt(colMeans(sweep(X, 2L, colMeans(X), "-")^2))
  used_set <- c(mandatory_core, candidate_set)
  if (any(!is.finite(whole_scale[used_set])) ||
      any(whole_scale[used_set] <= sqrt(.Machine$double.eps))) {
    stop("Every mandatory or candidate predictor must have positive variance.",
         call. = FALSE)
  }
  global_partial <- .cpss_partialize_mandatory(
    X, y, mandatory_core, candidate_set, rank_tolerance
  )
  eligible_candidate_set <- candidate_set[global_partial$eligible]
  ineligible_candidate_set <- candidate_set[!global_partial$eligible]
  if (!length(eligible_candidate_set)) {
    stop(
      "No candidate remains after residualizing against mandatory_core.",
      call. = FALSE
    )
  }
  complementary_pairs <- as.integer(complementary_pairs)
  base_selection_size <- as.integer(base_selection_size)
  path_points <- as.integer(path_points)
  lasso_maxit <- as.integer(lasso_maxit)
  mcp_max_iter <- as.integer(mcp_max_iter)
  if (length(complementary_pairs) != 1L || is.na(complementary_pairs) ||
      complementary_pairs < 1L) {
    stop("complementary_pairs must be a positive integer.", call. = FALSE)
  }
  if (length(base_selection_size) != 1L || is.na(base_selection_size) ||
      base_selection_size < 1L ||
      base_selection_size > length(eligible_candidate_set)) {
    stop(
      "base_selection_size must lie between one and the number of eligible ",
      "candidate predictors.", call. = FALSE
    )
  }
  if (length(path_points) != 1L || is.na(path_points) || path_points < 2L) {
    stop("path_points must be an integer of at least two.", call. = FALSE)
  }
  if (length(stability_threshold) != 1L ||
      !is.finite(stability_threshold) || stability_threshold <= 0 ||
      stability_threshold > 1) {
    stop("stability_threshold must lie in (0, 1].", call. = FALSE)
  }
  lambda_min_ratio <- .check_positive_scalar(
    lambda_min_ratio, "lambda_min_ratio"
  )
  if (lambda_min_ratio >= 1) {
    stop("lambda_min_ratio must be below one.", call. = FALSE)
  }
  coefficient_tolerance <- .check_positive_scalar(
    coefficient_tolerance, "coefficient_tolerance", allow_zero = TRUE
  )
  mcp_gamma <- .check_positive_scalar(mcp_gamma, "mcp_gamma")
  if (mcp_gamma <= 1) stop("mcp_gamma must exceed one.", call. = FALSE)
  if (!isTRUE(no_top_k_fallback)) {
    stop("no_top_k_fallback must be TRUE for the strict CPSS contract.",
         call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (length(seed) != 1L || is.na(seed)) {
      stop("seed must be NULL or one finite integer.", call. = FALSE)
    }
  }

  counts <- numeric(p)
  strength <- numeric(p)
  diagnostics <- vector("list", 2L * complementary_pairs)
  position <- 0L
  for (pair in seq_len(complementary_pairs)) {
    pair_seed <- if (is.null(seed)) NULL else
      as.integer((as.numeric(seed) + 1009 * pair) %%
                   (.Machine$integer.max - 1) + 1)
    halves <- .cpss_complementary_halves(n, strata, pair_seed)
    for (half in seq_len(2L)) {
      position <- position + 1L
      rows <- if (half == 1L) halves$first else halves$second
      started <- proc.time()[[3L]]
      half_partial <- .cpss_partialize_mandatory(
        X[rows, , drop = FALSE], y[rows], mandatory_core,
        eligible_candidate_set, rank_tolerance
      )
      fit <- .cpss_base_fit(
        half_partial$X, half_partial$y, selector,
        base_selection_size, path_points, lambda_min_ratio,
        coefficient_tolerance, lasso_maxit, mcp_gamma, mcp_max_iter
      )
      if (length(fit$support)) {
        selected_original <- eligible_candidate_set[fit$support]
        counts[selected_original] <- counts[selected_original] + 1
        strength[eligible_candidate_set] <-
          strength[eligible_candidate_set] + abs(fit$beta)
      }
      diagnostics[[position]] <- data.frame(
        selector = paste0("CPSS-", toupper(selector)),
        pair = pair,
        half = half,
        pair_seed = pair_seed %||% NA_integer_,
        n_half = length(rows),
        half_sample_index = paste(rows, collapse = ";"),
        selected_count = length(fit$support),
        mandatory_count = length(mandatory_core),
        mandatory_rank = half_partial$rank,
        mandatory_condition_number = half_partial$condition_number,
        eligible_candidate_count = sum(half_partial$eligible),
        ineligible_candidate_count = sum(!half_partial$eligible),
        selected_lambda = fit$lambda,
        selected_path_index = fit$selected_path_index,
        selected_path_locally_convex = fit$selected_path_locally_convex,
        converged = fit$converged,
        warning_count = length(fit$warnings),
        warning = if (length(fit$warnings)) {
          paste(fit$warnings, collapse = " | ")
        } else NA_character_,
        runtime_seconds = proc.time()[[3L]] - started,
        stringsAsFactors = FALSE
      )
    }
  }
  retention_frequency <- counts / (2L * complementary_pairs)
  retention_frequency[mandatory_core] <- 1
  optional_frequency <- rep(NA_real_, p)
  optional_frequency[candidate_set] <- counts[candidate_set] /
    (2L * complementary_pairs)
  optional_ranking <- eligible_candidate_set[order(
    -optional_frequency[eligible_candidate_set],
    -strength[eligible_candidate_set],
    feature_state$names[eligible_candidate_set],
    eligible_candidate_set,
    method = "radix"
  )]
  selected_extension <- optional_ranking[
    optional_frequency[optional_ranking] >= stability_threshold
  ]
  outside_candidate_set <- setdiff(seq_len(p), c(
    mandatory_core, candidate_set
  ))
  ranking <- c(
    mandatory_core, optional_ranking, ineligible_candidate_set,
    outside_candidate_set
  )
  core_set <- c(mandatory_core, selected_extension)
  pfer_applicable <- stability_threshold > 0.5
  threshold_regime <- if (pfer_applicable) {
    "above_half_stability_threshold_mb_pfer"
  } else {
    "exploratory_at_or_below_half_frequency_screen_no_mb_pfer"
  }
  stability_table <- data.frame(
    original_index = seq_len(p),
    feature = feature_state$names,
    stability_frequency = optional_frequency,
    cumulative_abs_standardized_coefficient = strength,
    selected_for_core = seq_len(p) %in% core_set,
    core_role = ifelse(
      seq_len(p) %in% mandatory_core, "mandatory",
      ifelse(
        seq_len(p) %in% selected_extension, "CPSS_extension",
        ifelse(
          seq_len(p) %in% eligible_candidate_set,
          "eligible_optional_not_selected",
          ifelse(
            seq_len(p) %in% ineligible_candidate_set,
            "ineligible_optional", "outside_candidate_universe"
          )
        )
      )
    ),
    stringsAsFactors = FALSE
  )
  role_order <- match(
    stability_table$core_role,
    c(
      "mandatory", "CPSS_extension", "eligible_optional_not_selected",
      "ineligible_optional", "outside_candidate_universe"
    )
  )
  stability_table <- stability_table[order(
    role_order,
    -stability_table$stability_frequency,
    -stability_table$cumulative_abs_standardized_coefficient,
    stability_table$feature,
    method = "radix"
  ), , drop = FALSE]
  rownames(stability_table) <- NULL
  output <- list(
    call = call,
    selector = paste0("CPSS-", toupper(selector)),
    base_selector = selector,
    core_set = as.integer(core_set),
    core_original_index = as.integer(core_set),
    core_feature_names = feature_state$names[core_set],
    mandatory_core = as.integer(mandatory_core),
    mandatory_core_original_index = as.integer(mandatory_core),
    mandatory_feature_names = feature_state$names[mandatory_core],
    candidate_set = as.integer(candidate_set),
    candidate_feature_names = feature_state$names[candidate_set],
    eligible_candidate_set = as.integer(eligible_candidate_set),
    eligible_candidate_feature_names =
      feature_state$names[eligible_candidate_set],
    ineligible_candidate_set = as.integer(ineligible_candidate_set),
    selected_extension = as.integer(selected_extension),
    selected_extension_original_index = as.integer(selected_extension),
    selected_extension_feature_names =
      feature_state$names[selected_extension],
    ranked_original_index = as.integer(ranking),
    selection_frequency = optional_frequency,
    stability_frequency = optional_frequency,
    retention_frequency = retention_frequency,
    cpss_selection_frequency = optional_frequency,
    stability_table = stability_table,
    base_fit_diagnostics = do.call(rbind, diagnostics),
    pair_diagnostics = do.call(rbind, diagnostics),
    feature_names = feature_state$names,
    feature_names_source = feature_state$source,
    dimensions = c(
      n = n, p = p, mandatory = length(mandatory_core),
      candidates = length(candidate_set),
      eligible_candidates = length(eligible_candidate_set),
      selected_extension = length(selected_extension),
      final_core = length(core_set)
    ),
    complementary_pairs = complementary_pairs,
    base_selection_size = base_selection_size,
    stability_threshold = stability_threshold,
    threshold_regime = threshold_regime,
    mb_pfer_applicable = pfer_applicable,
    pfer_upper_bound_mb = if (pfer_applicable) {
      base_selection_size^2 /
        ((2 * stability_threshold - 1) *
           length(eligible_candidate_set))
    } else NA_real_,
    pfer_bound_label = if (pfer_applicable) {
      "Meinshausen_Buehlmann_type_expression_requires_assumptions"
    } else {
      "not_applicable_at_or_below_one_half"
    },
    pfer_nonapplicability_reason = if (pfer_applicable) {
      NA_character_
    } else {
      "stability_threshold_not_greater_than_one_half"
    },
    path_lower_bound_semantics = if (identical(selector, "lasso")) {
      "glmnet_lambda_min_ratio"
    } else {
      "grpreg_lambda_min_ratio"
    },
    base_selector_engine = if (identical(selector, "lasso")) {
      "glmnet"
    } else {
      "grpreg_singleton_group_MCP"
    },
    mcp_local_convexity_diagnostic = if (identical(selector, "mcp")) {
      "not_available_from_grpreg"
    } else {
      "not_applicable"
    },
    conditional_on_mandatory_core = length(mandatory_core) > 0L,
    partialization = if (length(mandatory_core)) {
      "within_half_centered_FWL_economy_SVD"
    } else {
      "none"
    },
    rank_tolerance = rank_tolerance,
    pfer_family = "eligible_optional_candidates_only",
    no_top_k_fallback = TRUE,
    empty_core = !length(core_set),
    empty_extension = !length(selected_extension),
    selection_sample_independence_required = TRUE
  )
  class(output) <- c("cpss_core_selection", "list")
  output
}

.resolve_cpss_core <- function(X, core_set = NULL, selection = NULL) {
  if (inherits(core_set, "cpss_core_selection") && is.null(selection)) {
    selection <- core_set
    core_set <- NULL
  }
  if (!is.null(core_set) && !is.null(selection)) {
    stop("Supply either core_set or selection, not both.", call. = FALSE)
  }
  if (!is.null(selection)) {
    if (!inherits(selection, "cpss_core_selection") ||
        is.null(selection$core_set) || is.null(selection$feature_names)) {
      stop("selection must be returned by cpss_select_core().",
           call. = FALSE)
    }
    if (isTRUE(selection$empty_core) || !length(selection$core_set)) {
      stop("The strict CPSS selection produced an empty core.",
           call. = FALSE)
    }
    if (ncol(X) != length(selection$feature_names)) {
      stop("X does not contain the selection object's feature universe.",
           call. = FALSE)
    }
    if (identical(selection$feature_names_source, "column_names")) {
      current <- .cpss_feature_names(X)
      if (!setequal(current$names, selection$feature_names)) {
        stop("X feature names do not match the CPSS selection object.",
             call. = FALSE)
      }
      core_set <- match(selection$core_feature_names, current$names)
    } else {
      core_set <- selection$core_set
    }
  }
  if (is.null(core_set)) {
    stop("A nonempty core_set or CPSS selection object is required.",
         call. = FALSE)
  }
  core_set <- .validate_indices(core_set, ncol(X), "core_set")
  if (length(core_set) >= nrow(X)) {
    stop("length(core_set) must be strictly smaller than nrow(X).",
         call. = FALSE)
  }
  list(core_set = core_set, selection = selection)
}
