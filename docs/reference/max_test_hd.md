# Multiplier-bootstrap maximum test from debiased influence scores

The function requires only coordinatewise score variances. It does not
form or invert a `q` by `q` covariance matrix and therefore permits
`q > n`.

## Usage

``` r
max_test_hd(
  theta_tilde,
  target,
  psi_centered,
  alpha = 0.05,
  bootstrap_B = 999L,
  bootstrap_seed = NULL,
  multipliers = NULL,
  block_size = 256L
)
```

## Arguments

- theta_tilde:

  Debiased target estimates.

- target:

  Null target vector.

- psi_centered:

  Centered `n` by `q` influence-score matrix.

- alpha:

  Test level.

- bootstrap_B:

  Number of Gaussian multiplier draws.

- bootstrap_seed:

  Optional reproducible seed. The caller's RNG state is restored after
  multiplier generation.

- multipliers:

  Optional fixed `n` by `bootstrap_B` multiplier matrix.

- block_size:

  Number of draws processed together by the C++ kernel.

## Value

A list containing the observed maximum, bootstrap distribution, critical
value, p-value, separate null-moment calibrations, and the
backward-compatible `null_energy` alias for the null second moment.
