# Continuous parameter binomial distribution functions

Density, distribution function, and quantile function for the continuous
parameter binomial (CPB), and its zero-truncated form. The CPB has a
hard ceiling: King's CPB (the default) puts no probability above
\\\lceil \lambda/(1-\alpha) \rceil\\, the floor support
(`support = "floor"`, the definition before version 0.2.0) none above
\\\lfloor \lambda/(1-\alpha) \rfloor\\. Complements the simulator
[`rcpb()`](https://bagozzib.github.io/underdisp/reference/rcpb.md).

## Usage

``` r
dcpb(
  x,
  lambda,
  alpha,
  truncated = FALSE,
  log = FALSE,
  support = c("ceiling", "floor")
)

pcpb(
  q,
  lambda,
  alpha,
  truncated = FALSE,
  lower.tail = TRUE,
  log.p = FALSE,
  support = c("ceiling", "floor")
)

qcpb(
  p,
  lambda,
  alpha,
  truncated = FALSE,
  lower.tail = TRUE,
  log.p = FALSE,
  support = c("ceiling", "floor")
)
```

## Arguments

- x, q:

  Vector of quantiles (non-negative integers).

- lambda:

  Rate parameter (scalar or vector, recycled). It is the mean where
  `lambda / (1 - alpha)` is an integer; otherwise the mean lies above it
  on King's support and below it on the floor support.

- alpha:

  Shape parameter in (0, 1).

- truncated:

  If `TRUE`, use the zero-truncated CPB. On the floor support, where the
  implied ceiling `lambda / (1 - alpha)` is below 1 the CPB puts all its
  mass on 0, and its zero-truncated form is the point mass at 1, the
  zero-truncated distribution for every ceiling in \[1, 2); King's
  support always contains 1.

- log, log.p:

  If `TRUE`, probabilities are given as log.

- support:

  The CPB's support, `"ceiling"` (King's CPB, the default) or `"floor"`;
  see [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- lower.tail:

  If `TRUE` (default), probabilities are \\P(X \le x)\\.

- p:

  Vector of probabilities.

## Value

`dcpb` a density, `pcpb` a distribution function, `qcpb` a quantile.

## See also

[`rcpb()`](https://bagozzib.github.io/underdisp/reference/rcpb.md)

## Examples

``` r
dcpb(0:5, lambda = 3, alpha = 0.5)
#> [1] 0.015625 0.093750 0.234375 0.312500 0.234375 0.093750
pcpb(3, lambda = 3, alpha = 0.5)
#> [1] 0.65625
qcpb(0.9, lambda = 3, alpha = 0.5)
#> [1] 5
```
