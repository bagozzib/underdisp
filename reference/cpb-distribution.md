# Continuous parameter binomial distribution functions

Density, distribution function, and quantile function for the continuous
parameter binomial (CPB), and its zero-truncated form. The CPB has a
hard ceiling at \\\lfloor \lambda/(1-\alpha) \rfloor\\; probability
above it is zero. Complements the simulator
[`rcpb()`](https://bagozzib.github.io/underdisp/reference/rcpb.md).

## Usage

``` r
dcpb(x, lambda, alpha, truncated = FALSE, log = FALSE)

pcpb(q, lambda, alpha, truncated = FALSE, lower.tail = TRUE, log.p = FALSE)

qcpb(p, lambda, alpha, truncated = FALSE, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- x, q:

  Vector of quantiles (non-negative integers).

- lambda:

  Mean parameter (scalar or vector, recycled).

- alpha:

  Shape parameter in (0, 1).

- truncated:

  If `TRUE`, use the zero-truncated CPB.

- log, log.p:

  If `TRUE`, probabilities are given as log.

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
