# Generalized event count (Katz) distribution functions

Density, distribution, quantile, and random generation for the King
(1989) generalized event count model with rate `lambda` and Katz
dispersion `delta` (`< 1` underdispersed, `= 1` Poisson, `> 1`
overdispersed). On an unbounded support (`delta >= 1`) the mean is
`lambda` and the variance-to-mean ratio is `delta` exactly; for
`delta < 1` the support is finite and the renormalized distribution's
mean and variance equal those values only when `lambda/(1 - delta)` is
an integer (the exact moments are finite sums of the pmf, as
[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)'s
siblings report them). Consistent with the estimator
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md).

## Usage

``` r
dgec(x, lambda, delta, max.support = 500, log = FALSE)

pgec(q, lambda, delta, max.support = 500, lower.tail = TRUE, log.p = FALSE)

qgec(p, lambda, delta, max.support = 500, lower.tail = TRUE, log.p = FALSE)

rgec(n, lambda, delta, max.support = 500)
```

## Arguments

- x, q:

  Vector of quantiles (non-negative integers).

- lambda:

  Rate parameter (scalar or vector, recycled).

- delta:

  Katz dispersion parameter (scalar).

- max.support:

  Guard on the evaluated support.

- log, log.p:

  Return log probabilities.

- lower.tail:

  If `TRUE` (default), \\P(X \le x)\\.

- p:

  Vector of probabilities.

- n:

  Number of draws.

## Value

`dgec` a density, `pgec` a CDF, `qgec` a quantile, `rgec` a numeric
vector of count draws.

## See also

[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md)

## Examples

``` r
dgec(0:5, lambda = 3, delta = 0.7)
#> [1] 0.02824752 0.12106082 0.23347444 0.26682793 0.20012095 0.10291935
var(rgec(2000, lambda = 3, delta = 0.7)) / mean(rgec(2000, lambda = 3, delta = 0.7))
#> [1] 0.708527
```
