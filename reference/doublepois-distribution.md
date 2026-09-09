# Double Poisson distribution functions

Density, distribution, quantile, and random generation for Efron's
(1986) double Poisson distribution with location `mu` and dispersion
`theta`, normalized by the exact sum rather than Efron's closed-form
approximation. `theta > 1` gives underdispersion, `theta = 1` is the
Poisson, `theta < 1` overdispersion; the variance-to-mean ratio is
approximately `1/theta`, and the mean is close to (but not exactly)
`mu`. The exact mean and variance are what
`count_reg(family = "doublepois")` reports through
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html).

## Usage

``` r
ddoublepois(x, mu, theta, log = FALSE)

pdoublepois(q, mu, theta, lower.tail = TRUE, log.p = FALSE)

qdoublepois(p, mu, theta, lower.tail = TRUE, log.p = FALSE)

rdoublepois(n, mu, theta)
```

## Arguments

- x, q:

  Vector of quantiles (non-negative integers).

- mu:

  Location parameter (scalar or vector, recycled).

- theta:

  Dispersion parameter (scalar, positive).

- log, log.p:

  Return log probabilities.

- lower.tail:

  If `TRUE` (default), \\P(X \le x)\\.

- p:

  Vector of probabilities.

- n:

  Number of draws.

## Value

`ddoublepois` a density, `pdoublepois` a CDF, `qdoublepois` a quantile,
`rdoublepois` a numeric vector of count draws.

## References

Efron, B. (1986). Double exponential families and their use in
generalized linear regression. *Journal of the American Statistical
Association*, 81(395), 709-721.

## See also

[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)

## Examples

``` r
ddoublepois(0:5, mu = 3, theta = 2)
#> [1] 0.003568911 0.087311745 0.267005171 0.322575691 0.208081004 0.083404477
sum(ddoublepois(0:60, mu = 3, theta = 2))
#> [1] 1
```
