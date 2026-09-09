# Generalized Poisson distribution functions

Density, distribution, quantile, and random generation for the
Consul-Jain generalized Poisson distribution in its constant-dispersion
form: mean `mu` and dispersion `lambda` in (-1, 1), with \\P(Y = y) =
\theta(\theta + \lambda y)^{y-1} e^{-\theta - \lambda y}/y!\\ and
\\\theta = \mu(1 - \lambda)\\, so that the variance-to-mean ratio is
\\1/(1-\lambda)^2\\. `lambda < 0` gives underdispersion, `lambda = 0`
the Poisson, `lambda > 0` overdispersion. For `lambda < 0` the support
is finite (the terms are positive only while \\\theta + \lambda y \>
0\\) and the pmf is renormalized on that support, so it is a proper
distribution and the exact mean and variance (which then differ from
`mu` and `mu/(1-lambda)^2` only by the negligible mass beyond the
support) are what `count_reg(family = "genpois")` reports through
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html).

## Usage

``` r
dgenpois(x, mu, lambda, log = FALSE)

pgenpois(q, mu, lambda, lower.tail = TRUE, log.p = FALSE)

qgenpois(p, mu, lambda, lower.tail = TRUE, log.p = FALSE)

rgenpois(n, mu, lambda)
```

## Arguments

- x, q:

  Vector of quantiles (non-negative integers).

- mu:

  Mean parameter (scalar or vector, recycled).

- lambda:

  Dispersion parameter in (-1, 1) (scalar).

- log, log.p:

  Return log probabilities.

- lower.tail:

  If `TRUE` (default), \\P(X \le x)\\.

- p:

  Vector of probabilities.

- n:

  Number of draws.

## Value

`dgenpois` a density, `pgenpois` a CDF, `qgenpois` a quantile,
`rgenpois` a numeric vector of count draws.

## References

Consul, P. C. and Jain, G. C. (1973). A generalization of the Poisson
distribution. *Technometrics*, 15(4), 791-799. Consul, P. C. and Famoye,
F. (2006). *Lagrangian Probability Distributions*. Birkhauser.

## See also

[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)

## Examples

``` r
dgenpois(0:5, mu = 3, lambda = -0.3)
#> [1] 0.02024191 0.10656252 0.23734318 0.29125435 0.21495599 0.09781863
var(rgenpois(2000, mu = 3, lambda = -0.3)) / 3   # about 1/(1.3)^2
#> [1] 0.6026076
```
