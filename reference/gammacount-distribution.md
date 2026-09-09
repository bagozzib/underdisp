# Gamma-count distribution functions

Density, distribution, quantile, and random generation for Winkelmann's
(1995) gamma-count distribution: the number of events in a unit interval
when the waiting times between events are independent Gamma variables
with shape `alpha` and rate `alpha * mu`, so that `mu` is the long-run
event rate. `alpha > 1` (waiting times more regular than exponential)
gives underdispersion, `alpha = 1` is the Poisson, `alpha < 1`
overdispersion; the variance-to-mean ratio is approximately `1/alpha`.
The exact mean and variance are finite sums of incomplete-gamma terms
and are what `count_reg(family = "gammacount")` reports through
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html).

## Usage

``` r
dgammacount(x, mu, alpha, log = FALSE)

pgammacount(q, mu, alpha, lower.tail = TRUE, log.p = FALSE)

qgammacount(p, mu, alpha, lower.tail = TRUE, log.p = FALSE)

rgammacount(n, mu, alpha)
```

## Arguments

- x, q:

  Vector of quantiles (non-negative integers).

- mu:

  Rate parameter (scalar or vector, recycled).

- alpha:

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

`dgammacount` a density, `pgammacount` a CDF, `qgammacount` a quantile,
`rgammacount` a numeric vector of count draws.

## References

Winkelmann, R. (1995). Duration dependence and dispersion in count-data
models. *Journal of Business & Economic Statistics*, 13(4), 467-474.

## See also

[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)

## Examples

``` r
dgammacount(0:5, mu = 3, alpha = 2)
#> [1] 0.01735127 0.13385262 0.29447576 0.29830012 0.17209622 0.06383205
var(rgammacount(2000, mu = 3, alpha = 2)) / 3
#> [1] 0.5282695
```
