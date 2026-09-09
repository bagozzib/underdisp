# Conway–Maxwell–Poisson distribution functions

Density, distribution, quantile, and random generation for the
COM-Poisson with rate `lambda` (or mean `mu`) and dispersion `nu`
(`nu > 1` underdispersed, `nu = 1` Poisson, `nu < 1` overdispersed).
Complements the estimators `count_reg(..., family = "compois")` (rate
parameterization) and `count_reg(..., family = "mpcmp")` (mean
parameterization).

## Usage

``` r
dcompois(x, lambda, nu, log = FALSE, mu = NULL)

pcompois(q, lambda, nu, lower.tail = TRUE, log.p = FALSE, mu = NULL)

qcompois(p, lambda, nu, lower.tail = TRUE, log.p = FALSE, mu = NULL)

rcompois(n, lambda, nu, mu = NULL)
```

## Arguments

- x, q:

  Vector of quantiles (non-negative integers).

- lambda:

  Rate parameter (scalar or vector, recycled). Give `mu` instead to
  specify the distribution by its mean.

- nu:

  Dispersion parameter (scalar).

- log, log.p:

  Return log probabilities.

- mu:

  Optional mean (scalar or vector, recycled); when supplied, the rate
  `lambda` solving \\\mathrm{E}(Y) = \mu\\ is found numerically (Huang's
  2017 mean parameterization, the one `count_reg(family = "mpcmp")`
  fits) and `lambda` is ignored.

- lower.tail:

  If `TRUE` (default), \\P(X \le x)\\.

- p:

  Vector of probabilities.

- n:

  Number of draws.

## Value

`dcompois` a density, `pcompois` a CDF, `qcompois` a quantile,
`rcompois` a numeric vector of count draws.

## Examples

``` r
dcompois(0:5, lambda = 3, nu = 1.5)
#> [1] 0.10062763 0.30188288 0.32019515 0.18486475 0.06932428 0.01860166
mean(rcompois(1000, lambda = 3, nu = 1.5))
#> [1] 1.939
```
