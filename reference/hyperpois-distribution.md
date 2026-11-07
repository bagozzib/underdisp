# Hyper-Poisson distribution functions

Density, distribution, quantile, and random generation for the
hyper-Poisson distribution of Bardwell and Crow (1964) in the mean
parameterization of Saez-Castillo and Conde-Sanchez (2013): mean `mu`
and dispersion `gamma > 0`, with \\P(Y = y) = \lambda^y / ((\gamma)\_y\\
{}\_1F_1(1; \gamma; \lambda))\\, \\(\gamma)\_y = \Gamma(\gamma +
y)/\Gamma(\gamma)\\, and the rate \\\lambda\\ the unique solution of
\\E(Y) = \mu\\. `gamma < 1` gives underdispersion, `gamma = 1` the
Poisson, `gamma > 1` overdispersion. The underdispersion is bounded: as
`gamma -> 0` the distribution tends to 1 + Poisson(`mu` - 1) for `mu` \>
1 and to a Bernoulli(`mu`) below, so the variance is at least `mu - 1`
(`mu >= 1`) or `mu (1 - mu)` (`mu < 1`), and Var/Mean approaches one as
the mean grows (0.5 at a mean of 2, 0.9 at 10). The normalizer and the
moments are exact sums over the support.

## Usage

``` r
dhyperpois(x, mu, gamma, log = FALSE)

phyperpois(q, mu, gamma, lower.tail = TRUE, log.p = FALSE)

qhyperpois(p, mu, gamma, lower.tail = TRUE, log.p = FALSE)

rhyperpois(n, mu, gamma)
```

## Arguments

- x, q:

  Vector of quantiles (non-negative integers).

- mu:

  Mean (scalar or vector, recycled).

- gamma:

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

`dhyperpois` a density, `phyperpois` a CDF, `qhyperpois` a quantile,
`rhyperpois` a numeric vector of count draws.

## References

Bardwell, G. E. and Crow, E. L. (1964). A two-parameter family of
hyper-Poisson distributions. *Journal of the American Statistical
Association*, 59(305), 133-141. Sáez-Castillo, A. J. and Conde-Sánchez,
A. (2013). A hyper-Poisson regression model for overdispersed and
underdispersed count data. *Computational Statistics & Data Analysis*,
61, 148-157.

## See also

[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)

## Examples

``` r
dhyperpois(0:5, mu = 3, gamma = 0.5)
#> [1] 0.02867980 0.14422154 0.24174798 0.24313464 0.17466375 0.09759201
var(rhyperpois(2000, mu = 3, gamma = 0.5)) / 3
#> [1] 0.8686927
```
