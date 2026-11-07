# Exponentially weighted Poisson distribution functions

Density, distribution, quantile, and random generation for the
exponentially weighted Poisson of Ridout and Besbeas (2004), in its
three-parameter form (EWP3) and, with `beta2 = beta1`, its two-parameter
form (EWP2): a Poisson with centrality `lambda` reweighted by a
two-sided exponential in the distance from it, \$\$P(Y = y) \propto w_y
e^{-\lambda} \lambda^y / y!, \quad w_y = e^{-\beta_1 (\lambda - y)} \\
(y \le \lambda), \quad w_y = e^{-\beta_2 (y - \lambda)} \\ (y \>
\lambda).\$\$ `beta1` and `beta2` greater than zero thin the lower and
upper tails (underdispersion), values below zero thicken them, and
`beta1 = beta2 = 0` is the Poisson. `lambda` is not the mean; the
normalizer and the exact mean and variance are sums over the support (no
truncation of the normalizing sum at a fixed count), and the exact mean
is what `count_reg(family = "ewp3")` and `count_reg(family = "ewp2")`
report through [`predict()`](https://rdrr.io/r/stats/predict.html) and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html).

## Usage

``` r
dewp(x, lambda, beta1, beta2 = beta1, log = FALSE)

pewp(q, lambda, beta1, beta2 = beta1, lower.tail = TRUE, log.p = FALSE)

qewp(p, lambda, beta1, beta2 = beta1, lower.tail = TRUE, log.p = FALSE)

rewp(n, lambda, beta1, beta2 = beta1)
```

## Arguments

- x, q:

  Vector of quantiles (non-negative integers).

- lambda:

  Centrality parameter (scalar or vector, recycled).

- beta1, beta2:

  Lower- and upper-tail dispersion parameters (scalars); `beta2`
  defaults to `beta1`, the two-parameter EWP.

- log, log.p:

  Return log probabilities.

- lower.tail:

  If `TRUE` (default), \\P(X \le x)\\.

- p:

  Vector of probabilities.

- n:

  Number of draws.

## Value

`dewp` a density, `pewp` a CDF, `qewp` a quantile, `rewp` a numeric
vector of count draws.

## References

Ridout, M. S. and Besbeas, P. (2004). An empirical model for
underdispersed count data. *Statistical Modelling*, 4(1), 77-89.

## See also

[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)

## Examples

``` r
dewp(0:6, lambda = 3, beta1 = 1, beta2 = 1)
#> [1] 0.006081605 0.049594547 0.202217932 0.549685330 0.151663449 0.033476319
#> [7] 0.006157625
var(rewp(2000, lambda = 3, beta1 = 1, beta2 = 1))
#> [1] 0.8357929
```
