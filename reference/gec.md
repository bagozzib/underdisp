# Generalized event count (Katz-family) regression

Fits King's generalized event count model, a Katz-family count
regression whose single dispersion parameter `delta` (the
variance-to-mean ratio) is estimated freely and spans underdispersion
(`delta < 1`, a finite-support member; the continuous parameter binomial
is this cell), equidispersion (`delta = 1`, Poisson), and overdispersion
(`delta > 1`, the negative binomial). Unlike
[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md), which
fixes the direction of dispersion to under, `gec()` lets the data
choose. The Katz recursion delivers exact first and second moments—the
Winkelmann–Signorino–King correction realized directly on an unbounded
support—and the likelihood is evaluated in C++. For `delta < 1` the
support is finite and the renormalized distribution's mean and variance
equal `exp(x'b)` and `delta * exp(x'b)` only when `exp(x'b)/(1 - delta)`
is an integer, so `exp(x'b)` is the rate parameter of the recursion, the
fitted mean is computed exactly from the pmf, and `delta` is the Katz
dispersion parameter (the variance-to-mean ratio on an unbounded
support). The optimizer works on covariates scaled to unit standard
deviation and maps the coefficients back, so the fit does not depend on
the covariates' units.

## Usage

``` r
gec(
  formula,
  data,
  truncated = FALSE,
  se = c("none", "bootstrap"),
  B = 500,
  cluster = NULL,
  offset = NULL,
  weights = NULL,
  cores = 1L,
  max.support = 500,
  maxit = 20000,
  reltol = 1e-08
)
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame.

- truncated:

  Logical; if `TRUE`, fit the zero-truncated GEC (all `Y >= 1`), the
  intensity model of
  [`hurdle_gec()`](https://bagozzib.github.io/underdisp/reference/hurdle_gec.md).

- se:

  Coefficient inference: `"none"` (default; fast) or `"bootstrap"`.

- B:

  Bootstrap resamples when `se = "bootstrap"`.

- cluster:

  Optional cluster identifier (a column name in `data` or a vector) for
  a cluster/block bootstrap; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- offset:

  Optional offset on the log-mean scale (an exposure): a numeric vector
  or the name of a column in `data`.

- weights:

  Optional frequency weights (a numeric vector or a column name); see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- cores:

  Worker processes for the bootstrap; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- max.support:

  Guard on the maximum evaluated support.

- maxit, reltol:

  Optimizer controls.

## Value

An object of class `"gec"` with `coefficients`, `delta` (the estimated
dispersion), `loglik`, bootstrap standard errors, and bookkeeping.

## See also

[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md),
[`compare_dispersion()`](https://bagozzib.github.io/underdisp/reference/compare_dispersion.md),
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md)

## Examples

``` r
set.seed(1); x <- rnorm(400)
y <- rpois(400, exp(1 + 0.5 * x))
gec(y ~ x, data = data.frame(y = y, x = x), se = "none")
#> Generalized event count (Katz family) regression
#> Call:  gec(formula = y ~ x, data = data.frame(y = y, x = x), se = "none")
#> (Intercept)           x 
#>      0.9760      0.4751 
#> 
#> dispersion delta (Katz; Var/Mean on an unbounded support) = 1.131  [overdispersed]
#> logLik = -764.34,  n = 400
```
