# Zero-inflated GEC (Katz-family) regression

A structural-zero mixture with a GEC count component whose dispersion
`delta` is estimated freely. Estimated by direct maximum likelihood from
several starts, as in
[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md):
the single-equation GEC with a logit of the zero indicator for the
inflation, at four dispersions; the zero-inflated Poisson's
coefficients; and a low-inflation start. The best is polished, which
avoids the degenerate `pi -> 0` basin that coordinate-wise EM falls
into.

## Usage

``` r
zi_gec(
  formula,
  data,
  zero = NULL,
  zero_fe = NULL,
  se = c("none", "bootstrap"),
  B = 500,
  cluster = NULL,
  max.support = NULL,
  maxit = 200,
  tol = 1e-06,
  offset = NULL,
  weights = NULL,
  cores = 1L
)
```

## Arguments

- formula:

  Count formula, `y ~ x`.

- data:

  A data frame.

- zero:

  One-sided formula for the structural-zero logit; defaults to the count
  formula's right-hand side.

- zero_fe:

  Optional column name for fixed effects in the zero-inflation equation,
  entered as factor dummies (opt-in). The dummy of a unit with no
  positive count is not counted in the degrees of freedom (see
  [`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md),
  Degrees of freedom); in new data such a unit gets the structural-zero
  probability one, and one whose every count is positive zero.

- se:

  Coefficient inference: `"none"` (default) or `"bootstrap"`.

- B:

  Bootstrap resamples.

- cluster:

  Optional cluster (column name or vector) for the bootstrap.

- max.support:

  Guard on the evaluated support of the count component; `NULL`
  (default) sets `max(500, 10 * max(y))`. See
  [`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md).

- maxit, tol:

  Optimizer controls.

- offset:

  Optional exposure offset (log scale) for the count component: a
  numeric vector or the name of a column in `data`; the bootstrap
  carries it through every replicate.

- weights:

  Optional frequency weights (a numeric vector or a column name); see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- cores:

  Worker processes for the bootstrap; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

## Value

An object of class `"zi_gec"`.

## Separation

The fit checks its zero-inflation equation for separation as
[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md)
does (see its Separation section). It warns when the data put no useful
bound on a covariate coefficient, and stores the check in the component
`separation`.

## See also

[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md),
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md),
[`hurdle_gec()`](https://bagozzib.github.io/underdisp/reference/hurdle_gec.md)

## Examples

``` r
# \donttest{
set.seed(4); n <- 300; x <- rnorm(n); z <- rnorm(n)
y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.8 * z)) == 1, 0L, rpois(n, exp(1 + 0.3 * x)))
zi_gec(y ~ x, data.frame(y = y, x = x, z = z), zero = ~ z)
#> Zero-Inflated Generalized Event Count (Katz family)
#> Call:  zi_gec(formula = y ~ x, data = data.frame(y = y, x = x, z = z),     zero = ~z)
#> N: 300   logLik: -494.3
#> 
#> Count (GEC) coefficients:
#> (Intercept)           x 
#>      1.0339      0.2552 
#> Count dispersion delta: 0.914  [underdispersed (point estimate)]
#> 
#> Zero-inflation (logit link) -- positive coefficients raise P(structural zero), i.e. lower the
#> chance of a positive count (the opposite direction from a hurdle participation model):
#> (Intercept)           z 
#>     -0.5617      0.8045 
#> Mean structural-zero probability: 0.371 
# }
```
