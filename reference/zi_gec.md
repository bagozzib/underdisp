# Zero-inflated GEC (Katz-family) regression

A structural-zero mixture with a GEC count component whose dispersion
`delta` is estimated freely. Estimated by robust direct maximum
likelihood, seeded from separate fits (a GEC on the positive counts for
the count parameters and a logit of the zero indicator for the
inflation), exactly as in
[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md);
this avoids the degenerate `pi -> 0` basin that coordinate-wise EM falls
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
  max.support = 500,
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
  entered as factor dummies (opt-in).

- se:

  Coefficient inference: `"none"` (default) or `"bootstrap"`.

- B:

  Bootstrap resamples.

- cluster:

  Optional cluster (column name or vector) for the bootstrap.

- max.support:

  Guard on the maximum evaluated support.

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
#>      1.0341      0.2551 
#> Count dispersion delta: 0.914  [underdispersed (point estimate)]
#> 
#> Zero-inflation (logit link) -- positive coefficients raise P(structural zero), i.e. lower the
#> chance of a positive count (the opposite direction from a hurdle participation model):
#> (Intercept)           z 
#>     -0.5634      0.8040 
#> Mean structural-zero probability: 0.371 
# }
```
