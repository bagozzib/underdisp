# Hurdle GEC (Katz-family) regression

Joins a participation (hurdle) logit to a zero-truncated GEC intensity
on the positive counts. Unlike
[`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md),
whose intensity is fixed to the underdispersed CPB, the GEC intensity
estimates its own dispersion direction.

## Usage

``` r
hurdle_gec(
  formula,
  data,
  participation = NULL,
  part_fe = NULL,
  link = c("logit", "probit", "cloglog"),
  offset = NULL,
  cluster = NULL,
  se = c("none", "bootstrap"),
  B = 500,
  weights = NULL,
  cores = 1L,
  max.support = 500
)
```

## Arguments

- formula:

  Intensity formula, `y ~ x`.

- data:

  A data frame.

- participation:

  One-sided formula for the participation logit; defaults to the
  intensity's right-hand side.

- part_fe:

  Optional column name for fixed effects in the participation model
  (added as factor dummies).

- link:

  Link for the participation model: `"logit"` (default), `"probit"`, or
  `"cloglog"`.

- offset:

  Optional offset (log scale) for the intensity: a numeric vector or a
  column name in `data`.

- cluster:

  Optional cluster (column name or vector) for the intensity's cluster
  bootstrap.

- se:

  Intensity inference: `"none"` (default) or `"bootstrap"`.

- B:

  Bootstrap resamples.

- weights:

  Optional frequency weights (a numeric vector or a column name),
  applied to both margins; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- cores:

  Worker processes for the intensity bootstrap; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- max.support:

  Guard on the maximum evaluated support.

## Value

An object of class `"hurdle_gec"`.

## Limitation

Unlike
[`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md),
`hurdle_gec()` has no `fe` argument for a concentrated fixed-effects
intensity (no concentrated zero-truncated GEC is implemented). For
within-unit intensities, enter the unit factor as dummies in `formula`
(feasible for moderate unit counts), or use
[`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md).

## See also

[`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md),
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md),
[`zi_gec()`](https://bagozzib.github.io/underdisp/reference/zi_gec.md)

## Examples

``` r
# \donttest{
set.seed(3); n <- 300; x <- rnorm(n); z <- rnorm(n)
y <- ifelse(rbinom(n, 1, plogis(0.3 + 0.8 * z)) == 1, rpois(n, exp(1 + 0.3 * x)) + 1L, 0L)
hurdle_gec(y ~ x, data.frame(y = y, x = x, z = z), participation = ~ z)
#> Hurdle Generalized Event Count (Katz family)
#> Call:  hurdle_gec(formula = y ~ x, data = data.frame(y = y, x = x, z = z),     participation = ~z)
#> Units: 300 (166 participate, 55%)
#> 
#> Participation (logit link) -- positive coefficients raise P(Y > 0), i.e. participation:
#> (Intercept)           z 
#>      0.2704      0.8885 
#> 
#> Intensity (zero-truncated GEC) coefficients:
#> (Intercept)           x 
#>      1.3504      0.2550 
#> Intensity dispersion delta: 0.683  [underdispersed (point estimate)]
# }
```
