# Hurdle continuous parameter binomial regression

Fits a participation (hurdle) model joined to a zero-truncated CPB for
the positive counts. Because the hurdle log-likelihood factorizes, the
two parts are fit separately: a logistic regression of participation on
all units, and a zero-truncated CPB
([`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md)) on
the positive counts. The result is the natural model for a bounded,
underdispersed participation process — most units at zero, participants
carrying a tight count.

## Usage

``` r
hurdle_cpb(
  formula,
  data,
  participation = NULL,
  fe = NULL,
  part_fe = NULL,
  link = c("logit", "probit", "cloglog"),
  offset = NULL,
  cluster = NULL,
  se = c("none", "bootstrap"),
  B = 500,
  weights = NULL,
  cores = 1L,
  ...
)
```

## Arguments

- formula:

  Intensity model formula (`y ~ x`).

- data:

  A data frame.

- participation:

  Optional one-sided formula (`~ z`) for the participation (hurdle)
  model; defaults to the intensity model's right-hand side.

- fe:

  Optional column name for unit fixed effects in the intensity model,
  absorbed by a concentrated likelihood
  ([`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)).
  This is how the within-unit underdispersion is recovered.

- part_fe:

  Optional column name for fixed effects in the participation model
  (added as factors). Units with no within-unit variation in
  participation are uninformative for a fixed-effects logit, and the
  dummy logit's coefficients carry the incidental-parameters bias of
  order 1/T (in short panels the slopes are inflated, roughly toward
  twice their value at T = 2); the split-panel jackknife
  (`bias_correct`) corrects only the intensity. The hurdle's likelihood
  separates its two parts, so for the logit link a conditional logit on
  the same data
  ([`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html))
  gives participation slopes free of that bias, and an analytical
  correction (the bife package, for the logit and the probit) reduces it
  to second order; the hurdle's predictions and quantities of interest
  use the dummy logit's unit effects. On King's support the dummy of a
  unit that never participates, at minus infinity, is not counted in the
  degrees of freedom (see
  [`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md),
  Degrees of freedom). In new data
  ([`predict.hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/predict.hurdle_cpb.md),
  [`score()`](https://bagozzib.github.io/underdisp/reference/score.md),
  [`cv_score()`](https://bagozzib.github.io/underdisp/reference/cv_score.md))
  such a unit's rows get the participation probability at that limit,
  zero, and those of a unit whose every count is positive
  (`part_one_units`) one; a unit the fit has not seen, or every row of
  `newdata` without the `part_fe` column, enters at the reference unit
  of the units the data identify. The participation equation is checked
  for separation by
  [`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md)'s
  rule (Separation; component `part_separation`). On the floor support
  these dummies are counted in the degrees of freedom, new rows of these
  units are evaluated at the dummies' fitted values, and an unseen unit
  enters at the mean of all the unit dummies.

- link:

  Link for the participation model: `"logit"` (default), `"probit"`, or
  `"cloglog"`. Positive participation coefficients raise the probability
  of a positive count.

- offset:

  Optional offset on the log-mean scale for the intensity (an exposure):
  a numeric vector or the name of a column in `data`.

- cluster:

  Optional cluster identifier (a column name in `data` or a vector). The
  participation equation's standard errors are then cluster-robust (the
  sandwich, with the G/(G - 1) factor), and with `se = "bootstrap"` the
  intensity's bootstrap resamples whole clusters; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- se:

  Inference for the intensity coefficients: `"none"` or `"bootstrap"`.

- B:

  Bootstrap replicates when `se = "bootstrap"`.

- weights:

  Optional frequency weights (a numeric vector or a column name),
  applied to both margins; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- cores:

  Worker processes for the intensity bootstrap; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- ...:

  Passed to
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md) or
  [`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md),
  including `support` (King's CPB, `"ceiling"`, by default; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md)).

## Value

An object of class `"hurdle_cpb"` with elements `participation` (a
`glm`), `intensity` (a `cpb` or `cpb_fe`), and bookkeeping.

## Examples

``` r
set.seed(1)
n <- 800; x <- rnorm(n); z <- rnorm(n)
y <- rhurdle_cpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5,
                 p = plogis(-0.2 + 0.8 * z))
fit <- hurdle_cpb(y ~ x, data = data.frame(y = y, x = x, z = z),
                  participation = ~ z)
fit
#> Hurdle Continuous Parameter Binomial
#> Call:  hurdle_cpb(formula = y ~ x, data = data.frame(y = y, x = x, z = z),     participation = ~z)
#> Units: 800 (392 participate, 49%)
#> 
#> Participation (logit link) -- positive coefficients raise P(Y > 0), i.e. participation:
#> (Intercept)           z 
#>     -0.0437      0.7652 
#> 
#> Intensity (zero-truncated CPB) coefficients:
#> (Intercept)           x 
#>      1.2862      0.5043 
#> Intensity alpha (shape): 0.5093 
```
