# Predictions from a CPB fit

Predictions from a CPB fit

## Usage

``` r
# S3 method for class 'cpb'
predict(
  object,
  newdata = NULL,
  type = c("response", "rate", "link", "ceiling", "prob"),
  at = NULL,
  offset = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `"cpb"` object.

- newdata:

  Optional data frame of new covariate profiles; if omitted, the fitted
  data are used.

- type:

  One of `"response"` (the exact mean of the fitted distribution; for a
  zero-truncated fit the conditional mean E(Y \| Y \>= 1)), `"rate"`
  (the CPB rate parameter lambda = exp(x'b), which equals the mean only
  when lambda/(1-alpha) is an integer), `"link"` (the linear predictor,
  log lambda), `"ceiling"` (the implied ceiling lambda/(1-alpha), the
  CPB's continuous number of trials; the largest count with positive
  probability is its ceiling on King's support and its floor on the
  floor support), or `"prob"` (the probability that `Y` equals `at`).

- at:

  For `type = "prob"`, the count value(s) `y` whose probability is
  returned (length 1, or one per row of the prediction data).

- offset:

  Optional offset (log scale) for the `newdata` rows: a numeric vector
  or a column name in `newdata`. By default a fit made with an offset
  column reads that column from `newdata`, and one made with an offset
  vector requires this argument. For `newdata = NULL` the fit's own
  offset is used.

- se.fit, interval, level:

  The uncertainty of the predicted quantity at the rows of `newdata`:
  `se.fit = TRUE` adds standard errors, and `interval = "confidence"` a
  `level` confidence interval. They follow the rule of
  [`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md):
  the fit's stored bootstrap draws (a percentile interval, and the
  draws' standard deviation), otherwise the delta method over the fit's
  covariance, otherwise none. They need `newdata`. A row at a
  fixed-effects unit's own effect, or at its own zero-equation dummy,
  has none, since the draws do not carry that parameter; a row of a unit
  the fit has not seen, or `newdata` without the unit column, enters at
  the reference unit, which they carry.
  [`ame()`](https://bagozzib.github.io/underdisp/reference/ame.md) gives
  average marginal effects.

- ...:

  Unused.

## Value

A numeric vector. With `interval = "confidence"`, a matrix with columns
`fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
vector or the matrix), `se.fit` and `method` (each row's interval
source). Without `se.fit`, the source is the attribute `method`.

## Examples

``` r
set.seed(1); x <- rnorm(300)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "none")
predict(fit, newdata = data.frame(x = 0), type = "ceiling")
#>          
#> 10.22242 
predict(fit, newdata = data.frame(x = c(-1, 1)), type = "prob", at = 5)
#> [1] 0.09402889 0.06443154
```
