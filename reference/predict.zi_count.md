# Predict from a zero-inflated Poisson/NB fit

Predict from a zero-inflated Poisson/NB fit

## Usage

``` r
# S3 method for class 'zi_count'
predict(
  object,
  newdata = NULL,
  type = c("response", "count", "zero", "intensity"),
  offset = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `"zi_count"` object.

- newdata:

  Optional data frame of covariate profiles. A fixed-effects unit an
  equation has not seen, or every row when its column is absent, enters
  that equation at the reference unit. A `zero_fe` unit with no positive
  count gets the structural-zero probability one and one whose every
  count is positive zero, and an `fe` unit with no positive count a
  count mean of zero (see
  [`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md),
  Degrees of freedom).

- type:

  `"response"` (marginal E(Y)), `"count"` (equivalently `"intensity"`;
  the count component's mean E(Y) of that component), or `"zero"`
  (structural-zero probability). `"intensity"` is accepted as a synonym
  for `"count"` for consistency with
  [`predict.zi_cpb()`](https://bagozzib.github.io/underdisp/reference/predict.zi_cpb.md).

- offset:

  Optional offset (log scale) of the count component for the `newdata`
  rows; by default the fit's own offset column is read from `newdata`
  (see
  [`predict.cpb()`](https://bagozzib.github.io/underdisp/reference/predict.cpb.md)).

- se.fit, interval, level:

  Standard errors (`se.fit = TRUE`) and a `level` confidence interval
  (`interval = "confidence"`) of the predicted quantity at the rows of
  `newdata`; see
  [`predict.cpb()`](https://bagozzib.github.io/underdisp/reference/predict.cpb.md).

- ...:

  Unused.

## Value

A numeric vector. With `interval = "confidence"`, a matrix with columns
`fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
vector or the matrix), `se.fit` and `method` (each row's interval
source). Without `se.fit`, the source is the attribute `method`.
