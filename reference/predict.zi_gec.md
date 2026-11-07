# Predict from a zi_gec fit

Predict from a zi_gec fit

## Usage

``` r
# S3 method for class 'zi_gec'
predict(
  object,
  newdata = NULL,
  type = c("response", "zero", "intensity"),
  offset = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `"zi_gec"` object.

- newdata:

  Optional covariate profiles. A `zero_fe` unit the fit has not seen, or
  every row when that column is absent, enters at the zero equation's
  reference unit; a unit with no positive count gets the structural-zero
  probability one, and one whose every count is positive zero (see
  [`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md),
  Degrees of freedom).

- type:

  `"response"` (marginal E(Y)), `"zero"` (structural-zero probability),
  or `"intensity"` (the count mean E of the GEC component).

- offset:

  Optional offset (log scale) of the count component for the `newdata`
  rows; see
  [`predict.cpb()`](https://bagozzib.github.io/underdisp/reference/predict.cpb.md).

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
