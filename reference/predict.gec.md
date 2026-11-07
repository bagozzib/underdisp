# Predictions from a generalized event count fit

Predictions from a generalized event count fit

## Usage

``` r
# S3 method for class 'gec'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "prob"),
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

  A `"gec"` object.

- newdata:

  Optional covariate profiles.

- type:

  `"response"` (the mean), `"link"` (the log mean), or `"prob"` (the
  probability of the count `at`).

- at:

  Count value(s) for `type = "prob"`.

- offset:

  Optional offset (log scale) for the `newdata` rows: a numeric vector
  or a column name in `newdata`; by default the fit's own offset column
  is read from `newdata` (see
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
