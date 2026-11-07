# Predictions from a matched count-family fit

Predictions from a matched count-family fit

## Usage

``` r
# S3 method for class 'count_reg'
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

  A `"count_reg"` object.

- newdata:

  Optional data frame of covariate profiles.

- type:

  `"response"` (the mean E(Y); for a zero-truncated fit the conditional
  mean E(Y \| Y \> 0)), `"link"` (the linear predictor), or `"prob"`
  (P(Y = `at`); for a zero-truncated fit this is P(Y = `at` \| Y \> 0)).

- at:

  Count value for `type = "prob"`.

- offset:

  Optional offset (log scale) for `newdata`: a numeric vector or a
  column name in `newdata`.

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
