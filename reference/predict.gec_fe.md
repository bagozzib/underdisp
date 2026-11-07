# Predictions from a fixed-effects GEC fit

Predictions from a fixed-effects GEC fit

## Usage

``` r
# S3 method for class 'gec_fe'
predict(
  object,
  newdata = NULL,
  type = c("response", "link"),
  offset = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `"gec_fe"` object.

- newdata:

  Optional covariate profiles. A row whose `fe` column names a unit of
  the fit enters at that unit's estimated effect (as in
  [`score()`](https://bagozzib.github.io/underdisp/reference/score.md)
  and
  [`ceiling_exceedance()`](https://bagozzib.github.io/underdisp/reference/ceiling_exceedance.md));
  a row of an unseen unit, or every row when `newdata` has no `fe`
  column, enters at the reference unit, the mean fixed effect of the
  units whose effect the data identify (see Details of
  [`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)).

- type:

  `"response"` (the mean) or `"link"` (the log rate).

- offset:

  Optional offset (log scale) for the `newdata` rows; see
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
