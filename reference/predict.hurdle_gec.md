# Predict from a hurdle-GEC fit

Predict from a hurdle-GEC fit

## Usage

``` r
# S3 method for class 'hurdle_gec'
predict(
  object,
  newdata = NULL,
  type = c("response", "participation", "intensity"),
  offset = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `"hurdle_gec"` object.

- newdata:

  Optional covariate profiles (intensity and participation). A `part_fe`
  unit the fit has not seen, or every row when that column is absent,
  enters at the reference unit, and a unit at a limit at that limit (see
  `part_fe` in
  [`hurdle_gec()`](https://bagozzib.github.io/underdisp/reference/hurdle_gec.md)).

- type:

  `"response"` (marginal E(Y)), `"participation"` (P(Y\>0)), or
  `"intensity"` (E(Y \| Y\>0)).

- offset:

  Optional offset (log scale) of the intensity for the `newdata` rows;
  see
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
