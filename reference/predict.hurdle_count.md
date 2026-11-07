# Predict from a hurdle Poisson/NB fit

Predict from a hurdle Poisson/NB fit

## Usage

``` r
# S3 method for class 'hurdle_count'
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

  A `"hurdle_count"` object.

- newdata:

  Optional data frame of covariate profiles. A fixed-effects unit an
  equation has not seen, or every row when its column is absent, enters
  that equation at the reference unit, and a unit at a limit at that
  limit (see `fe,part_fe` in
  [`hurdle_count()`](https://bagozzib.github.io/underdisp/reference/hurdle_count.md)).

- type:

  `"response"` (marginal E(Y), zeros included), `"participation"` (P(Y
  \> 0)), or `"intensity"` (E(Y \| Y \> 0)).

- offset:

  Optional offset (log scale) for the intensity when predicting on
  `newdata`; a numeric vector or a column name in `newdata`.

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
