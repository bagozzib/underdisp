# Predict from a hurdle-CPB fit

Predict from a hurdle-CPB fit

## Usage

``` r
# S3 method for class 'hurdle_cpb'
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

  A `"hurdle_cpb"` object.

- newdata:

  Data frame of covariate profiles (the intensity and participation
  covariates). A row of a fixed-effects unit the fit has not seen, or
  every row when the `fe` or `part_fe` column is absent, enters that
  equation at the reference unit, and a unit at a limit (no positive
  count, or no zero) at that limit (see `part_fe` in
  [`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md)).

- type:

  `"response"` for the marginal expected count E(Y) (zeros included),
  `"participation"` for P(Y \> 0), or `"intensity"` for E(Y \| Y \> 0).

- offset:

  Optional offset (log scale) of the intensity for the `newdata` rows: a
  numeric vector or a column name in `newdata`; by default the fit's own
  offset column is read from `newdata` (see
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
