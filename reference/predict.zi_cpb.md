# Predict from a zi_cpb fit

Predict from a zi_cpb fit

## Usage

``` r
# S3 method for class 'zi_cpb'
predict(
  object,
  newdata = NULL,
  type = c("response", "zero", "intensity", "rate"),
  offset = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A `"zi_cpb"` object.

- newdata:

  Data frame of covariate profiles. A `zero_fe` unit the fit has not
  seen, or every row when that column is absent, enters at the zero
  equation's reference unit, and a unit at a limit at that limit (see
  `zero_fe` in
  [`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md));
  for the count equation's `fe`, see `offset`.

- type:

  `"response"` for the marginal expected count E(Y), `"zero"` for the
  structural-zero probability, `"intensity"` for the mean of the CPB
  component (on the floor support, its rate lambda), or `"rate"` for the
  component's rate lambda.

- offset:

  Optional offset (log scale) of the count component for the `newdata`
  rows; see
  [`predict.cpb()`](https://bagozzib.github.io/underdisp/reference/predict.cpb.md).
  With `fe`, a row whose `fe` column names a unit of the fit enters at
  that unit's count effect, a row of an unseen unit (or every row, when
  `newdata` has no `fe` column) at the reference unit, the mean effect
  of the units the data identify. A unit left out for having no positive
  count (see
  [`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md),
  `fe`) predicts a mean of zero, and a structural-zero probability of
  `NA` where its own `zero_fe` dummy would be needed.

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
