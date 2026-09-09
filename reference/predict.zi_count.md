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
  ...
)
```

## Arguments

- object:

  A `"zi_count"` object.

- newdata:

  Optional data frame of covariate profiles.

- type:

  `"response"` (marginal E(Y)), `"count"` (equivalently `"intensity"`;
  the count component's mean E(Y) of that component), or `"zero"`
  (structural-zero probability). `"intensity"` is accepted as a synonym
  for `"count"` for consistency with
  [`predict.zi_cpb()`](https://bagozzib.github.io/underdisp/reference/predict.zi_cpb.md).

- offset:

  Optional offset (log scale) for the count component when predicting on
  `newdata`; a numeric vector or a column name in `newdata`.

- ...:

  Unused.

## Value

A numeric vector.
