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

  Optional offset (log scale) for the count component when predicting on
  `newdata`; a numeric vector or a column name in `newdata`.

- ...:

  Unused.

## Value

A numeric vector.
