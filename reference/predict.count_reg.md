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

- ...:

  Unused.

## Value

A numeric vector.
