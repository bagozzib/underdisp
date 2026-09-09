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
  ...
)
```

## Arguments

- object:

  A `"hurdle_count"` object.

- newdata:

  Optional data frame of covariate profiles.

- type:

  `"response"` (marginal E(Y), zeros included), `"participation"` (P(Y
  \> 0)), or `"intensity"` (E(Y \| Y \> 0)).

- offset:

  Optional offset (log scale) for the intensity when predicting on
  `newdata`; a numeric vector or a column name in `newdata`.

- ...:

  Unused.

## Value

A numeric vector.
