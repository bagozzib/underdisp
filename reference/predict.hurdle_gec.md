# Predict from a hurdle-GEC fit

Predict from a hurdle-GEC fit

## Usage

``` r
# S3 method for class 'hurdle_gec'
predict(
  object,
  newdata = NULL,
  type = c("response", "participation", "intensity"),
  ...
)
```

## Arguments

- object:

  A `"hurdle_gec"` object.

- newdata:

  Optional covariate profiles (intensity and participation).

- type:

  `"response"` (marginal E(Y)), `"participation"` (P(Y\>0)), or
  `"intensity"` (E(Y \| Y\>0)).

- ...:

  Unused.

## Value

A numeric vector.
