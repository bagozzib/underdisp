# Predictions from a fixed-effects CPB fit

Predictions from a fixed-effects CPB fit

## Usage

``` r
# S3 method for class 'cpb_fe'
predict(
  object,
  newdata = NULL,
  type = c("response", "rate", "link", "ceiling"),
  ...
)
```

## Arguments

- object:

  A `"cpb_fe"` object.

- newdata:

  Optional covariate profiles. With `newdata`, predictions use the
  average unit (the mean fixed effect), since the panel's own unit
  effects do not apply to new rows.

- type:

  `"response"` (the mean; for a zero-truncated fit the conditional mean
  E(Y \| Y \>= 1)), `"rate"` (the CPB rate lambda), `"link"` (the log
  rate), or `"ceiling"`.

- ...:

  Unused.

## Value

A numeric vector.
