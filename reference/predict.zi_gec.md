# Predict from a zi_gec fit

Predict from a zi_gec fit

## Usage

``` r
# S3 method for class 'zi_gec'
predict(object, newdata = NULL, type = c("response", "zero", "intensity"), ...)
```

## Arguments

- object:

  A `"zi_gec"` object.

- newdata:

  Optional covariate profiles.

- type:

  `"response"` (marginal E(Y)), `"zero"` (structural-zero probability),
  or `"intensity"` (the count mean E of the GEC component).

- ...:

  Unused.

## Value

A numeric vector.
