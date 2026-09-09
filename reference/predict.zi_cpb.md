# Predict from a zi_cpb fit

Predict from a zi_cpb fit

## Usage

``` r
# S3 method for class 'zi_cpb'
predict(object, newdata = NULL, type = c("response", "zero", "intensity"), ...)
```

## Arguments

- object:

  A `"zi_cpb"` object.

- newdata:

  Data frame of covariate profiles.

- type:

  `"response"` for the marginal expected count E(Y), `"zero"` for the
  structural-zero probability, or `"intensity"` for the CPB mean.

- ...:

  Unused.

## Value

A numeric vector.
