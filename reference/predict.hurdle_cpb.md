# Predict from a hurdle-CPB fit

Predict from a hurdle-CPB fit

## Usage

``` r
# S3 method for class 'hurdle_cpb'
predict(
  object,
  newdata = NULL,
  type = c("response", "participation", "intensity"),
  ...
)
```

## Arguments

- object:

  A `"hurdle_cpb"` object.

- newdata:

  Data frame of covariate profiles (must contain both the intensity and
  participation covariates).

- type:

  `"response"` for the marginal expected count E(Y) (zeros included),
  `"participation"` for P(Y \> 0), or `"intensity"` for E(Y \| Y \> 0).

- ...:

  Unused.

## Value

A numeric vector.
