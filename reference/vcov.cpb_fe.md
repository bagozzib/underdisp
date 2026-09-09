# Bootstrap covariance for a fixed-effects CPB fit

The covariance of the covariate coefficients from the stored
pairs/cluster bootstrap (`se = "bootstrap"` at fit time); `NULL` when no
bootstrap was run.

## Usage

``` r
# S3 method for class 'cpb_fe'
vcov(object, ...)
```

## Arguments

- object:

  A `"cpb_fe"` object.

- ...:

  Unused.

## Value

A covariance matrix, or `NULL`.
