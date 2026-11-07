# Residuals of a fitted model

Response residuals (`y` minus the fitted mean) for every model class,
and for the single-equation classes Pearson residuals, the response
residual divided by the standard deviation of the fitted distribution
(its exact variance; the conditional one of a zero-truncated fit).

## Usage

``` r
# S3 method for class 'count_reg'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'gec'
residuals(object, type = c("response", "pearson"), ...)

# S3 method for class 'cpb_fe'
residuals(object, type = c("response", "pearson"), ...)
```

## Arguments

- object:

  A fitted model from this package.

- type:

  `"response"` (default) or, for `cpb`, `cpb_fe`, `gec`, `gec_fe` and
  `count_reg`, `"pearson"`.

- ...:

  Unused.

## Value

A numeric vector, in the rows of the data.
