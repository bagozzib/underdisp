# Augment data with CPB fitted values (broom method)

Augment data with CPB fitted values (broom method)

## Usage

``` r
# S3 method for class 'cpb'
augment(x, data = NULL, ...)
```

## Arguments

- x:

  A `"cpb"` object.

- data:

  Optional: the data frame the model was fit to. Its rows the fit used
  (those without a missing model variable) are returned with the new
  columns appended.

- ...:

  Unused.

## Value

A data frame with fitted values, response residuals, and the implied
per-observation ceiling (appended to `data` when it is given).
