# Glance at a fitted model (broom method)

One row of fit statistics for every model class of the package:
`logLik`, `AIC`, `BIC`, `df` and `nobs` (the weight total of a weighted
fit), the dispersion parameter as `dispersion` with its name
(`dispersion_parameter`) and the direction its estimate points
(`dispersion_direction`; the families' parameters point different ways,
so the CPB's `alpha` below 1 and the gamma-count's `alpha` above 1 are
both underdispersion), `coef_scale`, the scale of the count equation's
coefficients (see
[`tidy.cpb()`](https://bagozzib.github.io/underdisp/reference/tidy.cpb.md)),
and the class's own columns (`alpha` for the CPB family, `family` for
the matched count families, `n.units` for the fixed-effects fits).

## Usage

``` r
# S3 method for class 'cpb'
glance(x, ...)

# S3 method for class 'cpb_fe'
glance(x, ...)

# S3 method for class 'hurdle_cpb'
glance(x, ...)

# S3 method for class 'zi_cpb'
glance(x, ...)

# S3 method for class 'gec'
glance(x, ...)

# S3 method for class 'hurdle_gec'
glance(x, ...)

# S3 method for class 'zi_gec'
glance(x, ...)

# S3 method for class 'count_reg'
glance(x, ...)

# S3 method for class 'zi_count'
glance(x, ...)

# S3 method for class 'hurdle_count'
glance(x, ...)
```

## Arguments

- x:

  A fitted model from this package.

- ...:

  Unused.

## Value

A one-row data frame of fit statistics.
