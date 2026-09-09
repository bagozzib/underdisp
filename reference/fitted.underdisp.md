# Fitted means for the two-part and fixed-effects classes

Completes
[`stats::fitted()`](https://rdrr.io/r/stats/fitted.values.html) across
the family, returning the same quantity as
`predict(object, type = "response")` on the estimation data: the mean
for the single-equation fits and the marginal expected count (zeros
included) for the two-part fits. `gec_fe` inherits `fitted.gec`.

## Usage

``` r
# S3 method for class 'cpb_fe'
fitted(object, ...)

# S3 method for class 'hurdle_cpb'
fitted(object, ...)

# S3 method for class 'hurdle_gec'
fitted(object, ...)

# S3 method for class 'hurdle_count'
fitted(object, ...)

# S3 method for class 'zi_cpb'
fitted(object, ...)

# S3 method for class 'zi_gec'
fitted(object, ...)

# S3 method for class 'zi_count'
fitted(object, ...)
```

## Arguments

- object:

  A fitted model from this package.

- ...:

  Unused.

## Value

A numeric vector, one element per estimation observation.
