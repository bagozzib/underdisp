# Summarize a CPB fit

Summarize a CPB fit

## Usage

``` r
# S3 method for class 'cpb'
summary(object, ...)
```

## Arguments

- object:

  A `"cpb"` object.

- ...:

  Unused.

## Value

An object of class `"summary.cpb"` with the coefficient table, the
dispersion parameter and its profile-likelihood interval (first-order,
or calibrated when the fit went through
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md)),
the implied ceiling, fit statistics, and the likelihood-ratio test
against a (zero-truncated) Poisson.
