# Confidence intervals for a CPB fit

Coefficient intervals use the cold-multistart bootstrap percentile
method. In simulations from King's CPB (n = 300) the 95% interval for a
slope covers 0.91 at alpha 0.3, 0.95 at 0.5 and 0.94 at 0.7, and on the
floor support 0.95. The interval for `alpha` is a profile-likelihood
interval. It is first-order by default: on King's support it covers
close to its level in simulations from the CPB, and on the floor support
about 0.90 to 0.94, with most misses on the upper side (see
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md)).
It is calibrated by parametric bootstrap when the fit went through
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md).

## Usage

``` r
# S3 method for class 'cpb'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A `"cpb"` object. Coefficient intervals need `se = "bootstrap"` at fit
  time; the interval for `"alpha"` does not.

- parm:

  Optional subset of parameters (coefficient names and/or `"alpha"`).

- level:

  Confidence level (default 0.95).

- ...:

  Unused.

## Value

A matrix of lower/upper bounds.
