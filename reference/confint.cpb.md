# Confidence intervals for a CPB fit

Coefficient intervals use the cold-multistart bootstrap percentile
method (validated to nominal coverage). The interval for `alpha` is a
profile-likelihood interval: first-order by default (it covers about
0.88 to 0.95 in simulations from the CPB, with its misses on the upper
side; see
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md)
for the reason), and calibrated by parametric bootstrap when the fit
went through
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md).

## Usage

``` r
# S3 method for class 'cpb'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A `"cpb"` object fit with `se = "bootstrap"`.

- parm:

  Optional subset of parameters (coefficient names and/or `"alpha"`).

- level:

  Confidence level (default 0.95).

- ...:

  Unused.

## Value

A matrix of lower/upper bounds.
