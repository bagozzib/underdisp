# Confidence intervals for a CPB fit

Coefficient intervals use the cold-multistart bootstrap percentile
method (validated to nominal coverage); the interval for `alpha` uses
the profile-likelihood method, which is reliable except under strong
underdispersion, where `alpha` sits at the feasibility boundary and the
interval is one-sided.

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
