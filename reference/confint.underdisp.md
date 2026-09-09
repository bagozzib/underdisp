# Confidence intervals for every model class

[`confint()`](https://rdrr.io/r/stats/confint.html) methods for the
whole family, one interval rule per inference type: bootstrap percentile
intervals for
[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md)
(coefficients; the dispersion parameter gets its profile-likelihood
interval) and
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md)
(coefficients and `delta`); normal-approximation intervals from the
bootstrap standard errors for the fixed-effects fits
([`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md),
[`gec_fe()`](https://bagozzib.github.io/underdisp/reference/gec_fe.md))
and the zero-inflated mixtures
([`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md),
[`zi_gec()`](https://bagozzib.github.io/underdisp/reference/zi_gec.md));
Wald intervals from the analytic covariance for
[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)
and
[`zi_count()`](https://bagozzib.github.io/underdisp/reference/zi_count.md)
(the dispersion parameter's interval is mapped from its estimation scale
to the natural scale); and, for the hurdles, the participation model's
Wald intervals (prefixed `part:`) stacked over the intensity model's
intervals (prefixed `int:`). A fit without inference errors
informatively.

## Usage

``` r
# S3 method for class 'gec'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'cpb_fe'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'count_reg'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'zi_count'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'hurdle_cpb'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'hurdle_gec'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'hurdle_count'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'zi_cpb'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'zi_gec'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A fitted model from this package.

- parm:

  Optional subset of parameter names.

- level:

  Confidence level (default 0.95).

- ...:

  Unused.

## Value

A two-column matrix of lower and upper bounds.

## See also

[`confint.cpb()`](https://bagozzib.github.io/underdisp/reference/confint.cpb.md)
