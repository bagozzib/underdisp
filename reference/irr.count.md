# Rate and odds ratios for the matched count families

Rate ratios (`exp` of the count/intensity coefficients, column `IRR`)
and, for the two-part models, odds ratios of the binary stage (column
`OR`), by equation. Element names match the CPB family (`intensity` for
the count component; `participation` for the hurdle stage, `zero` for
the inflation stage).

## Usage

``` r
# S3 method for class 'count_reg'
irr(object, level = 0.95, ...)

# S3 method for class 'gec'
irr(object, level = 0.95, ...)

# S3 method for class 'gec_fe'
irr(object, level = 0.95, ...)

# S3 method for class 'cpb_fe'
irr(object, level = 0.95, ...)

# S3 method for class 'hurdle_gec'
irr(object, level = 0.95, ...)

# S3 method for class 'zi_gec'
irr(object, level = 0.95, ...)

# S3 method for class 'hurdle_count'
irr(object, level = 0.95, ...)

# S3 method for class 'zi_count'
irr(object, level = 0.95, ...)
```

## Arguments

- object:

  A fitted model from this package.

- level:

  Confidence level.

- ...:

  Unused.

## Value

A `"ud_irr"` data frame – the package-wide ratio contract shared by
every [`irr()`](https://bagozzib.github.io/underdisp/reference/irr.md)
method: columns `term`, `equation` (`"count"`, `"participation"` for a
hurdle's binary stage, or `"inflation"` for a zero-inflated model's
structural-zero stage; the two point in opposite directions, an odds
ratio above one raising P(Y \> 0) in the first and P(structural zero) in
the second), `ratio` (`"IRR"`, `"OR"`, or the parameter's ratio named
below), `estimate`, `lower`, `upper`, and `method` (the interval source;
`"none"` with `NA` bounds when no covariance is available). Two-part
models stack both equations. The count rows are \\\exp(\beta)\\, the
ratio of the family's natural parameter. That is a ratio of means for
the Poisson, the negative binomial, the mean-parameterized COM-Poisson
(`"mpcmp"`), the hyper-Poisson, and the generalized Poisson with
\\\lambda \ge 0\\ (`"IRR"`). It is not for the others, whose rows name
what it is, with a message: the rate-parameterized COM-Poisson
(`"compois"`, `"lambda ratio"`; its mean moves roughly as the rate to
the power \\1/\nu\\, so that \\\exp(\beta)\\ is roughly the mean ratio
to the power \\\nu\\), the EWP (`"centrality ratio"`), the gamma count
(`"long-run rate ratio"`), the double Poisson (`"location ratio"`) and
the generalized Poisson below zero (`"mu ratio"`). The CPB and the
underdispersed GEC keep `"IRR"`, a ratio of their rate, which their mean
follows only approximately (the rate renormalized over a finite
support). For all of these,
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
and [`ame()`](https://bagozzib.github.io/underdisp/reference/ame.md)
report the exact mean contrast.
