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
the second), `ratio` (`"IRR"` or `"OR"`), `estimate`, `lower`, `upper`,
and `method` (the interval source; `"none"` with `NA` bounds when no
covariance is available). Two-part models stack both equations. The
`IRR` rows are ratios of the rate parameter \\\exp(\beta)\\; for the CPB
and the underdispersed GEC the ratio of fitted means differs from it by
the support renormalization, and
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
reports the exact mean contrast.
