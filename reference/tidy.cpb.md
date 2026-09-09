# Tidy a CPB fit (broom method)

The two-part classes (`hurdle_*`, `zi_*`) return one row per coefficient
in each equation, with a `component` column (the broom convention for
multi-equation models) and the term prefixed by its component, so that
modelsummary lays the table out without a `shape` argument.

## Usage

``` r
# S3 method for class 'cpb'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'cpb_fe'
tidy(x, ...)

# S3 method for class 'hurdle_cpb'
tidy(x, ...)

# S3 method for class 'zi_cpb'
tidy(x, ...)

# S3 method for class 'gec'
tidy(x, ...)

# S3 method for class 'hurdle_gec'
tidy(x, ...)

# S3 method for class 'zi_gec'
tidy(x, ...)

# S3 method for class 'count_reg'
tidy(x, ...)

# S3 method for class 'zi_count'
tidy(x, ...)

# S3 method for class 'hurdle_count'
tidy(x, ...)
```

## Arguments

- x:

  A `"cpb"` object.

- conf.int:

  If `TRUE`, add `conf.low`/`conf.high` (requires a fit with bootstrap
  standard errors).

- conf.level:

  Confidence level for the interval.

- ...:

  Unused.

## Value

A data frame with one row per coefficient.
