# Tidy a fitted model (broom method)

One row per coefficient with its estimate, standard error, z statistic
and p-value, for every model class of the package. The two-part classes
(`hurdle_*`, `zi_*`) have a `component` column (the broom convention for
multi-equation models) and the term prefixed by its component, so that
modelsummary lays the table out without a `shape` argument. A term of a
zero-inflation or participation equation that the data do not bound
(separation; see
[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md))
has no standard error.

## Usage

``` r
# S3 method for class 'cpb'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'cpb_fe'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'gec'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'count_reg'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'hurdle_cpb'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'hurdle_gec'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'hurdle_count'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'zi_cpb'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'zi_gec'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'zi_count'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A fitted model from this package.

- conf.int:

  If `TRUE`, add `conf.low`/`conf.high`, the intervals of
  [`confint()`](https://rdrr.io/r/stats/confint.html) (`NA` where the
  fit carries no inference).

- conf.level:

  Confidence level for the interval.

- ...:

  Unused.

## Value

A data frame with one row per coefficient.

## Details

The `scale` column names what each coefficient multiplies: `"log mean"`
for the mean-parameterized count families (`"log untruncated mean"` for
a zero-truncated fit or a hurdle's intensity), the binary equation's
link for a participation or zero-inflation row, and otherwise the
family's own parameter (`"log rate"` for the CPB and GEC, `"log lambda"`
for the rate-parameterized COM-Poisson, `"log centrality"` for the EWP,
and so on; see
[irr.count](https://bagozzib.github.io/underdisp/reference/irr.count.md)).
Coefficients on different scales are not comparable across the columns
of a table;
[`ame()`](https://bagozzib.github.io/underdisp/reference/ame.md) puts
every model's effect on the mean.
