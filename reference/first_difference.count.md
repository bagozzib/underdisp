# First differences for the matched count families

The discrete-change effect on the expected count E(Y) as `variable`
moves from `from` to `to`, holding the other covariates at their sample
means (or at the values `at` gives them). The interval follows the
package's one rule (see
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)):
percentile intervals from the fit's stored bootstrap draws for the GEC,
fixed-effects CPB and zero-inflated GEC fits, the delta method over the
full covariance, the dispersion parameter included, for the matched
count families and the hurdles' margins. For `count_reg` this is a
single number; for the two-part models it is reported on each margin:
the binary stage's probability (participation, or the structural-zero
probability), the count stage's mean, and the marginal mean, their
product, each at `from` and `to`. The first two are different quantities
and their changes do not add up to the marginal change, which is
\\\Delta(pm) = \Delta p\\ m_0 + p_0 \Delta m + \Delta p\\ \Delta m\\.

## Usage

``` r
# S3 method for class 'gec'
first_difference(object, variable, from, to, level = 0.95, at = NULL, ...)

# S3 method for class 'gec_fe'
first_difference(object, variable, from, to, level = 0.95, at = NULL, ...)

# S3 method for class 'cpb_fe'
first_difference(object, variable, from, to, level = 0.95, at = NULL, ...)

# S3 method for class 'hurdle_gec'
first_difference(
  object,
  variable,
  from,
  to,
  level = 0.95,
  stage = c("both", "participation", "intensity", "zero", "count"),
  at = NULL,
  ...
)

# S3 method for class 'zi_gec'
first_difference(
  object,
  variable,
  from,
  to,
  level = 0.95,
  stage = c("both", "participation", "intensity", "zero", "count"),
  at = NULL,
  ...
)

# S3 method for class 'count_reg'
first_difference(object, variable, from, to, level = 0.95, at = NULL, ...)

# S3 method for class 'hurdle_count'
first_difference(
  object,
  variable,
  from,
  to,
  level = 0.95,
  stage = c("both", "participation", "intensity", "zero", "count"),
  at = NULL,
  ...
)

# S3 method for class 'zi_count'
first_difference(
  object,
  variable,
  from,
  to,
  level = 0.95,
  stage = c("both", "participation", "intensity", "zero", "count"),
  at = NULL,
  ...
)
```

## Arguments

- object:

  A fitted `count_reg`, `gec`, `gec_fe`, `cpb_fe`, `hurdle_count`,
  `zi_count`, `hurdle_gec`, or `zi_gec` model.

- variable:

  Name of the covariate to change (must be in the model).

- from, to:

  The two values of `variable`.

- level:

  Confidence level of the intervals.

- at:

  Optional named numeric vector fixing other model-matrix columns at
  given values instead of their means (see
  [`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md));
  in a two-part model a column of both equations is fixed in both.

- ...:

  Unused; unknown arguments error.

- stage:

  For the two-part models, which equation(s) the change is applied to.
  `"both"` (default) moves the variable wherever it appears;
  `"intensity"`/`"count"` or `"participation"`/`"zero"` moves it in only
  that equation, holding it at the reference in the other. For a
  variable in only one equation all options coincide; for a variable in
  both, `"both"` gives the total effect and the single-stage options the
  partial effect through that margin. All component rows are always
  returned.

## Value

A `"ud_fd"` data frame – the package-wide first-difference contract
shared by every
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
method: columns `component`, `from`, `to`, `diff`, `lower`, `upper`,
`method`. Single-equation fits return one row (`component = "mean"`);
two-part fits return one row per margin level (binary-stage probability,
count-stage mean, marginal mean). `method` records the uncertainty
source per row (`"delta"`, a bootstrap label, or `"none"` with `NA`
bounds).
