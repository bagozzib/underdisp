# First differences for the matched count families

The discrete-change effect on the expected count E(Y) as `variable`
moves from `from` to `to`, holding the other covariates at their sample
means. For `count_reg` this is a single number with a delta-method
interval; for the two-part models it is decomposed exactly into
extensive (participation / non-structural-zero) and intensive (count)
channels that sum to the total.

## Usage

``` r
# S3 method for class 'gec'
first_difference(object, variable, from, to, level = 0.95, ...)

# S3 method for class 'gec_fe'
first_difference(object, variable, from, to, level = 0.95, ...)

# S3 method for class 'cpb_fe'
first_difference(object, variable, from, to, level = 0.95, ...)

# S3 method for class 'hurdle_gec'
first_difference(
  object,
  variable,
  from,
  to,
  level = 0.95,
  stage = c("both", "participation", "intensity", "zero", "count"),
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
  ...
)

# S3 method for class 'count_reg'
first_difference(object, variable, from, to, level = 0.95, ...)

# S3 method for class 'hurdle_count'
first_difference(
  object,
  variable,
  from,
  to,
  level = 0.95,
  stage = c("both", "participation", "intensity", "zero", "count"),
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

  Confidence level for the delta-method intervals (all methods on this
  page that can compute one).

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
