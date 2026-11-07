# Extensive/intensive first-difference decomposition for a hurdle- or ZI-CPB

The effect of moving `variable` `from` one value `to` another on the
participation (extensive) margin, the intensity (intensive) margin, and
the marginal expected count, holding other covariates at their means (or
at the values `at` gives them). The two margins are reported side by
side; their changes do not add up to the marginal change (see
[first_difference.count](https://bagozzib.github.io/underdisp/reference/first_difference.count.md)).
With `B = 0` the intervals follow the package's one rule (see
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)):
the stored bootstrap draws of a `zi_cpb` fit, and for a hurdle the delta
method over its two margins' covariances (the participation model's, and
the intensity's stored bootstrap draws). `B > 0` refits the whole model
on resampled data instead.

## Usage

``` r
# S3 method for class 'hurdle_cpb'
first_difference(
  object,
  variable,
  from,
  to,
  B = 0,
  level = 0.95,
  data = NULL,
  stage = c("both", "participation", "intensity", "zero", "count"),
  cores = 1L,
  at = NULL,
  ...
)

# S3 method for class 'zi_cpb'
first_difference(
  object,
  variable,
  from,
  to,
  B = 0,
  level = 0.95,
  data = NULL,
  stage = c("both", "participation", "intensity", "zero", "count"),
  cores = 1L,
  at = NULL,
  ...
)
```

## Arguments

- object:

  A `"hurdle_cpb"` or `"zi_cpb"` object.

- variable:

  Name of a model-matrix column to vary (in either margin).

- from, to:

  The two values to contrast.

- B:

  Refit-bootstrap replicates for percentile intervals; `0` (default)
  uses the fit's own inference instead (see Details). For `zi_cpb` the
  refit bootstrap is slow (each replicate refits the mixture).

- level:

  Confidence level.

- data:

  The original data frame; required when `B > 0`. Each replicate refits
  the model with the fit's own settings (fixed effects, link, offset,
  weights, support, and the like), resampling whole fixed-effects units
  (or the fit's clusters) when it has them; an offset or weights given
  to the fit as a vector must then be a column of `data` instead.

- stage:

  Which equation(s) the covariate moves in: `"both"` (default), the
  binary stage only (`"participation"`/`"zero"`), or the count stage
  only (`"intensity"`/`"count"`), holding the covariate at its reference
  in the other equation. All three component rows are always returned.

- cores:

  Worker processes for the bootstrap refits (default 1); see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- at:

  Optional named numeric vector fixing other model-matrix columns at
  given values instead of their means (see
  [`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)).

- ...:

  Unused; unknown arguments error.

## Value

A `"ud_fd"` data frame (the package-wide first-difference contract):
columns `component`, `from`, `to`, `diff`, `lower`, `upper`, `method`,
one row per margin level.
