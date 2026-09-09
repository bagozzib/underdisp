# Extensive/intensive first-difference decomposition for a hurdle- or ZI-CPB

The effect of moving `variable` `from` one value `to` another,
decomposed into the change in the participation (extensive) margin, the
intensity (intensive) margin, and the marginal expected count, holding
other covariates at their means. This is the quantity the
single-equation defaults cannot separate.

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

  Bootstrap replicates for percentile intervals; `0` (default) returns
  point estimates only. For `zi_cpb` the bootstrap is slow (each
  replicate refits the mixture).

- level:

  Confidence level.

- data:

  The original data frame; required when `B > 0`.

- stage:

  Which equation(s) the covariate moves in: `"both"` (default), the
  binary stage only (`"participation"`/`"zero"`), or the count stage
  only (`"intensity"`/`"count"`), holding the covariate at its reference
  in the other equation. All three component rows are always returned.

- cores:

  Worker processes for the bootstrap refits (default 1); see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- ...:

  Unused; unknown arguments error.

## Value

A `"ud_fd"` data frame (the package-wide first-difference contract):
columns `component`, `from`, `to`, `diff`, `lower`, `upper`, `method`,
one row per margin level.
