# GEC (Katz-family) regression with high-dimensional unit fixed effects

Fits [`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md)
with a full set of unit fixed effects, concentrating (profiling) out the
unit intercepts by a one-dimensional inner maximization per unit, so the
outer optimizer handles only the covariate coefficients and the free
dispersion parameter `delta`. This is the
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)
concentration generalized to the whole Katz family: the panel can be
under-, equi-, or overdispersed and the direction is estimated, not
presumed. For `delta < 1` the Katz support is finite, but the likelihood
is continuous across an integer ceiling (the entering support point's
mass grows from zero) and only kinks there, so each unit's objective is
unimodal in its intercept and is maximized by golden section on a
bracket around the unit's Poisson intercept, expanded while the maximum
sits at an edge; a maximum at a kink is tracked in the analytic gradient
of the concentrated likelihood, which the outer BFGS uses. The
feasibility floor is exact: every count must lie in
`0..ceiling(mu/(1-delta))`.

## Usage

``` r
gec_fe(
  formula,
  data,
  fe,
  se = c("none", "bootstrap"),
  B = 500,
  cluster = NULL,
  offset = NULL,
  weights = NULL,
  cores = 1L,
  max.support = NULL,
  inner_it = 30L,
  maxit = 3000L,
  reltol = 1e-07,
  bias_correct = c("none", "jackknife")
)
```

## Arguments

- formula:

  A model formula for the covariates only (no unit factor, no intercept;
  the fixed effects absorb it).

- data:

  A data frame.

- fe:

  Name of the column holding the unit identifier.

- se:

  Inference for the covariate coefficients: `"none"` (default) or
  `"bootstrap"`, a pairs/cluster bootstrap over units.

- B:

  Bootstrap resamples when `se = "bootstrap"`.

- cluster:

  Cluster for the bootstrap; `NULL` (default) resamples the
  fixed-effects units. See
  [`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md).

- offset:

  Optional offset on the log-mean scale (an exposure): a numeric vector
  or the name of a column in `data`.

- weights:

  Optional frequency weights (a numeric vector or a column name); see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- cores:

  Worker processes for the bootstrap; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- max.support:

  Guard on the maximum evaluated support.

- inner_it:

  Golden-section iterations for each unit's inner maximization.

- maxit, reltol:

  Outer optimizer controls.

- bias_correct:

  `"none"` (default) or `"jackknife"`, the split-panel jackknife
  correction for the 1/T incidental-parameters bias (see Details).

## Value

An object of class `c("gec_fe", "gec")`.

## Details

Limitation: unlike
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md),
`gec_fe()` has no `truncated` argument – a concentrated zero-truncated
GEC is not currently implemented. For a zero-truncated GEC with fixed
effects, enter the unit factor as dummies in
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md)'s
formula (feasible for moderate unit counts).

Short panels: `delta` and the fixed effects carry the
incidental-parameters bias of nonlinear fixed-effects estimation, of
order 1/T; treat `delta` cautiously below about T = 30. The covariate
coefficients are not materially affected. `bias_correct = "jackknife"`
removes the leading 1/T term by the split-panel jackknife exactly as in
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md):
refit on each unit's temporal halves and report
`2 * full - mean(halves)`, with the unit effects and fitted values
re-concentrated at the corrected parameters and `logLik`/`AIC` kept at
the maximum-likelihood fit (uncorrected estimates in `$uncorrected`).
The same time-homogeneity validity gate applies: on panels where the two
halves do not estimate a common parameter the correction is REFUSED with
a warning and the maximum-likelihood fit is returned (see
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md),
Details).

## See also

[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md),
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)

## Examples

``` r
# \donttest{
set.seed(5)
d <- do.call(rbind, lapply(1:30, function(i) {
  x <- rnorm(10); N <- pmax(round(exp(1 + rnorm(1, 0, 0.4) + 0.3 * x) / 0.5), 1)
  data.frame(unit = i, x = x, y = rbinom(10, N, 0.5))
}))
gec_fe(y ~ x, data = d, fe = "unit")   # delta ~ 0.5: underdispersed
#> GEC (Katz family) regression with 30 unit fixed effects (concentrated likelihood)
#> Coefficients:
#>      x 
#> 0.2897 
#> 
#> dispersion delta (Katz; Var/Mean on an unbounded support) = 0.488  [underdispersed],  n = 300
#> Note: delta is subject to incidental-parameters bias for short panels; see ?gec_fe.
# }
```
