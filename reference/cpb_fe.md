# CPB regression with high-dimensional unit fixed effects

Fits the CPB with a full set of unit fixed effects by concentrating
(profiling) out the unit intercepts. Each unit's intercept is solved by
a one-dimensional inner maximization, so the outer optimizer handles
only the covariate coefficients and the dispersion parameter. This
avoids an explicit unit-dummy design matrix and scales to thousands of
units, where dummy-based fitting exhausts memory. The concentrated
log-likelihood equals the full-dummy log-likelihood to the tolerance of
the inner search, which is solved exactly. Each unit's objective is a
saw-tooth in its intercept: whenever the intercept crosses
\\\log(k(1-\alpha)) - o\_{it}\\ for an integer \\k\\ the support of
observation \\t\\ gains the point \\k\\, its normalizing constant grows,
and the objective drops by about \\(1-\alpha)^k\\; between such
breakpoints it is smooth and unimodal. Its supremum is therefore
attained at an interior critical point of one tooth or as the left limit
at a breakpoint, and the inner search evaluates every breakpoint within
a window of the envelope's peak (located on a coarse grid around the
unit's Poisson intercept) and refines inside whole teeth by golden
section; every evaluation of the concentrated likelihood searches each
unit afresh, so the objective is a function of the parameters alone.
Solving the inner problem to its supremum makes the concentrated
likelihood smooth in `alpha` and continuous in the slopes, so the outer
optimizer (BFGS with the analytic envelope gradient from the Poisson
slopes and two dispersion starts, polished by Nelder-Mead) is reliable;
the solution is re-evaluated with every unit searched from cold before
it is returned.

## Usage

``` r
cpb_fe(
  formula,
  data,
  fe,
  truncated = FALSE,
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

  A model formula for the covariates only. Do not include the unit
  factor; supply it via `fe`. The intercept is absorbed by the fixed
  effects.

- data:

  A data frame.

- fe:

  Name of the column holding the unit identifier.

- truncated:

  Logical; fit the zero-truncated CPB if `TRUE`. **Note the sibling
  default differs:** `cpb_fe()` defaults to `FALSE` (whole-panel fits
  including zeros), while
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md)
  defaults to `TRUE` (positives-only fits). State `truncated` explicitly
  when moving a specification between the two.

- se:

  Inference for the covariate coefficients: `"none"` (default) or
  `"bootstrap"`, a pairs/cluster bootstrap.

- B:

  Bootstrap resamples when `se = "bootstrap"`.

- cluster:

  Cluster for the bootstrap. `NULL` (default) resamples the
  fixed-effects units themselves (unit-clustered inference, the natural
  panel default); a column name or vector resamples those clusters while
  preserving the unit fixed effects. Duplicated blocks are relabeled so
  the concentrated likelihood treats them as distinct units.

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

  Guard on the maximum evaluated support; `NULL` (default) sets
  `max(500, 10 * max(y))`. See
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- inner_it:

  Golden-section iterations refining each unit's inner maximization
  after the grid scans (see Details).

- maxit, reltol:

  Outer optimizer controls.

- bias_correct:

  `"none"` (default) or `"jackknife"`, the split-panel jackknife
  correction for the 1/T incidental-parameters bias (see Details).

## Value

An object of class `"cpb_fe"` with `coefficients`, `alpha`, `loglik`,
`fe` (the estimated unit intercepts), `fitted.values`, and the implied
per-observation `ceiling`.

## Details

Short panels: the dispersion parameter `alpha` and the fixed effects are
subject to the incidental-parameters bias of nonlinear fixed-effects
estimation. The bias in `alpha` is downward and of order 1/T, where T is
the per-unit number of observations. It is small for T at least about 30
and should be treated cautiously for short panels. The covariate
coefficients are not materially affected. `bias_correct = "jackknife"`
removes the leading 1/T term by the split-panel (half-panel) jackknife:
the model is refit on the first and second temporal halves of every
unit's series (rows are split in the order supplied, which should be
temporal order) and the corrected estimate is `2 * full - mean(halves)`.
In the package's fixed-effects bias Monte Carlo the alpha bias at T =
6/10/20/40 falls from -0.147/-0.096/-0.056/-0.033 to
-0.038/-0.022/-0.016/-0.009. With correction on, `coefficients` and
`alpha` are the corrected estimates (the maximum-likelihood values are
kept in `$uncorrected`), the unit effects and fitted values are
re-concentrated at the corrected parameters, and `logLik`/`AIC` continue
to refer to the maximum-likelihood fit.

Validity gate: the split-panel identity requires the two half-panels to
estimate the same pseudo-true parameter (time-homogeneity; Dhaene &
Jochmans 2015, Review of Economic Studies). On trending or
time-heterogeneous panels the correction is invalid, so the function
REFUSES it – returning the uncorrected maximum-likelihood fit with a
warning naming the failed check – whenever the corrected dispersion
leaves its feasible space or the two half-panel dispersion estimates
disagree beyond what the 1/T bias can explain. A refusal is diagnostic
information about the panel, not an error. The gate is deliberately
powered over sized: in the package's calibration Monte Carlo it refuses
about 9 percent of genuinely time-homogeneous panels (a conservative
nuisance; the returned fit is exactly the ordinary maximum-likelihood
estimate) while catching 98 percent of dispersion regime changes and 100
percent of smooth unmodeled trends – the cases where an uncaught
correction would be silently wrong.

One set of fixed effects is concentrated out. For two-way (e.g. unit and
time) fixed effects, add `factor(time)` to `formula` (entered as
dummies), or use
[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md),
which supports two-way fixed effects natively.

## See also

[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md),
[`ud_screen()`](https://bagozzib.github.io/underdisp/reference/ud_screen.md)

## Examples

``` r
# \donttest{
set.seed(1)
d <- do.call(rbind, lapply(1:60, function(i) {
  x <- rnorm(20); lam <- exp(rnorm(1, 0, 0.5) + 0.5 * x)
  N <- pmax(round(lam / 0.5), 1)
  data.frame(unit = i, x = x, y = rbinom(20, N, 0.5))
}))
fit <- cpb_fe(y ~ x, data = d, fe = "unit")
fit
#> CPB regression with 60 unit fixed effects (concentrated likelihood)
#> Coefficients:
#>      x 
#> 0.3635 
#> 
#> alpha (shape parameter): 0.5447   median implied bound: 2.35 
#> Note: alpha is subject to incidental-parameters bias for short panels; see ?cpb_fe.
# }
```
