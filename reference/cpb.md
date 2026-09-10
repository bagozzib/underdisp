# Fit a Continuous Parameter Binomial (CPB) regression

Fits the underdispersed continuous parameter binomial model of King
(1989), in which the conditional variance is a fraction of the
conditional mean, \\\mathrm{Var}(Y\mid x) = \alpha\\\mathrm{E}(Y\mid
x)\\ with \\0 \< \alpha \< 1\\, and each observation has an endogenous
ceiling \\\lambda_i/(1-\alpha)\\. A zero-truncated variant (the default)
conditions on \\Y \ge 1\\, appropriate when underdispersion lives among
the positive counts of an otherwise zero-inflated outcome.

## Usage

``` r
cpb(
  formula,
  data,
  truncated = TRUE,
  se = c("none", "bootstrap"),
  B = 500,
  cluster = NULL,
  offset = NULL,
  weights = NULL,
  cores = 1L,
  alpha.start = 0.5,
  max.support = NULL,
  maxit = 20000,
  reltol = 1e-08
)
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame.

- truncated:

  Logical; if `TRUE` (default) fit the zero-truncated CPB (requires all
  `Y >= 1`); if `FALSE` fit the untruncated CPB on `Y >= 0`. **Note the
  sibling default differs:**
  [`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)
  defaults to `truncated = FALSE`, because its typical call fits a whole
  panel including zeros, while `cpb()`'s typical call fits the positive
  counts of a zero-inflated outcome. State `truncated` explicitly when
  moving a specification between the two.

- se:

  Inference method: `"none"` (default; fast, no standard errors) or
  `"bootstrap"` (cold-multistart pairs/cluster bootstrap; slower).

- B:

  Number of bootstrap resamples (default 500).

- cluster:

  Optional cluster identifier for cluster-robust inference: a column
  name in `data` or a vector aligned to its rows. When supplied, the
  bootstrap resamples whole clusters (a block bootstrap), giving
  cluster-robust standard errors and intervals; when `NULL` (default) it
  resamples observations, giving heteroskedasticity-robust errors. This
  is the appropriate route to robust inference here because the CPB's
  parameter-dependent support makes a Hessian-based sandwich unreliable.

- offset:

  Optional offset on the log-mean scale (an exposure): a numeric vector
  or the name of a column in `data`, making the model a rate model.

- weights:

  Optional frequency weights: a numeric vector or the name of a column
  in `data`. A weight of \\w_i\\ is equivalent to \\w_i\\ copies of row
  \\i\\; the bootstrap resamples rows together with their weights, and
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the weight
  total.

- cores:

  Number of worker processes for the bootstrap (default 1). With
  `cores > 1` the replicates run on a socket cluster whose random
  streams are seeded from the calling session, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  reproduces a parallel run (parallel and serial draws differ).

- alpha.start:

  Starting value for the dispersion parameter (default 0.5); the fit
  also multi-starts over a spread of alpha values.

- max.support:

  Guard on the maximum evaluated support: parameter values whose implied
  ceiling exceeds it are treated as infeasible (they lie in the `alpha`
  near 1 region where the support explodes). The default `NULL` sets it
  to `max(500, 10 * max(y))`. A fit whose largest ceiling reaches the
  guard is flagged with a warning and in `$support_binding`, because
  `alpha` is then bounded by the guard rather than by the data.

- maxit, reltol:

  Optimizer controls passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

An object of class `"cpb"`: a list with `coefficients`, `alpha`,
bootstrap `se.beta`/`vcov`/`ci.beta`, `loglik`, `loglik.null`,
`fitted.values`, `ceiling` (the observation-specific implied ceiling),
the model frame pieces, and the original `call`.

## Details

The mean is modelled log-linearly, \\\lambda_i = \exp(x_i'\beta)\\.
Because the support \\0, \ldots, \lfloor \lambda_i/(1-\alpha) \rfloor\\
depends on the parameters, the log-likelihood is discontinuous: whenever
a parameter move carries an observation's ceiling across an integer its
normalizing constant jumps, by about \\(1-\alpha)^{k}\\ at ceiling
\\k\\. The maximizer is found by a deterministic multistart, BFGS from
the Poisson solution at each value of `alpha.start` polished by
Nelder-Mead with feasibility-repaired restarts, on covariates scaled to
unit standard deviation (the coefficients are mapped back, so the fit
does not depend on the covariates' units), and the reported maximum is
the maximum of the fit's own profile in `alpha`: the profile is traced
by continuation on both sides of the estimate (the slopes re-maximized
at each step from the neighbouring solution) until it has dropped four
log-likelihood units, and a trace point above the multistart's value
restarts the fit from there. The same trace gives the profile interval
of
[`confint.cpb()`](https://bagozzib.github.io/underdisp/reference/confint.cpb.md).
Two parameterizations of the same design (say, treatment and sum
contrasts) can still settle on different teeth of this surface, with
log-likelihoods differing by the order of the jumps; unit intercepts
belong in the fixed-effects estimator
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md),
which solves each of them exactly (a pooled fit with many dummy columns
can stop short of it). The numerical Hessian is unreliable on such a
surface, so inference uses a cold-multistart bootstrap for the
coefficients (validated to nominal coverage) and a profile-likelihood
interval for \\\alpha\\ (see
[`confint.cpb()`](https://bagozzib.github.io/underdisp/reference/confint.cpb.md)).

Runtime: a fit with `se = "none"` takes a few seconds at 500 rows and
about fifteen at 2,000 on one core, the nine starts, the profile trace,
and the scan in `alpha` included. The bootstrap replicates are maximized
by the multistart without the profile trace (about a third of a second
each at 500 rows), so a replicate can sit on a different tooth from the
point estimate by the order of the jumps; `cores` runs them in parallel.

## References

King, G. (1989). Variance specification in event count models. *American
Journal of Political Science*, 33(3), 762-784.

## See also

[`ud_screen()`](https://bagozzib.github.io/underdisp/reference/ud_screen.md),
[`confint.cpb()`](https://bagozzib.github.io/underdisp/reference/confint.cpb.md),
[`predict.cpb()`](https://bagozzib.github.io/underdisp/reference/predict.cpb.md),
[`implied_ceiling()`](https://bagozzib.github.io/underdisp/reference/implied_ceiling.md),
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md)

## Examples

``` r
set.seed(1)
n <- 300; x <- rnorm(n)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
y <- rbinom(n, N, 0.5)                 # underdispersed (var/mean approx 0.5)
d <- data.frame(y = y, x = x)
fit <- cpb(y ~ x, data = d[d$y > 0, ], se = "none")
summary(fit)
#> 
#> Continuous Parameter Binomial regression (zero-truncated)
#> N = 298    inference: none 
#> 
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  1.58998         NA      NA       NA
#> x            0.50639         NA      NA       NA
#> 
#> alpha = 0.4862   (profile 95% CI: 0.457 to 0.591)
#> Implied ceiling lambda/(1-alpha): median 9.37   range 2.21 to 36.5 
#> logLik = -559.53    AIC = 1125.06 
#> LR vs ZT-Poisson (H0: alpha = 1): 51.09, p 4.4046e-13
```
