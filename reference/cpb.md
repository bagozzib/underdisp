# Fit a Continuous Parameter Binomial (CPB) regression

Fits the underdispersed continuous parameter binomial model of King
(1989): a binomial count with a continuous number of trials,
\\\lambda_i/(1-\alpha)\\, and success probability \\1-\alpha\\, so that
the conditional variance is close to a fraction of the conditional mean,
\\\mathrm{Var}(Y\mid x) \approx \alpha\\\mathrm{E}(Y\mid x)\\ with \\0
\< \alpha \< 1\\, and each observation has an endogenous ceiling,
\\\lceil \lambda_i/(1-\alpha) \rceil\\. A zero-truncated variant (the
default) conditions on \\Y \ge 1\\, appropriate when underdispersion
lives among the positive counts of an otherwise zero-inflated outcome.

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
  reltol = 1e-08,
  support = c("ceiling", "floor")
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
  resamples observations, giving heteroskedasticity-robust errors (a
  Hessian-based sandwich is unreliable here; see Details). Robust errors
  protect inference about the best-fitting CPB's parameters. The mean
  effects computed from them are the data's only when the CPB is the
  data's process: on underdispersed counts from another process
  (binomial counts with varying trials, a COM-Poisson) the CPB's average
  marginal effects can be off by several percent at any sample size,
  while the robust Poisson's
  ([`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)
  with `se = "robust"`) are consistent whenever the mean is log-linear;
  the two fits'
  [`ame()`](https://bagozzib.github.io/underdisp/reference/ame.md) side
  by side is the check. Nor do robust errors protect the CPB's
  distributional claims (its ceiling, its predicted probabilities) when
  the CPB is wrong: King and Roberts (2015) read a gap between robust
  and model-based errors as a sign of misspecification to fix.
  [`ceiling_exceedance()`](https://bagozzib.github.io/underdisp/reference/ceiling_exceedance.md),
  [`rootogram()`](https://bagozzib.github.io/underdisp/reference/rootogram.md),
  [`pit_hist()`](https://bagozzib.github.io/underdisp/reference/pit_hist.md)
  and
  [`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md)
  check those claims.

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

  Guard on the evaluated support. On King's support the pmf is summed
  over the points within 36 log units of its mode, and a distribution
  with more such points than `max.support` is treated as infeasible (the
  guard caps the work for very large rates, binding at 500 where
  \\\alpha\lambda\\ exceeds about 880; it does not bind as `alpha`
  approaches 1). On the floor support parameter values whose implied
  ceiling exceeds it are infeasible (they lie in the `alpha` near 1
  region where that support explodes), and a fit whose largest ceiling
  reaches the guard is flagged with a warning and in `$support_binding`,
  because `alpha` is then bounded by the guard rather than by the data.
  The default `NULL` sets it to `max(500, 10 * max(y))`.

- maxit, reltol:

  Optimizer controls passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- support:

  The CPB's support. `"ceiling"` (the default) is King's (1989) CPB as
  Winkelmann, Signorino and King (1995) normalize it, on \\0, \ldots,
  \lceil \lambda_i/(1-\alpha) \rceil\\; `"floor"` is the package's
  definition before version 0.2.0, on \\0, \ldots, \lfloor
  \lambda_i/(1-\alpha) \rfloor\\, kept so that earlier results
  reproduce. The two coincide when the ceiling is an integer.

## Value

An object of class `"cpb"`: a list with `coefficients`, `alpha`,
bootstrap `se.beta`/`vcov`/`ci.beta`, `loglik`, `loglik.null`,
`fitted.values`, `ceiling` (the observation-specific implied ceiling),
the model frame pieces, and the original `call`.

## Details

The rate is modeled log-linearly, \\\lambda_i = \exp(x_i'\beta)\\. The
probabilities follow the Katz recursion \\P(Y = k+1)/P(Y = k) =
(\lambda_i - k(1-\alpha))/(\alpha(k+1))\\ on the support \\0, \ldots,
\lceil \lambda_i/(1-\alpha) \rceil\\, King's normalization as
Winkelmann, Signorino and King (1995) give it. The CPB is therefore the
generalized event count model of
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md) at
\\\delta = \alpha\\, and the two functions return the same fit wherever
the GEC's estimate of \\\delta\\ is below one, given the same
`truncated` (`cpb()` defaults to the zero-truncated model,
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md) to the
untruncated one). Where \\\lambda_i/(1-\alpha)\\ is an integer the
distribution is the binomial with that many trials, with mean
\\\lambda_i\\ and variance \\\alpha\lambda_i\\; elsewhere its mean is at
least \\\lambda_i\\. In King's notation `alpha` is \\\sigma^2\\, the
variance-to-mean ratio at an integer ceiling; Winkelmann, Signorino and
King write the ceiling \\\lambda_i/(1-\alpha)\\ as \\\alpha\\ and
\\1-\alpha\\ as \\p\\. Because the support moves with the parameters,
the log-likelihood is continuous with kinks wherever an observation's
ceiling crosses an integer, and it tends to the Poisson log-likelihood
as \\\alpha \to 1\\; a fit that reaches that limit (\\\hat\alpha \>
0.999\\) warns, and
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md)
gives it the boundary value of the test. The maximizer is found by a
deterministic multistart, Nelder-Mead from the Poisson solution at
`alpha.start` and seven further dispersions, the best polished by
Nelder-Mead with feasibility-repaired restarts and a final BFGS pass, on
covariates scaled to unit standard deviation (the coefficients are
mapped back, so the fit does not depend on the covariates' units), and
the reported maximum is the maximum of the fit's own profile in `alpha`:
the profile is traced by continuation on both sides of the estimate (the
slopes re-maximized at each step from the neighboring solution) until it
has dropped four log-likelihood units, and a trace point above the
multistart's value restarts the fit from there. The same trace gives the
profile interval of
[`confint.cpb()`](https://bagozzib.github.io/underdisp/reference/confint.cpb.md).
Unit intercepts belong in the fixed-effects estimator
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md),
which solves each of them exactly (a pooled fit with many dummy columns
can stop short of it).

The likelihood is not regular at the crossings. Just above an integer
the top support point carries probability proportional to the ceiling's
distance \\\epsilon\\ past it, so the information an observation carries
about its rate grows like \\1/\epsilon\\ (about \\0.28/\epsilon\\ for
\\\log\lambda\\ at \\\alpha = 0.5\\ and a ceiling of \\3 + \epsilon\\,
about \\1.8/\epsilon\\ at \\\alpha = 0.3\\; King 1989, footnote 7, notes
that asymptotic normality needs non-integer ceilings). With continuous
covariates the ceilings' fractional parts are spread and the expected
information per observation is unbounded, so standard errors from the
numerical Hessian are unreliable (at small `alpha` they fall well below
the sampling variability of the estimates), and inference uses a
cold-multistart bootstrap for the coefficients and a profile-likelihood
interval for \\\alpha\\. In simulations from King's CPB at n = 300, the
bootstrap's 95% interval for a slope covers 0.91 at alpha 0.3, 0.95 at
0.5 and 0.94 at 0.7 (0.92 at alpha 0.3 with n = 1,000); on the floor
support it covers 0.95. The interval for \\\alpha\\ is first-order by
default and calibrated by parametric bootstrap after
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md)
(see
[`alpha_confint()`](https://bagozzib.github.io/underdisp/reference/alpha_confint.md)).

`support = "floor"` gives the package's definition before version 0.2.0,
on \\0, \ldots, \lfloor \lambda_i/(1-\alpha) \rfloor\\, whose mean is at
most \\\lambda_i\\. Its log-likelihood is discontinuous: whenever a
parameter move carries an observation's ceiling across an integer its
normalizing constant jumps, by about \\(1-\alpha)^{k}\\ at ceiling
\\k\\. The fit then adds a scan of the breakpoints in `alpha` next to
the estimate; two parameterizations of the same design (say, treatment
and sum contrasts) can settle on different teeth of the surface, with
log-likelihoods differing by the order of the jumps; and ceilings above
`max.support` are infeasible, so `alpha` can be bounded by that guard
rather than by the data.

Runtime: a fit with `se = "none"` takes about a second at 500 rows and
four at 2,000 on one core, the multistart and the profile trace included
(about four and fifteen on the floor support). The bootstrap replicates
are maximized by the multistart without the profile trace (about 0.1
seconds each at 500 rows; on the floor support a replicate can therefore
sit on a different tooth from the point estimate, by the order of the
jumps); `cores` runs them in parallel.

## References

King, G. (1989). Variance specification in event count models: From
restrictive assumptions to a generalized estimator. *American Journal of
Political Science*, 33(3), 762-784. Winkelmann, R., Signorino, C. S. and
King, G. (1995). A correction for an underdispersed event count
probability distribution. *Political Analysis*, 5, 215-228. King, G. and
Roberts, M. E. (2015). How robust standard errors expose methodological
problems they do not fix, and what to do about it. *Political Analysis*,
23(2), 159-179.

## See also

[`ud_screen()`](https://bagozzib.github.io/underdisp/reference/ud_screen.md),
[`confint.cpb()`](https://bagozzib.github.io/underdisp/reference/confint.cpb.md),
[`predict.cpb()`](https://bagozzib.github.io/underdisp/reference/predict.cpb.md),
[`implied_ceiling()`](https://bagozzib.github.io/underdisp/reference/implied_ceiling.md),
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md)

## Examples

``` r
set.seed(1)
n <- 200; x <- rnorm(n)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
y <- rbinom(n, N, 0.5)                 # underdispersed (var/mean approx 0.5)
d <- data.frame(y = y, x = x)
fit <- cpb(y ~ x, data = d[d$y > 0, ], se = "none")
summary(fit)
#> 
#> Continuous Parameter Binomial regression (zero-truncated)
#> N = 199    inference: none 
#> 
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  1.61082         NA      NA       NA
#> x            0.51992         NA      NA       NA
#> 
#> alpha = 0.5147   (first-order profile 95% CI: 0.434 to 0.629)
#> Implied ceiling lambda/(1-alpha): median 10.08   range 3.26 to 35.96 
#> logLik = -377.49    AIC = 760.97 
#> LR vs ZT-Poisson (H0: alpha = 1): 33.10, p 4.3802e-09
#> Note: the p-value is asymptotic, which over-rejects in finite samples at this boundary; dispersion_test() gives the parametric-bootstrap p-value.
```
