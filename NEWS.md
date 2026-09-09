# underdisp 0.1.1

## New families

* `count_reg()`, `hurdle_count()`, and `zi_count()` gain four underdispersed
  families through the same interface, each inheriting fixed effects,
  zero-truncation, hurdle and zero-inflated forms, offsets, weights, the
  sandwiches, proper scores, DHARMa simulation, and the QoI methods:
  `"mpcmp"` (the COM-Poisson in Huang's 2017 mean parameterization, with the
  rate solving E(Y) = mu found in C++), `"genpois"` (the Consul-Jain generalized
  Poisson with constant dispersion lambda, renormalized on its finite support
  when lambda < 0), `"gammacount"` (Winkelmann's 1995 renewal-process model
  with gamma waiting times), and `"doublepois"` (Efron's 1986 double Poisson
  with the exact normalizing constant). Every family's pmf was checked against
  an independent implementation (rmutil, gamlss.dist, VGAM) to machine
  precision, and the generalized-Poisson and mean-parameterized COM-Poisson
  log-likelihoods reproduce glmmTMB's `genpois` and `compois` fits exactly.
* `dgammacount()`, `ddoublepois()`, `dgenpois()` and their p/q/r companions;
  the COM-Poisson functions accept `mu =` for the mean parameterization.

## New diagnostics

* `dispersion_test()`: the likelihood-ratio test of a fitted family's
  dispersion parameter against its Poisson value on the fitted design (fixed
  effects included, the unit effects profiled out of the Poisson null), with
  two-sided or directional alternatives by the signed root and the Self-Liang
  boundary mixture where the Poisson value is on the boundary (CPB, negative
  binomial); and the Cameron-Trivedi auxiliary-regression test for a Poisson fit.
* `dispersion_profile()`: the empirical conditional variance-to-mean ratio by
  bins of the fitted mean against the ratio each fitted family implies, with a
  `plot()` method, so the mechanism behind an underdispersed outcome (a hard
  ceiling, regular event timing, a soft tail) can be read off the data.

## Weights, parallel bootstraps, and column scaling

* `weights =` (frequency weights) on every estimator; a weight of w is
  equivalent to w copies of the row in the likelihood, the information, the
  score products, and the bootstrap, and `nobs()` returns the weight total.
* `cores =` on every bootstrap, the parametric-bootstrap screen threshold,
  `cv_score()`, and the refit first differences: replicates run on a socket
  cluster seeded from the calling session.
* Every optimizer works on covariates scaled to unit standard deviation and
  maps the coefficients back, and the Nelder-Mead fits restart from their own
  solution until the improvement stops, so a fit no longer depends on the
  units of a covariate. The numerical Hessians and score matrices of the count
  families are differenced in the same scaled parameterization.

## The CPB mean

* `fitted()`, `predict(type = "response")`, `residuals()`, and
  `first_difference(quantity = "mean")` for `cpb`, `cpb_fe`, and the hurdle-CPB
  intensity now report the exact mean of the fitted distribution (a finite sum
  of the renormalized pmf, and the conditional mean E(Y | Y >= 1) for
  zero-truncated fits), as `gec()` and the matched families already did. The
  rate parameter exp(x'b) equals the mean only when lambda/(1-alpha) is an
  integer; it remains available as `$rate` and `predict(type = "rate")`.
  The same holds for the generalized event count: `delta` is the Katz
  dispersion parameter, which is the variance-to-mean ratio on an unbounded
  support, and `exp(x'b)` is the rate of the recursion.
* Zero-truncated `count_reg()` fits report E(Y | Y >= 1).

## Methods

* `confint()` for every model class (bootstrap percentile, normal-approximation
  from the bootstrap, or Wald from the analytic covariance, one rule per
  inference type); `vcov()` for the hurdle and zero-inflated classes;
  `implied_ceiling()` for `cpb_fe`; `irr()` and `first_difference()` on a `cpb`
  fit without a bootstrap return point estimates with `method = "none"`;
  `vcov()` returns `NULL` rather than an error when no inference was requested.
* `irr()` labels the binary stage `"participation"` (hurdle) or `"inflation"`
  (zero-inflated), which point in opposite directions.
* texreg `extract()` methods for `count_reg`, `hurdle_count`, `zi_count`; the
  methods for every class now attach to texreg's generic (in 0.1.0 the
  registration silently failed and tables fell back to broom), and the
  count families label their dispersion row by the parameter.
* `tidy()` for the two-part classes prefixes each term by its component, so
  `modelsummary()` lays the table out without a `shape` argument.
* `converged` is stored by every class and reported by `print()`/`summary()`
  when it is `FALSE`.

## Estimation and numerics

* The CPB profile-likelihood interval for `alpha` is traced by continuation
  with a feasibility-repaired start at each `alpha`, so the lower limit is the
  re-maximized profile's crossing point; the interval is one-sided only when
  the profile has not dropped to the cut by `alpha = 0.005`.
* `max.support` defaults to `max(500, 10 * max(y))`; a fit whose largest
  implied ceiling reaches the guard warns and records `$support_binding`; a
  fit whose every start is infeasible errors instead of returning the starting
  values; a CPB fit that sits below its own Poisson nest reports the negative
  likelihood-ratio statistic with no p-value.
* The concentrated fixed-effects fits solve each unit's intercept exactly.
  The unit objective is a saw-tooth in the intercept (the support floor moves
  as the rate crosses integer multiples of `1 - alpha`, and the objective drops
  at every such breakpoint), so its supremum sits at an interior critical point
  of a tooth or at the left limit of a breakpoint; the inner search enumerates
  the breakpoints within a window of the envelope's peak (as wide as the
  widest tooth in the unit) and refines inside whole teeth by golden section;
  every evaluation searches each unit afresh, so the concentrated likelihood is
  a function of the parameters alone, and the outer optimizer is BFGS with the analytic
  envelope gradient (the breakpoint a unit sits at is tracked in the gradient),
  polished by Nelder-Mead. The Katz family below `delta = 1` has a finite
  support too, but its likelihood only kinks at an integer ceiling (the
  entering point's mass grows from zero), so `gec_fe()` maximizes each unit's
  unimodal objective by golden section on an expanding bracket, with the
  feasibility floor exact (support `0..ceiling(mu/(1-delta))`). This makes the
  concentrated likelihood smooth in the dispersion parameter. The compiled
  code asserts the
  unit index against the design, and the fits use only the estimation rows, so
  rows dropped for missing values no longer misalign the unit blocks.
* `cpb()` documents the discontinuity of the CPB log-likelihood and the
  deterministic multistart that maximizes it (BFGS from the Poisson solution
  at each `alpha.start`, Nelder-Mead polish, feasibility-repaired restarts);
  the mean-parameterized COM-Poisson tabulates its log-factorials.
* Design matrices are rank-checked (an empty factor level or an exact linear
  combination is refused with the column named; a covariate constant within
  units is refused by the fixed-effects fits); a wrong-length vector `offset`,
  `weights`, or `cluster` is refused; `cluster` supplied without a bootstrap is
  announced; a cluster bootstrap with one cluster is refused and with fewer
  than thirty warned.
* `score()` raises `kmax` to the largest held-out count, and its predicted
  distributions carry the fit's offset (the in-sample log score equals
  `-logLik/n` for every class); the hurdle intensities are zero-truncated over
  the full support; `first_difference()` holds an offset at its mean and refuses
  a variable that enters an interaction; the Katz recursion runs at least to
  the observed count before its tail is dropped, so far-tail counts under an
  unbounded member get their small probability rather than an infeasibility
  sentinel; integer ceilings are feasible (floor with a tolerance) in the C++
  and the R pmf alike; the cluster-robust sandwich carries the G/(G-1) factor.
* The two-part models align both equations on the complete cases of every
  model variable; `zi_cpb()` records the user's call; `predict()` on `gec` and
  `gec_fe` uses the stored factor levels; the exported d/p/q/r functions
  return `numeric(0)` for zero-length input and refuse a vector dispersion
  parameter.

## Package

* `weights`, `cores`, and unknown-argument checks on the matched families
  (misspelled arguments are errors); `inst/CITATION`; a pkgdown site at
  https://bagozzib.github.io/underdisp/ and GitHub Actions running the full
  test suite (including the batteries skipped on CRAN) on five platforms.

# underdisp 0.1.0

First public release. `underdisp` provides diagnostics and a unified family of
estimators for underdispersed count data (conditional variance below the
conditional mean), the case the Poisson and negative binomial defaults cannot
represent. The likelihood core is implemented in C++.

## Screening

* `ud_screen()`: at-risk dispersion diagnostic that benchmarks a fitted count
  against a zero-truncated Poisson under unit and period fixed effects, with a
  marginal statistic, a zero-truncated-Poisson at-risk statistic, and an
  over-conditioning guard. The at-risk threshold is available in two modes,
  `ztp_threshold = "calibrated"` (a fast n-dependent calibrated cutoff) and
  `ztp_threshold = "bootstrap"` (a parametric bootstrap of the fitted
  zero-truncated-Poisson null on the data's own design, `ztp_boot_B`), and an
  optional native COM-Poisson comparator (`run_comp`).

## Estimators

* `cpb()`: continuous parameter binomial (CPB) regression and its zero-truncated
  variant, with an interpretable observation-specific bound (the implied
  ceiling). C++ multistart estimation.
* `cpb_fe()`: high-dimensional fixed-effects CPB via a concentrated likelihood
  that matches the full-dummy likelihood to machine precision.
* `hurdle_cpb()` and `zi_cpb()`: two-part (participation logit + zero-truncated
  CPB intensity) and true zero-inflated mixture forms, each with fixed effects
  available in both equations.
* `gec()`: King's generalized event count / Katz-family regression, a single
  model whose dispersion is estimated freely and spans underdispersion,
  equidispersion, and overdispersion, with `gec_fe()`, `hurdle_gec()`, and
  `zi_gec()` variants.
* `count_reg()`, `hurdle_count()`, `zi_count()`: matched Poisson, negative
  binomial, and COM-Poisson families sharing one interface, with fixed effects
  (including two-way), offsets, and truncated/hurdle/zero-inflated forms, for
  like-for-like comparison against the CPB and GEC families.

## Inference

* Coefficient inference by a cold-multistart pairs (or cluster) bootstrap with
  percentile intervals; dispersion inference by a profile-likelihood interval
  for the shape parameter, which respects the feasibility boundary where the
  bootstrap distribution is skewed (`alpha_confint()`).
* Robust and clustered standard errors throughout, via a block bootstrap rather
  than a Hessian sandwich, since the CPB's parameter-dependent support makes the
  Hessian unreliable. Fixed-effects bootstraps use normal-approximation
  intervals, because replicates share the incidental-parameters bias and
  percentile intervals would be off-center.

## Bias-corrected fixed-effects estimation

* `bias_correct = "jackknife"` on the concentrated fixed-effects estimators
  (`cpb_fe()`, `gec_fe()`, and the `fe =` intensity of `hurdle_cpb()`): the
  split-panel jackknife (Dhaene & Jochmans 2015) removes the leading 1/T
  incidental-parameters bias of the dispersion parameter, validated by the
  package's fixed-effects bias Monte Carlo (alpha bias at T = 6/10/20/40 falls
  from -.147/-.096/-.056/-.033 to -.038/-.022/-.016/-.009). The correction
  carries a VALIDITY GATE for the method's time-homogeneity requirement: the
  time-split disagreement is compared against exchangeable unit-split placebo
  noise, and on panels where the two halves do not estimate a common parameter
  the correction is refused with an informative warning and the
  maximum-likelihood fit is returned. A refusal is diagnostic information
  about the panel, not an error.
* For separation in dummy-heavy participation equations, the hurdle
  factorizes: fit the participation stage with a dedicated bias-reduction
  package (`logistf`, `brglm2`) alongside the zero-truncated intensity from
  this package.

## Comparison and testing

* `compare_dispersion()` and `compare_models()`: information-criteria and proper
  scoring-rule comparison across Poisson, negative binomial, COM-Poisson, and
  CPB (and their hurdle and zero-inflated forms), with type-safe rejection of
  non-package or mismatched objects.
* `zi_test()`: boundary-corrected likelihood-ratio test for zero-inflation
  (the plain model nested in the mixture at the boundary), for both the CPB and
  GEC families.

## Quantities of interest

* `predict()`, `implied_ceiling()`, `alpha_confint()`, `irr()`, and
  `first_difference()`. Every `first_difference()` method in the package
  returns one unified `"ud_fd"` data frame (`component`, `from`, `to`, `diff`,
  `lower`, `upper`, `method`): single-equation fits report the mean, and the
  two-part/mixture models report the exact extensive/intensive decomposition
  (binary-stage probability, count-stage mean, marginal mean), with a `stage`
  argument on every two-part method to move the covariate in a single equation
  and a per-row `method` label recording the uncertainty source (stored
  bootstrap for `cpb`, refit bootstrap via `B=` for the CPB two-parts,
  delta-method where an analytic covariance exists). Unknown arguments error
  rather than being silently ignored.

## Distributions and simulators

* Density, distribution, quantile, and random-generation functions for the
  package's own families: `dcpb`/`pcpb`/`qcpb`/`rcpb`,
  `dgec`/`pgec`/`qgec`/`rgec`, `dcompois`/`pcompois`/`qcompois`/`rcompois`,
  plus the mixture simulators `rhurdle_cpb()` and `rzicpb()`. The standard
  zero-truncated, hurdle, and zero-inflated Poisson/negative-binomial forms are
  deliberately internal: CRAN already serves those families (e.g. `VGAM`), and
  the exported surface is reserved for what the package adds.
* `simulate()` methods for every model class (`cpb`, `cpb_fe`, `gec`, `gec_fe`,
  `count_reg`, `hurdle_count`, `zi_count`, `hurdle_cpb`, `hurdle_gec`, `zi_cpb`,
  `zi_gec`), with `fitted()` completed across the two-part and fixed-effects
  classes, so any fit plugs into `DHARMa`'s simulated-residual diagnostics via
  `DHARMa::createDHARMa()`.

## Panel tools and calibration

* `mundlak()`: correlated-random-effects device that augments a formula with
  unit means, model-agnostic across the estimators.
* `score()` and `cv_score()` (log score and ranked probability score, in sample
  and by cross-validation), `rootogram()` (Tukey hanging), and `pit_hist()`
  (non-randomized PIT histogram).

## Numerical robustness

The estimators are hardened for the difficult likelihoods this package targets
(near-boundary dispersion, zero-heavy and dummy-heavy panels):

* Multistart estimation with per-start restart-polishing, and, for the
  zero-inflated CPB, an expectation-maximization cross-check run alongside the
  direct maximum-likelihood multistart so the better optimum is kept, avoiding
  the degenerate zero-mass-collapse basin.
* Loud, estimator-backed fallbacks that degrade gracefully rather than returning
  `NA`: a Fisher-scoring fallback for the zero-truncated-Poisson screen when the
  default IRLS fit errors or silently diverges, and a ridge-regularized Fisher
  step in the bootstrap null on singular information matrices. Each fallback
  warns and reports the method actually used.
* Fitted-rate sanity guards that key on the pathology (diverged rates) rather
  than its correlates (iteration counts), and floored mixture likelihoods so a
  transient zero-probability iterate cannot produce a non-finite objective.
* Fit-rank-exact degrees of freedom in the at-risk screen and the mixture
  comparisons, so zero-heavy designs are neither over- nor under-penalized.

## Data and integration

* Bundled `peacekeeping` dataset (UN peacekeeping contributions, a
  state-year roster with structural zeros).
* `broom` (`tidy`/`glance`/`augment`) methods for all model classes, and
  `modelsummary` / `texreg` table support.
