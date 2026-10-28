# underdisp 0.2.0

## The continuous parameter binomial follows King's definition

* Every CPB function now defaults to King's (1989) continuous parameter binomial
  as Winkelmann, Signorino, and King (1995) normalize it: the pmf on
  `0, ..., ceiling(lambda / (1 - alpha))`, which is the generalized event count
  of `gec()` at `delta = alpha < 1`. Earlier versions put the mass on
  `0, ..., floor(lambda / (1 - alpha))`, a distinct distribution that the
  documentation attributed to King. The two coincide when the ceiling is an
  integer. Under King's definition the realized mean sits at or above `lambda`
  (Winkelmann, Signorino, and King's equation 5), the log-likelihood is
  continuous in the parameters (the support point that enters as a ceiling
  crosses an integer carries mass that grows from zero), and a fit approaches
  the Poisson continuously as `alpha` approaches 1, so `max.support` no longer
  bounds `alpha`.
* The earlier definition remains available through `support = "floor"` in
  `cpb()`, `cpb_fe()`, `hurdle_cpb()`, `zi_cpb()`, `dcpb()`, `pcpb()`, `qcpb()`,
  `rcpb()`, `rhurdle_cpb()`, `rzicpb()`, `compare_dispersion()`, and
  `ud_screen()`, and reproduces 0.1.x results bit for bit (the zero-inflated fits
  can reach higher maxima; see below). A fit stores its `support` and every
  method reads it; an object saved by an earlier version is treated as a
  floor-support fit.
* On King's support the first-order profile interval for `alpha` covers
  between 0.935 and 0.956 in simulations from the CPB. On the floor support it
  covers 0.90 to 0.94, with most misses on the upper side. `calibrate_alpha()`
  remains available for either support. The documentation now presents it as
  the remedy for floor-support fits and as an optional refinement on King's.
* On King's support the pairs bootstrap's 95% percentile interval for a slope
  covers 0.92 at `alpha` = 0.3 and 0.94 at 0.5 and 0.7, in simulations at
  n = 300. Neither a parametric bootstrap (0.922) nor 999 replicates (0.932)
  brings it to the nominal level.
  The documentation reports these rates in place of "validated to nominal
  coverage", which held on the floor support (0.95).
* In the fixed-effects bias simulations (150 units, `alpha` = 0.5) the bias of
  `cpb_fe()`'s `alpha` under King's definition is -0.122 at T = 6 and -0.018 at
  T = 40 (floor support: -0.153 and -0.033), and `bias_correct = "jackknife"`
  removes it to within 0.009 (on the floor support it leaves 23 to 30 percent).
  The jackknife's validity gate refuses 9 to 13 percent of time-homogeneous
  panels (floor support: about 7 percent).
* A King's-support fit whose `alpha` reaches the equidispersed limit (above
  0.999) warns that the data show no underdispersion, in place of the support
  guard's warning, and `dispersion_test()` and `summary()` report the boundary
  value 0 there.
* At the other end, a fit whose `alpha` falls below 0.001 prints a note: as
  `alpha` approaches 0 the pmf at each rate puts its mass on at most two
  adjacent counts, which a series taking only two neighbouring values fits
  exactly, and the interval and tests for `alpha` are at that boundary. The
  maximum is legitimate, so this is a note, not a warning.
* `predict()` and `fitted()` of a King's-support `zi_cpb()` fit use the exact
  mean of the count component, as `zi_gec()` does.

## New families

* `count_reg()`, `hurdle_count()`, and `zi_count()` gain the hyper-Poisson
  (`"hyperpois"`, Saez-Castillo and Conde-Sanchez 2013, in its mean
  parameterization) and the exponentially weighted Poisson of Ridout and Besbeas
  (2004) with one or two tail weights (`"ewp2"`, `"ewp3"`), with
  `dhyperpois()`/`dewp()` and their p/q/r companions. `dispersion_test()` tests
  the two weights of `"ewp3"` jointly.

## Optimization

* `zi_cpb()`, `zi_gec()`, and `zi_count()` start from points that span the
  structural share of zeros, where the mixture likelihood's optima differ: the
  zero-inflated Poisson's coefficients and a low-inflation start (the
  single-equation fit with a constant inflation probability of 0.02) beside the
  separate-stage seeds, followed by a BFGS pass. On the pooled conflicts panel of
  the companion paper the earlier starts stopped 25 log-likelihood points short.
* `zi_cpb(fe = , method = "em")` follows its classification EM, which maximizes
  the classification likelihood rather than the mixture likelihood it reports,
  with block coordinate ascent on the observed-data likelihood.
* `gec()` polishes its Nelder-Mead solution by BFGS.
* `count_reg()` evaluates its likelihood and its covariance over the distinct
  rows of the data (the same count, covariates and offset, within a cluster for
  a clustered covariance), each weighted by its count, and so do the intensity
  of `hurdle_count()` and the starts of `zi_count()`. The results are those of
  the full data up to the order of summation; a design with repeated rows runs
  faster, by a factor of 5 to 20 on a three-level factor in the package's
  checks, and an intercept-only zero-truncated EWP on 5,414 rows of seven
  distinct counts fits in about a second instead of eleven minutes.
* `zi_cpb()`, `zi_gec()`, and `zi_count()` check their zero-inflation equation
  for separation. They warn when the data put no useful bound on a covariate
  coefficient.
  - The check uses the profile likelihood of each coefficient, with the count
    component held at its estimate. A coefficient is flagged when its profile
    falls by less than `qchisq(0.95, 1) / 2` over 100 standard deviations of its
    covariate.
  - The result is stored in the component `separation`, and the print methods
    note it. They also note when the check finds a higher point of the
    likelihood than the fit.
  - When the check finds a point at least 0.01 higher, the fit restarts its
    final polish there, keeps the result only if it gains, and checks again, up
    to three times.
  - A fit near the no-inflation boundary (data with few or no structural zeros)
    leaves every zero-equation coefficient unbounded for that reason alone. It
    is not reported as separated; the print methods note the boundary instead.
  - In a mixture, separation needs no covariate that splits the zeros from the
    positives, because the count component's probability of a zero covers the
    zeros on the wrong side. A hurdle on the same covariates can therefore have
    a finite maximum where the mixture has none.
  - Earlier versions muffled the warnings of the zero equation's internal logit
    fits and reported nothing.

## The ceiling as a quantity and as a test

* `implied_ceiling()` gains `max_count`, the largest count with positive
  probability (the integer ceiling of `lambda / (1 - alpha)` on King's
  support, its floor on the floor support); `ceiling` stays the continuous
  number of trials.
* New `ceiling_exceedance(fit, newdata)` counts held-out observations above the
  fitted bound, for `cpb()`, `cpb_fe()`, `hurdle_cpb()` and `zi_cpb()` fits.
* `compare_dispersion()` and `ud_screen()` no longer report a ceiling-exceedance
  share: computed on the data the CPB was fit to, it was zero for every
  feasible fit.
* `predict.zi_cpb()` gains `type = "rate"`.

## Quantities of interest and their intervals

* Every `first_difference()` method takes its interval by one rule. A fit with
  stored bootstrap draws gives the percentile interval of the quantity computed
  at every draw, jointly over the coefficients, the dispersion parameter and,
  for a fixed-effects fit, the reference unit's effect; a fixed-effects fit's
  draws are first recentered on its reported estimate, since its replicates
  share the incidental-parameters bias (this also centres a jackknife-corrected
  fit's intervals on the corrected estimate). Other fits use the delta method
  over the full covariance, the dispersion parameter included, and a hurdle the
  independent covariances of its two margins. The intervals of `gec()`,
  `gec_fe()`, `cpb_fe()` and `count_reg()` used to hold the dispersion
  parameter fixed, and the zero and marginal rows of `zi_count()` and the
  decompositions of `hurdle_cpb()`, `zi_cpb()`, `hurdle_gec()` and `zi_gec()`
  without `B` had none.
* `first_difference()` gains `at =`, a named vector that fixes other
  model-matrix columns (a factor's dummies, say) instead of holding them at
  their means.
* `predict()` gains `se.fit`, `interval = "confidence"` and `level` for every
  class, by the same rule: the predicted quantity at each row of `newdata` is
  computed at every stored draw, or by the delta method over the fit's
  covariance. A row at a fixed-effects unit's own effect, or at its own
  zero-equation dummy, has no interval, since the draws do not carry that
  parameter; a row of an unseen unit, or `newdata` without the unit column,
  enters at the reference unit, which they do. These arguments used to be
  ignored.
* New `ame()`: average marginal effects, the average over the rows of `newdata`
  (weighted by `weights =`, such as a tabulated fit's frequency weights) of the
  derivative of a predicted quantity in a numeric variable, or of its contrast
  between two values of any variable, with the interval of the same rule. The
  variable enters through the model's formula, so interactions and
  transformations of it are handled, which `first_difference()` refuses. For a
  concentrated fixed-effects fit, rows at their own unit effects (the sample's
  average partial effect) have no interval from the stored draws; rows at the
  reference unit (a typical, geometric-mean unit) do. `ame(B =)` gives the
  average partial effect its interval: it refits the model on resampled units
  (or the fit's clusters, or rows), averages each replicate over its own rows at
  its own unit effects, and takes the percentile interval, an interval for the
  average partial effect in the population of units. The derivative of a
  floor-support fit, whose mean jumps at integer ceilings, is refused.
* `predict(newdata)` on a fixed-effects fit (`cpb_fe()`, `gec_fe()`,
  `zi_cpb(fe =)`) enters a unit the fit estimated at its own effect, as
  `score()` and `ceiling_exceedance()` do, and an unseen unit, or every row when
  `newdata` has no unit column, at the reference unit. Every row used to enter
  at the reference unit.
* `irr()` takes its intervals as `confint()` does (for `gec()`, the exponential
  of the percentile interval). A term the data do not bound has `NA` ratios,
  with a message, and the fixed-effects dummies are left out of the table.
* `implied_ceiling()`'s `lower` and `upper` are renamed
  `ceiling_at_alpha_lower` and `ceiling_at_alpha_upper`: they hold the rate at
  its estimate and are not a confidence interval for the ceiling.
* `confint(fit, "alpha")` gives the profile interval of a `cpb()` fit without a
  bootstrap; coefficient intervals still need one.

## One bootstrap

* Every class bootstraps through one engine: the same resampling (rows,
  clusters, or fixed-effects units, a unit drawn twice relabeled as two), one
  layout of the stored draws (named columns, the dispersion parameter on its
  natural scale, and a fixed-effects fit's reference unit effect), replicates
  dropped when their refit fails or does not report convergence, and the
  attributes `nboot_ok`, `B` and `cluster`. The draws of `cpb_fe()`, `gec_fe()`
  and `zi_cpb()` now carry the dispersion parameter, and `zi_gec()` stores
  `delta` rather than its logarithm. A `zi_cpb()` or `zi_gec()` replicate runs
  the separation check's restart as the fit does, and relabels the zero
  equation's unit dummies with the count equation's units. Objects saved by
  earlier versions are read in their own layout. The intervals of `confint()`
  follow the same rules as before.
* A serial bootstrap or parametric bootstrap (`dispersion_test()`,
  `calibrate_alpha()`, `ud_screen()`) that will take over a minute says so, with
  the time its first replicate took. Fewer than 100 bootstrap replicates bring a
  note on the Monte Carlo error of the standard errors. `print()` of
  `dispersion_test()` says when its p-value sits at the floor 1/(B + 1).
* With `bias_correct = "jackknife"` and `se = "bootstrap"`, `cpb_fe()` and
  `gec_fe()` (and a `hurdle_cpb()` intensity) correct every bootstrap replicate
  the same way, by its two half-panel refits, with the validity gate checked once
  on the data, so the standard errors include the correction's own variance.
  The replicates used to be uncorrected fits whose interval was recentred on
  the corrected estimate. On the floor support they still are. In the
  package's simulations (40 units) the slope interval's coverage rose from 0.81
  to 0.88 at T = 6 and from 0.89 to 0.94 at T = 12; `?cpb_fe` also reports
  what the correction costs in the slopes' variance.
* `cpb_fe()` and `gec_fe()` sort a `cluster` given as a vector together with
  the data. It was not reordered when the rows were sorted by unit, which put
  the wrong clusters on any panel not supplied unit by unit.
* A `cluster` given as a factor with levels the data do not hold no longer
  draws those levels as empty clusters, so a replicate has as many clusters
  with rows as the data. Such a factor (a subset of a larger data set, or rows
  dropped for missing values) used to shrink and vary the replicates' size.

## Separation, the hurdles' participation equation, and `zi_test()`

* `zi_cpb(fe =)` leaves a unit with no positive count out of the estimation,
  as fixest leaves such units out of a Poisson fixed-effects fit, and says how
  many it left out. At the maximum such a unit's count effect is minus infinity,
  so its zeros are certain whatever its structural-zero probability and add
  nothing to the likelihood. Its rows stay in the fit at that limit, with
  fitted means of zero, and a structural-zero probability of `NA` where
  `zero_fe` would give the unit its own dummy. Before, the unit's count effect
  was filled with a value of the fit's choosing (the smallest estimated effect
  minus 2 under `method = "ml"`, the mean effect under `method = "em"`), whose
  zeros then shaped the inflation estimates, and the two methods filled it
  differently. On the floor support the fill is kept, as before 0.2.0.
* The participation equation of `hurdle_cpb()`, `hurdle_gec()` and
  `hurdle_count()` is checked for separation by the zero-inflated fits' rule.
  With `cluster` its standard errors are cluster-robust (`hurdle_count()`
  follows its `se`, robust or clustered), in `summary()`, `vcov()`,
  `confint()`, `tidy()`, `irr()`, the tables and the first differences. They
  were the glm's model-based ones in every case.
* A separated term has no standard error in `summary()`, `tidy()` and the
  texreg tables. `compare_models()` gains the columns `truncated` and
  `separated`, and warns when a separated fit, or a zero-truncated likelihood
  next to untruncated ones, enters the ranking.
* `zi_test()` refuses covariate-dependent inflation, a fixed-effects
  `zi_cpb()`, and a separated zero equation, where its boundary p-value has no
  reference distribution. It used to warn about the first and accept the others.

## Defaults, tables and messages

* The help pages say what the coefficients are on each family's scale.
  `?irr.count` lists the families whose `exp(beta)` is a ratio of means (the
  Poisson, the negative binomial, the mean-parameterized COM-Poisson, the
  hyper-Poisson, the generalized Poisson at `lambda >= 0`) and those whose is
  not (the CPB, the underdispersed GEC, the rate-parameterized COM-Poisson,
  whose `exp(beta)` is roughly the mean ratio to the power `nu`, the EWP, the
  gamma count, the double Poisson); `print()` of a `count_reg()` fit names the
  scale for the latter, and `irr()` names their ratio for what it is
  (`"lambda ratio"` for the rate-parameterized COM-Poisson, `"centrality ratio"`
  for the EWP, `"long-run rate ratio"` for the gamma count, `"location ratio"`
  for the double Poisson, `"mu ratio"` for the generalized Poisson below zero),
  with a message pointing to `ame()`; the CPB and the GEC keep `"IRR"`, a ratio
  of their rate. `?count_reg` says that the robust Poisson is consistent
  for the conditional mean whatever the distribution (Gourieroux, Monfort and
  Trognon 1984; Wooldridge 1999) and is the benchmark for mean effects, and its
  summary of a Poisson fit reads the gap between robust and model-based errors
  as a misspecified variance. The generalized Poisson's mean departs from
  `exp(x'b)` when underdispersed at small means; the hyper-Poisson is the
  constant-dispersion model of \enc{Sáez-Castillo}{Saez-Castillo} and
  \enc{Conde-Sánchez}{Conde-Sanchez} (2013) (their model with dispersion
  covariates is `DGLMExtPois::glm.hP()`), and its underdispersion is bounded
  (Var/Mean at least 0.5 at a mean of 2, 0.9 at 10), with its `gamma -> 0`
  limit stated for means below one. The families' sources are cited, and
  `dispersion_profile()` is said to show which variance function the data
  follow, not the process behind them. The fixed-effects pages no longer say
  the covariate coefficients are unaffected by the incidental-parameters bias,
  and `?hurdle_cpb` notes that bias in a participation logit with `part_fe` and
  points to a conditional logit (`survival::clogit()`) or the bife package for
  consistent participation slopes.
* `?cpb` and the vignette's new section "Mean effects and the robust Poisson"
  say that the CPB's mean effects are consistent only when the CPB is the data's
  process. On underdispersed counts from other processes (binomial counts with
  varying trials, a COM-Poisson) the package's checks put its average marginal
  effects 5 to 7 percent off at n = 5,000, where the robust Poisson's were
  within 1 percent.
* `?implied_ceiling` says that the ceiling is the fitted distribution's
  support end, set by its variance-mean link, and a bound on the counts only
  when the CPB is the data's process (on binomial counts out of six trials the
  fitted ceiling ran from 3 to 12); the README describes it the same way, and
  `?ceiling_exceedance` and `?compare_dispersion` say that the held-out
  exceedance check is one-sided: it detects a ceiling that is too low, not one
  that is too high.
* `tidy()` gains `scale`, what each coefficient multiplies (`"log mean"`, the
  family's own parameter such as `"log centrality"`, or a binary equation's
  link), and `glance()` gains `coef_scale`; a texreg column is named by its
  coefficient scale unless `custom.model.names` replaces it. A table row can
  hold coefficients on different scales; `ame()` puts every model's effect on
  the mean.

* `bias_correct = "jackknife"` needs `time =` in `cpb_fe()` and `gec_fe()` (and
  `hurdle_cpb()`, which passes it on), the column that orders each unit's rows;
  the half-panels and the residual trend check follow it. Rows not supplied in
  time order used to be split as they came, which defeated the validity gate.
* `max.support = NULL`, meaning `max(500, 10 * max(y))`, is now the default of
  `gec()`, `hurdle_gec()`, `zi_gec()`, `zi_cpb()` and `compare_dispersion()`, as
  it was of `cpb()`, `cpb_fe()` and `gec_fe()`; they used 500, and `gec()`
  failed on data with a count far above it. On the floor support `zi_cpb()` and
  `compare_dispersion()` keep 500, as before.
* `tidy()` takes `conf.int` for every class. `glance()` and the texreg tables
  report the dispersion parameter under its own name and with the side of its
  Poisson value that is underdispersion (the CPB's `alpha` below 1 and the
  gamma-count's `alpha` above 1 are both underdispersion), and `glance()` adds
  BIC. `augment()` takes `data =`.
* `residuals(type = "pearson")` for `cpb_fe()`, `gec()`, `gec_fe()` and
  `count_reg()`, from the exact variance of the fitted distribution.
* `count_reg()` with `se = "robust"` or `"cluster"` also keeps the model-based
  standard errors, and `summary()` prints the ratio of the two: King and Roberts
  (2015) read a large gap as a misspecification to fix. `?cpb` and `?count_reg`
  say what robust standard errors protect and what they do not.
* Failures are reported, not masked: a covariance the numerical derivatives
  cannot form, or a singular information, warns; `compare_dispersion()` has a
  `converged` column, and `ud_screen()` says which comparators it did not fit
  and why; a non-finite gradient of the concentrated fixed-effects objective
  falls back to central differences (it was replaced by zeros, which stopped
  BFGS where it stood); the matched families' fits report convergence when the
  BFGS stage or its polish converged (the polish used to set it regardless);
  and a texreg registration failure is reported at load.
* `summary()` counts the rows dropped for a missing value; `print()` of a CPB
  fit says that `alpha` is the variance-to-mean ratio at an integer ceiling;
  `count_reg(family = poisson)` takes the family object; a `|` inside `I()` is
  no longer read as a two-part formula; `cpb_fe()` and `gec_fe()` say where
  two-way fixed effects go.
* `?cpb` gives the notation of King (1989) and of Winkelmann, Signorino and
  King (1995), and why standard errors from the Hessian fail: an observation
  whose ceiling lies a distance `eps` above an integer carries information of
  order `1/eps` about its rate.

## Fixes

* On King's support, the profile interval for `alpha` follows the profile
  below 0.005, in geometric steps down to 1e-6, and its limit is found to a
  thousandth of its value. A lower side that never reaches its cut reaches
  alpha's limit, 0 (the two-point distribution), which `summary()` says; an
  upper side that never does ends at 0.999 or at an estimate beyond it. The trace
  used to stop at 0.005, so an estimate between 0.005 and about 0.02 got 0.005
  as its lower limit with the note "lower limit at the parameter bound", and an
  estimate below 0.005 an interval that did not contain it. The floor support
  keeps its trace and its interval.

* The standard errors of the exponentially weighted Poisson (`"ewp3"`,
  `"ewp2"`; `count_reg()`, the hurdle intensity, `zi_count()`) come from the
  expected information, summed over each observation's support. Its
  log-likelihood is continuous but kinked (its weight switches at
  \eqn{y = \lambda} and its normalizer at every integer \eqn{\lambda}), and the
  finite-difference Hessian used before measured only the local curvature of
  the pieces: with a covariate, the standard errors were erratic and mostly far
  too small (a slope's by up to an order of magnitude). The robust and
  clustered sandwiches use the same information and scores differenced at a
  small step. The message for a covariance that is not positive definite no
  longer blames fixed effects alone.

* Delta-method intervals (`first_difference()`, `predict(se.fit =, interval =)`,
  `ame()` and the hurdles' participation equation) differentiate along each
  parameter in steps of its standard error. The steps were relative to the
  parameter's value, so a coefficient on a covariate in small units (a
  population in persons) gave intervals of 1e80 or `NA`, and one in moderate
  units (an income in dollars) intervals off by a few percent; they now do not
  depend on the covariates' units.

* `zi_cpb(fe =, method = "ml")`, the default, is a joint maximum of the mixture
  likelihood on King's support: its plug-in fit (the unit effects of a
  first-stage fit to the positive counts, with a common shift now free) starts
  the observed-data polish of the EM path, which frees every unit's count
  effect. The plug-in alone was not a maximum. A unit whose positive counts are
  all 1 got an arbitrary effect from the first stage's flat set, the reported
  log-likelihood (which `compare_models()` ranks) could lie tens of points below
  the mixture's maximum, and the fit depended on where the covariates were
  centred. On a 40-unit panel the polish raised the log-likelihood from -394.2
  to -387.5, the EM path's value, and three centrings of the covariates agreed
  to 0.03 instead of 4. The floor support keeps the plug-in fit of earlier
  versions.

* A fixed-effects unit with no positive count sits 30 log units below its
  largest linear predictor (`cpb_fe()` on King's support, `gec_fe()`, and the
  left-out units of `zi_cpb(fe =)`), so its rates are at most `exp(-30)` whatever
  the covariates' scale. The sentinel was an absolute -30: with an uncentred
  covariate (a year trend) such a unit's rates were as large as `exp(53)`, its
  fitted values `NA`, and the analytic gradient unusable, which moved the
  estimates in the third digit. `gec_fe()`'s search for a unit's effect at
  `delta >= 1` also starts relative to the unit's Poisson intercept, not at an
  absolute -30 that could lie above the optimum.

* `gec_fe()` refuses a likelihood term whose sum the support guard cuts while
  its terms are still large, as `gec()` does, instead of renormalizing over the
  cut support (which overstated the likelihood toward large `delta`), and finds
  no ceiling crossings past the guard.

* On King's support, `zi_cpb(method = "em")` reports convergence as it happened
  (the EM's tolerance met, the observed-data polish within its round cap), so
  the bootstrap drops a replicate that did not converge; it used to report
  `TRUE` always, as it still does on the floor support.

* The gamma count's log-probabilities are computed from log tails, so they no
  longer underflow to `-Inf` at large means (a zero at a mean of 400 has log
  probability about -793); its moments refuse at the support cap, as the other
  families do, instead of reporting a mean summed over half the mass, and its
  variance no longer overflows past a mean of about 44,000.

* The variances of the COM-Poisson, the double Poisson, the generalized Poisson
  and the EWP are computed about the mean in a second pass. The one-pass sums
  cancelled for concentrated distributions (100% errors at a COM-Poisson
  `nu` of 1e4 and a double-Poisson `theta` of 100). An EWP sum that the support
  cap stops while its terms still rise is refused (`NA`) instead of returned as
  complete.

* The mean-parameterized COM-Poisson finds its rate at a mean below 1 and a
  large `nu` (a Bernoulli's root, from which it now starts); it returned `NA`
  there from `nu` of about 2000.

* `compare_models()` warns when the models were fit to different responses,
  not only to different numbers of observations.

* Parallel bootstraps (`cores > 1`) load the package on the workers from the
  library the session loaded it from.

* `first_difference(B > 0)` on `hurdle_cpb()` and `zi_cpb()` fits refit each
  replicate with the fit's own settings (fixed effects, link, offset, weights,
  support, estimation method), and resamples whole fixed-effects units, or the
  fit's clusters, when it has them. It used to refit the formulas alone.
* The COM-Poisson's `nu` and the hyper-Poisson's `gamma` no longer run off
  toward a limit of the family (the geometric, the Bernoulli, or 1 + Poisson)
  on data more extreme than the family. The search holds them in `[1e-4, 1e4]`,
  a fit on that edge is compared with one at the limit itself (`1e-8` or
  `1e8`), and `print()` notes a value at or past 1e-4 or 1e4.
* The mean-parameterized COM-Poisson (`family = "mpcmp"`) no longer evaluates
  a truncated normalizing sum. Its root search refuses a mean that the
  1e5-term support cap cannot represent, and the likelihood treats such a point
  as infeasible. An optimizer step to such a mean used to run 200 Newton steps
  of 1e5 terms per observation, at an overstated likelihood: a fit to
  zero-heavy overdispersed data that did not finish in half an hour now takes
  about a minute, and fits that never reach the cap are unchanged. For
  `nu < 1` the root search starts between the Poisson and the geometric rates
  rather than from Shmueli et al.'s approximation, whose correction diverges as
  `nu` approaches 0 (an overdispersed fit: 136 s to 4 s, the same maximum).
* The other families whose normalizer is a sum over the support (the
  rate-parameterized COM-Poisson, the double Poisson, the generalized Poisson
  below zero, the EWP) refuse it in the same way: where the 1e5-term cap would
  cut off part of the sum, a mean of order 1e5 or more, the normalizer is `NA`
  and the likelihood treats the point as infeasible. They used to return the
  truncated sum, and so an overstated likelihood. `dcompois()` and its siblings
  return `NA` there. Everywhere else the values are unchanged.
* The generalized Poisson with `lambda >= 0` uses its closed-form normalizer,
  1, in place of a sum over the capped support, which fell short where the
  tail is long: at `lambda = 0.99` the log-likelihood per observation was
  overstated by 1e-6 at a mean of 100 and by 0.03 at a mean of 50,000.
* The generalized Poisson's predicted probabilities no longer fail for a
  slightly negative `lambda` at large means, where the end of the support
  passed the integer range.
* A fit made with an offset carries it to new rows. `predict()`, `score()`,
  `cv_score()` and `ceiling_exceedance()` read the fit's offset column from
  `newdata`, or take `offset =`, and refuse rows without one. They used to
  evaluate new rows at exposure 1. Every `predict()` method now takes
  `offset =`, and the first differences of `hurdle_gec()`, `zi_gec()`,
  `hurdle_count()`, `zi_count()` and `zi_cpb()` hold the offset at its mean, as
  the other classes did.
* In first differences and in predictions for new rows, the reference unit of a
  fixed-effects fit averages only the units whose effect the data identify.
  Excluded are units with no positive count, whose effect is minus infinity,
  and, under King's zero-truncated support, units whose counts are all 1. Their
  sentinel or arbitrary effects used to move the reference, by an order of
  magnitude on panels with a few all-zero units. The same rule holds the unit
  dummies of `count_reg()`, `hurdle_count()`, `zi_count()` and the hurdle
  participation equations.
* `cpb_fe()`, `gec_fe()` and `zi_cpb(fe =)` return `fitted()`, `residuals()`,
  `predict()` and `simulate()` in the data's row order. They used to return
  them sorted by unit.
* `first_difference()` refuses a variable that also enters another column of the
  model (`I(x^2)`, `log(x)`, `poly(x, 2)`, a spline basis, an interaction). Only
  interactions were refused before, and the others gave a silently wrong
  contrast. `hurdle_gec()` and `zi_gec()` now refuse interactions too.
* `first_difference()` on a `hurdle_cpb()` fit uses the participation link it
  was fitted with. It used the logit under every link.
* `irr()` reports odds ratios only for a logit participation or inflation
  equation. Under probit or cloglog, exp(coefficient) is not an odds ratio and
  is no longer reported as one.
* `count_reg(fe =)`, `hurdle_count(fe =)` and `zi_count()` keep a transformed
  response (`pmin(y, 3) ~ x`) and accept `y ~ .`. They used to replace the
  response by its first variable, or fail.
* `alpha_confint()` refuses a `cpb_fe()` fit, whose profile would need its unit
  effects re-concentrated. It used to return an interval far too narrow.
* `zi_cpb(fe =)` stops when its first-stage fixed-effects fit fails. It used to
  continue with every unit effect at zero. Its two-step path counts in `df`
  only the unit effects the first stage estimated, not the fill values of units
  with no positive count.
* `zi_cpb(method = "em")` reports the observed-data log-likelihood at its
  estimate. Its fixed-effects polish maximizes that likelihood, with the
  probabilities floored only against underflow. A floor at 1e-12 on each
  probability used to overstate the log-likelihood whenever an observation's
  probability fell below it.
* `zi_cpb()`'s warning that the count component cannot cover the data no longer
  fires on legitimate fits whose zeros are all structural at large rates.
* `zi_test()` refuses a nest fitted on different rows, a zero-truncated nest, or
  a nest on the other CPB support (an object saved before 0.2.0 reads as the
  floor support). It used to warn, or to say nothing, and return a p-value.
* `rgec()` and `qgec()` return `NA` with a warning where the support guard
  binds, as `dgec()` does. They used to return the cut point. A non-finite rate
  gives `NA` in the GEC's moments, not a point mass at zero.
* On King's support, the R-side CPB distribution functions, predictions of
  `type = "prob"`, `simulate()` and the draws of `calibrate_alpha()` no longer
  cut the support at `max.support`, which there bounds only the width of the
  likelihood's pmf. Past the 1e5-term cap they return `NA` with a warning,
  rather than renormalizing over a cut support.
* Character covariates predict on new data that hold only some of their values.
* `glance()` and the BIC of `compare_models()` count the weight total for a
  weighted fit, as `nobs()` does. `max.support` is limited to one below the
  integer maximum, where the compiled loops would overflow.
* The zero equation's reference in the first differences of `zi_cpb()` and
  `zi_gec()`, and the participation equation's in `hurdle_gec()`, hold the
  `zero_fe` or `part_fe` dummies at their means over the units the data
  identify, as the other two-part classes do.
* `print()` of a CPB fit shows the median largest feasible count, in place of
  the continuous ceiling. The help pages lose no text to escaped percent signs.

# underdisp 0.1.2

## Bug fixes

* Converting a fitted ceiling to an integer is now guarded against values outside
  the range of `int`. The CPB's ceiling is `lambda / (1 - alpha)` and the GEC's is
  `mu / (1 - delta)`; both diverge as the dispersion parameter approaches its
  equidispersed limit, which an optimizer reaches routinely, and the conversion
  was undefined behaviour there. CRAN's macOS sanitizer build reported it in
  `cpb_pmf.h` and `cpb_fe.cpp`. All six conversions of this kind, including two the
  sanitizer did not reach, now cap the ceiling at `max.support + 1` first, through
  a single rule in `src/support_cap.h`.

  A ceiling past `max.support` was already treated as infeasible, on both of the
  platforms whose out-of-range conversions differ, so fitted values are unchanged:
  every cross-sectional and fixed-effects CPB fit, every GEC fit without fixed
  effects, and the calibrated interval reproduce bit for bit. The exception is
  `gec_fe()` on data whose dispersion is driven to the Katz limit, where the
  breakpoint scan now runs from a defined bound and the optimizer can settle a few
  parts in 10^8 away.

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
  two-sided or directional alternatives by the signed root and the boundary
  mixture where the Poisson value is on the boundary (CPB, negative binomial);
  and the Cameron-Trivedi auxiliary-regression test for a Poisson fit. The
  p-value is a parametric bootstrap under the fitted Poisson (`B`, `cores`)
  whenever the Poisson value is on the boundary of the family's parameter space
  (CPB, negative binomial) or the fit has unit fixed effects; where it is
  interior, the asymptotic p-value is used unless the first-order bias from the
  null's estimated mean parameters, p / sqrt(2n) standard deviations, pushes a
  5% test's size above 6%. `summary()` of a CPB fit flags its asymptotic p-value.
* `ud_screen()`: in `ztp_threshold = "bootstrap"` mode the NB-vs-Poisson p-value
  is a parametric bootstrap under the fitted Poisson (`nb_boot`) when the mean
  model is within the comparator gates; the fast mode's asymptotic value is
  labeled conservative. The bootstrap restores the random-number state, so the
  at-risk threshold and every other result are unchanged.
* `dispersion_profile()`: the empirical conditional variance-to-mean ratio by
  bins of the fitted mean against the ratio each fitted family implies, with a
  `plot()` method, so the mechanism behind an underdispersed outcome (a hard
  ceiling, regular event timing, a soft tail) can be read off the data.
* `calibrate_alpha()`: an opt-in calibration of the interval for the CPB
  dispersion parameter. The support of the CPB moves with the parameters
  through the log-ceiling coefficients, which behave like endpoint parameters,
  and the estimate of `alpha` inherits a bias toward zero from them; the
  first-order (chi-square) profile interval therefore covers about 0.90 to
  0.94 in simulations from the CPB, with its misses on the upper side.
  `calibrate_alpha(fit, B, cores, seed)` stores the parametric-bootstrap
  distribution of the signed root of the profile likelihood ratio in the fit,
  and `alpha_confint()`, `confint()`, `implied_ceiling()` and `summary()` then
  report the interval cut at its quantiles (0.92 to 0.96 in the same
  simulations, misses balanced). The responses are drawn once on the calling
  process and the refits use no random numbers, so the result does not depend
  on `cores`. It is refused, with the reason, for a fit at the `max.support`
  guard, weights that are not whole numbers, and a fit saved by 0.1.0. With few
  distinct covariate patterns the first-order interval already covers at its
  level and the calibration is not needed. Without it the accessors report the
  first-order interval, and `attr(, "method")` and `summary()` say which one is
  shown.

## Weights, parallel bootstraps, and column scaling

* `weights =` (frequency weights) on every estimator; a weight of w is
  equivalent to w copies of the row in the likelihood, the information, the
  score products, and the bootstrap, and `nobs()` returns the weight total.
  In the parametric bootstraps (`dispersion_test()`, `calibrate_alpha()`) a row
  of weight w contributes w independent simulated responses; weights that are
  not whole numbers have no such replicate data, so `dispersion_test()` then
  reports the asymptotic p-value with a note.
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
* `predict()`, `score()`, and the QoI methods rebuild a `newdata` model matrix
  in the fit's own factor coding (stored levels and contrasts) and align it to
  the coefficients by name; a factor whose reference level was set through its
  contrasts attribute (as in AER's `NMES1988`) is no longer re-coded on
  `newdata`, and a coding mismatch is refused instead of multiplied by position.
* A formula `offset()` term and a pscl-style `y ~ x | z` formula are refused
  with a message naming `offset =`, `participation =`, and `zero =`; rows whose
  unit identifier is missing are dropped like other incomplete rows; the
  screen's generalized-Poisson comparator is the package's own family, so it
  can report underdispersion, and it is gated by `comp_max_par`/`comp_max_n`
  like the COM-Poisson (a full fit, not a diagnostic, on a dummy-heavy design); `print.gec()` labels the dispersion direction by
  the likelihood-ratio test; `count_reg(fe =)` hides the absorbed dummies and
  warns when `cluster` is the fixed-effects variable; standard errors from a
  Hessian that is not positive definite are `NA` with a warning; the
  generalized Poisson's support is capped for negative dispersion; `score()`
  honours frequency weights; the auxiliary test is weighted; the parallel
  driver draws its cluster seed eagerly, so consecutive `cores > 1` bootstraps
  differ.
* `converged` is stored by every class and reported by `print()`/`summary()`
  when it is `FALSE`; a matched-family fit whose optimizer never left the
  infeasible plateau (a fitted mean outside `[1e-12, 1e8]`) is reported as not
  converged with a warning, and its log-likelihood is `-Inf`, instead of
  returning the sentinel as a likelihood.
* The predictive distribution of a zero-truncated `gec()` or `count_reg()`
  fit, as `score()`, `rootogram()`, `pit_hist()`, and `cv_score()` use it, is
  the fit's own zero-truncated pmf normalized over the full support (in 0.1.0
  the Katz pmf was scored untruncated and the count families renormalized over
  the observed range), so the in-sample log score equals `-logLik/n` for every
  class; `fitted()` on a `hurdle_gec` is the exact conditional mean that
  `predict()` returns.
* `cv_score()`, `rootogram()`, `pit_hist()`, and `dispersion_profile()` apply
  frequency weights as `score()` does, so a row of weight w counts as w copies
  of itself: the rootogram's observed and expected frequencies, the PIT
  histogram, the dispersion profile's bins (equal weight, rows kept whole) and
  bin means, and, through a new `weights` argument, the mean held-out scores of
  `cv_score()` and `score(newdata =)`. `rootogram()` no longer fails when `kmax`
  exceeds the largest count.
* `pit_hist()` read the lower bound F(y - 1) out of row order wherever a zero
  count preceded a positive one, so its histogram was wrong for data with
  zeros; the bound is now looked up row by row.
* `ud_screen()` reports the NB-vs-Poisson likelihood ratio at its boundary
  value 0 when `glm.nb()` stops short of the Poisson limit, and prints "at the
  Poisson boundary" there rather than an overdispersion reading.
* The zero-truncated CPB at an implied ceiling below 1 is the point mass at 1
  in the pmf, the moments, the first differences, and the simulator alike (the
  zero-truncated distribution for every ceiling in [1, 2)), so the hurdle's
  predictions, scores, and rootograms are finite at rows whose covariates put
  the rate below a unit's feasibility floor.
* One rule names the coefficients of the two-part classes everywhere:
  `participation:`/`intensity:` for the hurdles and `count:`/`zero:` for the
  zero-inflated models, in `vcov()`, `confint()`, `tidy()`, and the texreg
  tables, so a term read off a table indexes the covariance and the interval.
* A vector `cluster` is checked for its length and for missing values by every
  estimator before it is used (the fixed-effects and two-part estimators
  recycled or padded a short vector in 0.1.0), and a row dropped for a missing
  unit identifier drops its entry of a vector `offset`, `weights`, or
  `cluster` with it; `hurdle_count()` selects the cluster-robust covariance
  when `cluster` is supplied, as its siblings do, and accepts a vector
  `cluster` aligned to the data; `zi_cpb()` and `zi_gec()` refuse a
  non-integer response instead of flooring it; `zi_cpb()`, `zi_count()`, and
  `zi_gec()` work on scaled covariates with a Nelder-Mead polish like every
  other estimator (the unscaled joint fit of `zi_count()` could stop at its
  starting values: on lme4's `grouseticks` it stopped 4.4 log-likelihood units
  below pscl's optimum); `summary.zi_count()` prints the
  stored inflation standard errors; every bootstrapping summary states how
  many replicates converged; `count_reg()` and `zi_count()` store
  `$fitted.values` on the mean scale and the natural parameter as `$mu`.
* The parallel workers of `cores > 1` attach the package, so a `cv_score()`
  fit function that names `cpb()` or `count_reg()` unqualified runs on them.
* The analytic, robust, and cluster covariances of the matched families are
  differenced at an absolute step of 1e-3 in the scaled parameterization
  (numDeriv's default relative step left the family's support at large fitted
  means and returned `NA` standard errors).
* `predict()` with `newdata` returns `NA` for a row whose covariate or offset is
  missing, as `predict.glm()` does, instead of failing inside the distribution
  routines. A column holding only missing values (which R reads as logical) is
  treated as a missing value of the fitted variable's type, and a variable
  supplied with a different type than in the fitting data is refused with a
  message that says so.

## Estimation and numerics

* The CPB profile-likelihood interval for `alpha` is traced by continuation
  with a feasibility-repaired start at each `alpha`, and the profile at each
  `alpha` is maximized with the fit's own multistart (BFGS, Nelder-Mead
  polish, perturbed restarts), so the lower limit is the re-maximized
  profile's crossing point; the interval is one-sided only when the profile
  has not dropped to the cut by `alpha = 0.005`. Each limit is the outermost
  crossing of its cut (the profile is a saw-tooth, so the set above a cut need
  not be connected), and a cut deeper than the stored trace, as at a level
  above 0.995, re-runs the trace instead of returning the parameter bounds.
* `rcpb(truncated = TRUE)` warns when a rate's ceiling `lambda / (1 - alpha)` is
  below 1: the zero-truncated CPB does not exist there, and the draw of 1 it
  returns has probability zero under the stated parameters.
* `max.support` defaults to `max(500, 10 * max(y))`; a fit whose largest
  implied ceiling reaches the guard warns and records `$support_binding`; a
  fit whose every start is infeasible errors instead of returning the starting
  values; a CPB fit that the guard stops short of its Poisson limit reports the
  likelihood-ratio statistic at its boundary value 0.
* The concentrated fixed-effects fits solve each unit's intercept exactly.
  The unit objective is a saw-tooth in the intercept (the support floor moves
  as the rate crosses integer multiples of `1 - alpha`, and the objective drops
  at every such breakpoint), so its supremum sits at an interior critical point
  of a tooth or at the left limit of a breakpoint; the inner search enumerates
  the breakpoints within a window of the envelope's peak (as wide as the
  widest tooth in the unit) and refines inside whole teeth by golden section;
  every evaluation searches each unit afresh, so the concentrated likelihood is
  a function of the parameters alone. It is continuous in the dispersion
  parameter but can jump in the slopes where a unit sits on its feasibility
  floor, so the outer optimizer is a deterministic multistart: BFGS with the
  analytic envelope gradient (the breakpoint a unit sits at is tracked in the
  gradient), a Nelder-Mead polish, and two perturbed restarts; `?cpb_fe`
  states the size of the alternative optima. The Katz family below `delta = 1` has a finite
  support too, but its likelihood only kinks at an integer ceiling (the
  entering point's mass grows from zero), so `gec_fe()` maximizes each unit's
  unimodal objective by golden section on an expanding bracket, with the
  feasibility floor exact (support `0..ceiling(mu/(1-delta))`). The compiled
  code asserts the
  unit index against the design, and the fits use only the estimation rows, so
  rows dropped for missing values no longer misalign the unit blocks.
* `cpb()` documents the discontinuity of the CPB log-likelihood and the
  deterministic multistart that maximizes it (BFGS from the Poisson solution
  at each `alpha.start`, Nelder-Mead polish, feasibility-repaired restarts);
  the mean-parameterized COM-Poisson tabulates its log-factorials.
* The CPB pmf is normalized by a recursion outward from the mode of its terms
  (consecutive terms differ by the factor `(n - k)/(k + 1) (1 - alpha)/alpha`,
  so the sum runs over the body of the distribution and stops where the terms
  fall below relative precision on either side), shared by the pooled
  likelihoods and the concentrated fixed-effects likelihood and gradient. A
  likelihood evaluation costs the width of the pmf rather than its ceiling: on
  counts in the hundreds one evaluation of the concentrated objective takes
  0.04 s instead of 2.6 s, and the pooled objective is 45 times cheaper. The
  R implementation of the pmf keeps the direct lgamma form as the independent
  check (agreement to 5e-10 across alpha from 0.01 to 0.999 and rates from
  0.05 to 2,000).
* The reported maximum of a pooled `cpb()` fit is the maximum of the fit's own
  profile in `alpha`: the profile is traced by continuation on both sides of the
  estimate (the slopes re-maximized at each step from the neighbouring solution
  and from the Poisson slopes, with a finer pass next to the estimate) until it
  has dropped four log-likelihood units, and a trace point above the
  multistart's value restarts the fit from there; the multistart itself runs
  from nine values of `alpha` and ends with an exact breakpoint search in
  `alpha` at the fitted slopes, since the teeth of the surface in `alpha` are
  thinner than any grid. `confint()` and `summary()`
  read the profile interval off the stored trace. On small samples the
  multistart alone could stop on a lower tooth with `alpha` off by 0.1; the
  remaining differences against a dense grid of starts are of the size of the
  jumps (below 0.005 log-likelihood units on 12 samples of 72 zero-truncated
  counts). Unit intercepts belong in `cpb_fe()`, which solves them exactly; a
  pooled fit with many dummy columns can stop short of it (0.14
  log-likelihood units on pscl's `prussian` with 14 corps dummies).
* The inner search of `cpb_fe()` locates the peak of a unit's smooth part to
  1e-6 (in 0.1.1's first candidate the twelve golden-section steps left an
  error of 1.6e-3, wider than a tooth at ceilings above a few hundred, so the
  wrong tooth could be searched) and searches the neighbouring tooth when the
  interior optimum sits at a tooth's end.
* The GEC support guard refuses instead of renormalizing: a rate whose recursion
  has not reached the tail by `max.support` is infeasible in the likelihood, a
  fit whose fitted support reaches the guard warns and records
  `$support_binding`, and `predict()`, `dgec()`, and `pgec()` return `NA` with
  a warning for such a rate (in 0.1.0 the pmf was silently renormalized over
  `0..max.support`, so a prediction far outside the fitted range was too small).
* The COM-Poisson support rule no longer clips the rate at 1e10 (with `mu =` at
  strong underdispersion the rate exceeds it and the pmf was wrong), and the
  d/p/q/r functions return `NA` for a non-finite rate.
* `.ud_nm_polish()` chains its restarts until a restart gains less than the
  tolerance, as documented (the first candidate stopped after one).
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
* Every estimator drops unused factor levels from its model frame, as `lm()`
  and `glm()` do, so data subset from a larger set no longer trip the
  rank-deficiency refusal on an empty level.
* A response whose counts are all zero still fits, with a warning that the rate
  runs to its lower bound and the dispersion parameter is not identified.
* `max.support`, `inner_it`, `alpha.start`, and `B` (with `se = "bootstrap"`)
  are checked on entry and refused by name before any fitting: a support below
  1 or beyond the integer range, a fractional or negative count, fewer than two
  bootstrap replicates, a start outside (0, 1).
* The support guard's warning comes once, from the fit the user asked for: the
  refits of a bootstrap, a bootstrap first difference, or the split-panel
  jackknife, the start values of the zero-inflated fits, the nest that
  `zi_test()` fits from a formula, and the hurdle arm of `compare_dispersion()`
  no longer repeat it, and a `gec()` or `gec_fe()` fit no longer issues it
  twice.

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
* `broom` `tidy()`/`glance()` methods for all model classes (`augment()` for `cpb`), and
  `modelsummary` / `texreg` table support.
