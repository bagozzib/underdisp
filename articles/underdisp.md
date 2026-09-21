# Detecting and Modeling Underdispersed Counts

`underdisp` provides tools for detecting and modeling *underdispersion*
in count data: the case where the conditional variance is below the
conditional mean, so counts cluster more tightly around their
expectation than a Poisson allows. The Poisson and negative binomial
defaults cannot represent it; the negative binomial in particular
collapses onto the Poisson when the data are underdispersed.

``` r

library(underdisp)
```

## Simulating an underdispersed count

We generate a count with a conditional variance-to-mean ratio of about
one half.

``` r

n <- 400
x <- rnorm(n)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
y <- rbinom(n, N, 0.5)
d <- data.frame(y = y, x = x)
c(mean = mean(y), var = var(y), ratio = var(y) / mean(y))
#>      mean       var     ratio 
#>  5.527500 10.435332  1.887894
```

## Screening

[`ud_screen()`](https://bagozzib.github.io/underdisp/reference/ud_screen.md)
returns a marginal verdict, and, for zero-inflated outcomes, an at-risk
verdict benchmarked against a *zero-truncated* Poisson (which is what
separates genuine underdispersion from the artifact of conditioning on
positive counts).

``` r

ud_screen(y ~ x, data = d, run_cpb = FALSE)
#> 
#> === Underdispersion screen ===
#> Formula: y ~ x 
#> N=400  mean=5.527  max=22  %zero=1.8  (unconditional var/mean=1.89)
#> 
#> MARGINAL verdict: UNDERDISPERSED 
#>    Pearson=0.564  prop.slope=-0.391 (p=<2e-16)
#>    NB vs Poisson LR = 0, at the Poisson boundary (p= 0.5 asymptotic, conservative at this boundary )
#> 
#> AT-RISK (y>0) verdict:UNDERDISPERSED  [n_pos=393]
#>    ZTP-Pearson = 0.569  (underdispersed if < 0.885, the calibrated 5% threshold)
#> 
#> Model comparison (log-lik):
#>  Poisson       NB       GP     COMP      CPB 
#> -798.901 -798.904 -770.154 -772.010       NA 
#> 
#> COM-Poisson (full data): nu = 1.84  (> 1 = underdispersed, soft tail)
```

## Fitting the continuous parameter binomial

[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md) fits
the CPB, with `truncated = TRUE` for the common case in which
underdispersion lives among the positive counts of a zero-inflated
outcome.

``` r

fit <- cpb(y ~ x, data = d[d$y > 0, ], se = "none")
summary(fit)
#> 
#> Continuous Parameter Binomial regression (zero-truncated)
#> N = 393    inference: none 
#> 
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  1.58279         NA      NA       NA
#> x            0.49215         NA      NA       NA
#> 
#> alpha = 0.5256   (first-order profile 95% CI: 0.487 to 0.613)
#> Implied ceiling lambda/(1-alpha): median 10.22   range 2.48 to 37.8 
#> logLik = -745.05    AIC = 1496.1 
#> LR vs ZT-Poisson (H0: alpha = 1): 58.66, p 9.3805e-15
#> Note: the p-value is asymptotic, which over-rejects in finite samples at this boundary; dispersion_test() gives the parametric-bootstrap p-value.
```

The dispersion parameter `alpha` summarizes the compression, and each
observation carries an implied ceiling `lambda / (1 - alpha)`. The CPB
is a binomial with a non-integer number of trials, renormalized over its
finite support, so its mean is not exactly the rate `exp(x'b)`:
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
`predict(type = "response")` report the exact mean of the fitted
distribution (the conditional mean given `Y >= 1` for a zero-truncated
fit), and `predict(type = "rate")` the rate.

## Quantities of interest

Predicted means, the implied ceiling, and first differences are
available for user-specified covariate profiles.

``` r

predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "response")
#> [1] 3.027411 4.875001 7.964198
implied_ceiling(fit, newdata = data.frame(x = 0))
#>     lambda  ceiling    lower    upper
#> 1 4.868501 10.26202 9.491782 12.57789
first_difference(fit, "x", from = -1, to = 1)
#>  component  from    to  diff
#>       mean 3.027 7.964 4.937
```

## Testing equidispersion

The Poisson is the equidispersed member of every family in the package
(`alpha = 1` for the CPB, `delta = 1` for the GEC, `nu = 1` for the
COM-Poisson, and so on), so
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md)
tests a fitted model’s dispersion parameter against its Poisson value by
a likelihood-ratio test on the fitted design, with fixed effects carried
into the null. Where the Poisson value sits on the boundary of the
parameter space (the CPB, the negative binomial) the test is one-sided,
and its p-value comes by default from a parametric bootstrap under the
fitted Poisson, because the asymptotic Self–Liang mixture over-rejects
in finite samples; `B = 0` gives the quick asymptotic value shown here.
Elsewhere the alternative can be two-sided or directional. For a plain
Poisson fit the Cameron–Trivedi auxiliary regression is available as
`method = "auxiliary"`.

``` r

dispersion_test(fit, B = 0)   # asymptotic value; the default is the bootstrap p-value
#> 
#> Likelihood-ratio test of equidispersion (CPB vs Poisson; boundary null; asymptotic 1/2 chi2_0 + 1/2 chi2_1 mixture)
#> 
#> data: y ~ x
#> LR = 58.6579  (logLik: fitted family = -745.05, Poisson = -774.38),  p-value = 9.381e-15
#> alternative hypothesis: underdispersion
#> note: the p-value is asymptotic; at a boundary null the statistic is skewed toward rejection in finite samples, and the default (B = NULL) calibrates it by parametric bootstrap
#> estimate: alpha = 0.5256  (Poisson value 1)
dispersion_test(count_reg(y ~ x, data = d, family = "poisson"),
                method = "auxiliary", alternative = "under")
#> 
#> Cameron-Trivedi auxiliary regression test of equidispersion (Var = mu + a mu^2)
#> 
#> data: y ~ x
#> t = -9.3813,  p-value = < 2.2e-16
#> alternative hypothesis: underdispersion
#> estimate: a = -0.0619  (Poisson value 0)
```

## Bootstrap inference

Because the CPB’s support depends on its parameters, Hessian-based
standard errors are unreliable; coefficient inference uses a
cold-multistart pairs bootstrap (validated to nominal coverage in the
companion paper), and the dispersion parameter carries a
profile-likelihood interval. That interval is first-order by default:
the estimate of `alpha` is biased toward zero (the log-ceiling
coefficients behave like endpoint parameters), so it covers about 0.88
to 0.95 in simulations, with its misses on the upper side.
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md)
is the opt-in remedy: it stores the parametric-bootstrap distribution of
the signed root of the profile likelihood ratio in the fit
(`fit <- calibrate_alpha(fit, B = 199, cores = 4, seed = 1)`, about a
third of a fit per replicate), after which
[`alpha_confint()`](https://bagozzib.github.io/underdisp/reference/alpha_confint.md),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`implied_ceiling()`](https://bagozzib.github.io/underdisp/reference/implied_ceiling.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) report the
calibrated interval. Both intervals are model-based and assume
independent observations. `cores =` runs the replicates on a socket
cluster seeded from the session.

``` r

fit_b <- cpb(y ~ x, data = d[d$y > 0, ], se = "bootstrap", B = 49)   # a small B keeps the vignette quick
summary(fit_b)
#> 
#> Continuous Parameter Binomial regression (zero-truncated)
#> N = 393    inference: bootstrap 
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) 1.582786   0.017365  91.150 < 2.2e-16 ***
#> x           0.492152   0.017873  27.537 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> alpha = 0.5256   (first-order profile 95% CI: 0.487 to 0.613)
#> Implied ceiling lambda/(1-alpha): median 10.22   range 2.48 to 37.8 
#> logLik = -745.05    AIC = 1496.1 
#> LR vs ZT-Poisson (H0: alpha = 1): 58.66, p 9.3805e-15
#> Note: the p-value is asymptotic, which over-rejects in finite samples at this boundary; dispersion_test() gives the parametric-bootstrap p-value.
#> (49 bootstrap resamples converged)
irr(fit_b)             # rate ratios with percentile intervals
#>         term equation ratio estimate lower upper             method
#>  (Intercept)    count   IRR    4.869 4.732 5.056 bootstrap (stored)
#>            x    count   IRR    1.636 1.581 1.689 bootstrap (stored)
confint(fit_b)         # coefficients (percentile) and alpha (first-order profile likelihood)
#>                  2.5%     97.5%
#> (Intercept) 1.5543877 1.6206530
#> x           0.4579820 0.5243323
#> alpha       0.4870825 0.6129318
```

## The free-dispersion GEC

[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md) fits
King’s generalized event count (Katz) model, whose dispersion `delta` is
estimated freely, so the data choose the direction of dispersion rather
than the analyst presuming it. On the underdispersed count above it
recovers `delta` well below one; on a Poisson outcome it sits at one.

``` r

gec(y ~ x, data = d, se = "none")                                  # delta ~ 0.5
#> Generalized event count (Katz family) regression
#> Call:  gec(formula = y ~ x, data = d, se = "none")
#> (Intercept)           x 
#>      1.5768      0.4964 
#> 
#> dispersion delta (Katz; Var/Mean on an unbounded support) = 0.553  [underdispersed (LR p = 2.3e-14)]
#> logLik = -769.78,  n = 400
gec(y ~ x, data = data.frame(y = rpois(n, exp(1 + 0.4 * x)), x = x),
    se = "none")                                                   # delta ~ 1
#> Generalized event count (Katz family) regression
#> Call:  gec(formula = y ~ x, data = data.frame(y = rpois(n, exp(1 + 0.4 *     x)), x = x), se = "none")
#> (Intercept)           x 
#>      0.9971      0.4110 
#> 
#> dispersion delta (Katz; Var/Mean on an unbounded support) = 0.946  [equidispersion not rejected (LR p = 0.42)]
#> logLik = -739.98,  n = 400
```

The GEC carries the same zero-truncated, hurdle
([`hurdle_gec()`](https://bagozzib.github.io/underdisp/reference/hurdle_gec.md)),
zero-inflated
([`zi_gec()`](https://bagozzib.github.io/underdisp/reference/zi_gec.md)),
and fixed-effects
([`gec_fe()`](https://bagozzib.github.io/underdisp/reference/gec_fe.md))
variants as the CPB.

## High-dimensional fixed effects

Underdispersion is typically a *within-unit* phenomenon that pooled
analyses hide.
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)
absorbs a full set of unit fixed effects by concentrating them out of
the likelihood, so it scales to thousands of units.

``` r

panel <- do.call(rbind, lapply(1:30, function(i) {
  xx <- rnorm(12); NN <- pmax(round(exp(rnorm(1, 0, 0.4) + 0.4 * xx) / 0.5), 1)
  data.frame(unit = i, x = xx, y = rbinom(12, NN, 0.5))
}))
fe_fit <- cpb_fe(y ~ x, data = panel, fe = "unit")
fe_fit
#> CPB regression with 30 unit fixed effects (concentrated likelihood)
#> Coefficients:
#>     x 
#> 0.296 
#> 
#> alpha (shape parameter): 0.4596   median implied bound: 2.3 
#> Note: alpha is subject to incidental-parameters bias for short panels; see ?cpb_fe.
dispersion_test(fe_fit, B = 0)   # the statistic against a Poisson with the same unit effects
#> 
#> Likelihood-ratio test of equidispersion (CPB vs Poisson; boundary null; asymptotic 1/2 chi2_0 + 1/2 chi2_1 mixture)
#> 
#> data: y ~ x
#> LR = 151.9428  (logLik: fitted family = -356.55, Poisson = -432.52),  p-value = NA
#> alternative hypothesis: underdispersion
#> note: unit fixed effects bias the dispersion estimate toward underdispersion, so the asymptotic distribution does not apply; B > 0 gives the parametric-bootstrap p-value
#> estimate: alpha = 0.4596  (Poisson value 1)
```

With unit fixed effects the asymptotic distribution of that statistic
does not apply, so its p-value comes from a parametric bootstrap under
the fitted Poisson. Every replicate refits the model, which makes the
bootstrap about 199 times as slow as the fit; it is not run here:

``` r

dispersion_test(fe_fit, cores = 2)   # parametric-bootstrap p-value, 199 replicates
```

## Comparing the family

[`compare_dispersion()`](https://bagozzib.github.io/underdisp/reference/compare_dispersion.md)
fits the Poisson, negative binomial, the native soft-tail COM-Poisson,
the free-dispersion GEC, and the hard-ceiling CPB, and reports a fit
comparison plus the CPB’s ceiling-exceedance share.

``` r

compare_dispersion(y ~ x, data = d)$table
#>             df    logLik      AIC      BIC logscore       rps    zero_fit
#> Poisson      2 -798.9009 1601.802 1609.785 1.997252 0.9855061 0.024243957
#> NegBinomial  3 -798.9042 1603.808 1615.783 1.997261 0.9855095 0.024244780
#> COM-Poisson  3 -772.0096 1550.019 1561.994 1.930024 0.9631112 0.009312565
#> CPB          3 -769.2022 1544.404 1556.379 1.923006 0.9642686 0.010295571
#> GEC          3 -769.7769 1545.554 1557.528 1.924442 0.9644015 0.010538440
```

## The wider underdispersed family

Underdispersion has more than one mechanism, and
[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)
fits the distributions that express them through one interface: the
classical COM-Poisson (`"compois"`, rate-parameterized) and Huang’s
mean-parameterized COM-Poisson (`"mpcmp"`, a soft tail), the Consul–Jain
generalized Poisson (`"genpois"`, constant variance-to-mean ratio,
finite support when underdispersed), Winkelmann’s gamma-count
(`"gammacount"`, events with more regular timing than a Poisson
process), and Efron’s double Poisson (`"doublepois"`). Each inherits
fixed effects, zero-truncation, hurdle and zero-inflated forms, offsets,
frequency weights, robust and cluster-robust standard errors, and every
method in the package, so
[`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md)
places them next to the CPB on one footing.

``` r

dw <- d[1:150, ]                         # part of the sample keeps the seven fits quick
fits <- list(
  CPB          = cpb(y ~ x, data = dw, truncated = FALSE, se = "none"),
  Poisson      = count_reg(y ~ x, data = dw, family = "poisson"),
  NB           = count_reg(y ~ x, data = dw, family = "negbin"),
  `COM-Poisson`= count_reg(y ~ x, data = dw, family = "mpcmp"),
  GenPoisson   = count_reg(y ~ x, data = dw, family = "genpois"),
  GammaCount   = count_reg(y ~ x, data = dw, family = "gammacount"),
  DoublePois   = count_reg(y ~ x, data = dw, family = "doublepois")
)
do.call(compare_models, fits)
#>             df    logLik      AIC      BIC logscore      rps
#> CPB          3 -297.1212 600.2425 609.2744 1.980808 1.018955
#> GenPoisson   3 -297.4525 600.9050 609.9369 1.983017 1.018599
#> GammaCount   3 -299.0185 604.0369 613.0688 1.993456 1.018395
#> COM-Poisson  3 -299.2068 604.4135 613.4454 1.994712 1.019137
#> DoublePois   3 -299.9026 605.8053 614.8372 1.999351 1.019015
#> Poisson      2 -305.3268 614.6535 620.6748 2.035512 1.029690
#> NB           3 -305.3268 616.6535 625.6854 2.035512 1.029689
```

On underdispersed data the negative binomial collapses onto the Poisson,
while the underdispersed families capture the compression and win on AIC
and the proper scores. The families differ in how their dispersion moves
with the mean, and
[`dispersion_profile()`](https://bagozzib.github.io/underdisp/reference/dispersion_profile.md)
sets the empirical conditional variance-to-mean ratio, by bins of the
fitted mean, against the curve each family implies:

``` r

dispersion_profile(CPB = fits$CPB, NB = fits$NB, GammaCount = fits$GammaCount,
                   GenPoisson = fits$GenPoisson)
```

![](underdisp_files/figure-html/profile-1.png)

The correlated-random-effects device
([`mundlak()`](https://bagozzib.github.io/underdisp/reference/mundlak.md))
and matching `d`/`p`/`q`/`r` functions
([`dcpb()`](https://bagozzib.github.io/underdisp/reference/cpb-distribution.md),
[`dgec()`](https://bagozzib.github.io/underdisp/reference/gec-distribution.md),
[`dcompois()`](https://bagozzib.github.io/underdisp/reference/compois-distribution.md),
[`dgammacount()`](https://bagozzib.github.io/underdisp/reference/gammacount-distribution.md),
[`dgenpois()`](https://bagozzib.github.io/underdisp/reference/genpois-distribution.md),
[`ddoublepois()`](https://bagozzib.github.io/underdisp/reference/doublepois-distribution.md))
round out the family.

## Excess zeros: hurdle and zero-inflated models

Many count outcomes mix a participation process (most units at zero)
with a tight positive count. The bundled peacekeeping panel – the number
of UN operations each state contributes troops to per year – shows the
package’s central move: marginally the count looks overdispersed, but
conditioning on country fixed effects and benchmarking the positive
counts against a zero-truncated Poisson, the at-risk process is
underdispersed.

``` r

data(peacekeeping)
ud_screen(contributions ~ democracy + lgdppc + lpop + milper + factor(iso3),
          data = peacekeeping, run_cpb = FALSE, run_gp = FALSE)
#> 
#> === Underdispersion screen ===
#> Formula: contributions ~ democracy + lgdppc + lpop + milper + factor(iso3) 
#> N=4448  mean=3.191  max=18  %zero=41.5  (unconditional var/mean=4.56)
#> 
#> MARGINAL verdict: OVERDISPERSED 
#>    Pearson=1.086  prop.slope=0.191 (p=<2e-16)
#>    NB vs Poisson LR = 18.11 (p= 1e-05 asymptotic, conservative at this boundary ; sig => overdispersion)
#> 
#> AT-RISK (y>0) verdict:UNDERDISPERSED  [n_pos=2602]
#>    ZTP-Pearson = 0.841  (underdispersed if < 0.955, the calibrated 5% threshold)
#> 
#> Model comparison (log-lik):
#>   Poisson        NB        GP      COMP       CPB 
#> -6645.816 -6636.761        NA        NA        NA
```

That is the case for a two-part model with an underdispersed intensity.
[`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md)
joins a participation logit to a zero-truncated CPB, and
[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md)
fits the structural-zero mixture;
[`zi_test()`](https://bagozzib.github.io/underdisp/reference/zi_test.md)
and
[`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md)
adjudicate between them.

``` r

z <- rnorm(n)
yh <- rhurdle_cpb(n, lambda = exp(1.2 + 0.3 * x), alpha = 0.5,
                  p = plogis(0.3 + 0.8 * z))
dh <- data.frame(y = yh, x = x, z = z)
h <- hurdle_cpb(y ~ x, data = dh, participation = ~ z)
zi <- zi_cpb(y ~ x, data = dh, zero = ~ z)
compare_models(hurdle = h, mixture = zi)
#>         df    logLik      AIC      BIC logscore       rps
#> mixture  5 -601.2459 1212.492 1232.449 1.503115 0.9731981
#> hurdle   5 -602.5343 1215.069 1235.026 1.506336 0.9737676
```

The hurdle’s
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
separates the extensive and intensive margins exactly – which channel a
covariate moves, not just the blended marginal effect. One practical
note: because the hurdle factorizes, its participation stage is an
ordinary logistic regression; if a dummy-heavy participation equation
separates, fit that stage with a dedicated bias-reduction package
(`logistf`, `brglm2`) alongside this package’s zero-truncated intensity.

## Weights

Every estimator takes frequency `weights` (a weight of `w` is equivalent
to `w` copies of the row, in the likelihood, the information, and the
bootstrap), so tabulated data fit exactly as the expanded data would:

``` r

w <- sample(1:3, n, TRUE)
c(weighted = count_reg(y ~ x, data = d, family = "gammacount", weights = w)$loglik,
  expanded = count_reg(y ~ x, data = d[rep(seq_len(n), w), ], family = "gammacount")$loglik)
#>  weighted  expanded 
#> -1564.823 -1564.823
```

## Short panels: bias-corrected fixed effects

The concentrated fixed-effects dispersion estimate carries the
incidental- parameters bias of order 1/T: with few observations per
unit, `alpha` is biased *downward* (the panel looks more underdispersed
than it is). `bias_correct = "jackknife"` removes the leading bias term
by the split-panel jackknife of Dhaene and Jochmans (2015), refitting on
each unit’s temporal halves.

``` r

short <- do.call(rbind, lapply(1:16, function(i) {
  xx <- rnorm(8); NN <- pmax(round(exp(1.0 + rnorm(1, 0, 0.4) + 0.3 * xx) / 0.5), 1)
  data.frame(unit = i, x = xx, y = rbinom(8, NN, 0.5))
}))
ml <- cpb_fe(y ~ x, data = short, fe = "unit")
jk <- cpb_fe(y ~ x, data = short, fe = "unit", bias_correct = "jackknife")
c(ml = ml$alpha, jackknife = jk$alpha)   # truth is 0.5; ML is biased downward
#>        ml jackknife 
#> 0.4142559 0.5079890
```

The correction is only valid when the two half-panels estimate the same
parameter (the method’s time-homogeneity requirement), so it carries a
validity gate: the panel is also split cross-sectionally by units – a
placebo that is exchangeable under any time pattern – and if the
temporal halves disagree beyond that placebo noise, the correction is
*refused* with a warning naming the failed assumption and the
maximum-likelihood fit is returned. On a trending or regime-changing
panel, the refusal is the correct answer. The gate is deliberately
powered over sized: in calibration it refuses about 9% of genuinely
homogeneous panels (you keep the ordinary ML fit) while catching 98% of
dispersion regime changes and all smooth unmodeled trends.

## Simulated-residual diagnostics with DHARMa

Every fitted model in the package has a
[`simulate()`](https://rdrr.io/r/stats/simulate.html) method, so the
whole family plugs into `DHARMa`’s simulated-residual diagnostics.

``` r

sims <- simulate(h, nsim = 100, seed = 1)
res <- DHARMa::createDHARMa(simulatedResponse = as.matrix(sims),
                            observedResponse  = dh$y,
                            fittedPredictedResponse = fitted(h),
                            integerResponse = TRUE)
plot(res)
```

![](underdisp_files/figure-html/dharma-1.png)
