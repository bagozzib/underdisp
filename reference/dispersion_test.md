# Regression-adjusted test of equidispersion

Tests whether a fitted count model's dispersion departs from the
Poisson, conditional on the covariates. For a fit from a family with a
dispersion parameter (`method = "lr"`, the default), the test is the
likelihood-ratio test of that parameter against its Poisson value, with
the Poisson refit on the same design (fixed effects included); its
p-value comes from the asymptotic distribution or from a parametric
bootstrap under the fitted Poisson, by the rule in Details. For a
Poisson fit, `method = "auxiliary"` runs the Cameron–Trivedi (1990)
auxiliary regression, which needs no alternative family.

## Usage

``` r
dispersion_test(
  object,
  alternative = c("two.sided", "under", "over"),
  method = c("lr", "auxiliary"),
  B = NULL,
  cores = 1L
)
```

## Arguments

- object:

  A fitted `cpb`, `cpb_fe`, `gec`, `gec_fe`, or `count_reg` model (for
  `method = "auxiliary"`, a `count_reg` fit with `family = "poisson"`).

- alternative:

  `"two.sided"` (default where the null is interior), `"under"`
  (underdispersion), or `"over"` (overdispersion).

- method:

  `"lr"` (likelihood ratio against the Poisson refit; default) or
  `"auxiliary"` (the Cameron–Trivedi auxiliary regression for a Poisson
  fit).

- B:

  Parametric-bootstrap replicates for `method = "lr"`: `NULL` (the
  default) applies the calibration rule in Details (199 replicates when
  the bootstrap is needed), `0` requests the asymptotic distribution,
  and a positive whole number forces the bootstrap.

- cores:

  Worker processes for the bootstrap replicates.

## Value

An object of class `c("dispersion_test", "htest")` with the statistic,
the p-value, the estimated dispersion parameter and its Poisson value,
the two log-likelihoods, and the alternative; for `method = "lr"` also
`calibration` (`"asymptotic"` or `"parametric bootstrap"`), `B` and
`B_ok` (replicates requested and completed), `boot` (the simulated
statistics), `first_order_size` (the asymptotic 5\\ size),
`asymptotic_ok` (whether the calibration rule admits the asymptotic
distribution), and `note` where a statistic or p-value needs one.

## Details

The Poisson is nested in every family at one value of the dispersion
parameter: `alpha = 1` (CPB), `delta = 1` (GEC), `nu = 1` (both
COM-Poisson parameterizations), `lambda = 0` (generalized Poisson),
`alpha = 1` (gamma-count), `theta = 1` (double Poisson), and
`1/theta = 0` (negative binomial). When that value is interior to the
parameter space the likelihood-ratio statistic is asymptotically
\\\chi^2_1\\ under the null; a directional alternative (`"under"` or
`"over"`) uses the signed root \\r = \pm\sqrt{LR}\\, positive when the
estimate lies on the underdispersed side of the Poisson value, which is
asymptotically standard normal. When the Poisson value is on the
boundary (the CPB, whose `alpha` lives in (0, 1); the negative binomial,
whose overdispersion parameter is non-negative) the asymptotic null
distribution is the \\\tfrac12\chi^2_0 + \tfrac12\chi^2_1\\ mixture
(Self and Liang 1987; Andrews 2001) and the test is one-sided by
construction (`"under"` for the CPB, `"over"` for the negative
binomial). The CPB qualifies although its support moves with its
parameters: near `alpha = 1` its ceiling `lambda / (1 - alpha)`
diverges, the log-likelihood is regular on that side, and its score
there is the classical dispersion score \\-\[(y - \lambda)^2 - y\] /
(2\lambda)\\; the zero-truncated model reaches the same mixture through
its efficient score. Because the Poisson is the CPB's limit as `alpha`
approaches 1, the statistic cannot be negative in principle; when the
`max.support` guard stops `alpha` short of that limit and the fitted
log-likelihood lands just below the Poisson's, the statistic takes its
boundary value 0. Any other negative statistic means the optimizer fell
short, and no p-value is reported.

**Calibration.** Where the Poisson value is on the boundary of the
family's parameter space (the CPB, the negative binomial), or the fit
has unit fixed effects, the p-value is by default a parametric bootstrap
with 199 replicates. The one-sided boundary statistic is skewed toward
rejection in finite samples, and unit intercepts bias the dispersion
estimate toward underdispersion by an amount of order 1/T; in both cases
the asymptotic distribution over-rejects. Where the Poisson value is
interior, the asymptotic distribution is used unless the estimated mean
parameters shift it too far: to first order they move the signed root
toward underdispersion by \\p/\sqrt{2n}\\ standard deviations for \\p\\
mean parameters on \\n\\ observations (the weight total), the bias Dean
and Lawless (1989) remove from the score test, and the bootstrap takes
over when that shift pushes the first-order size of a 5\\ simulates
responses from the Poisson (zero-truncated Poisson) fitted to the same
design, offset, weights, and unit effects, refits both models to each,
and reports the share of simulated statistics at least as extreme as the
observed one, counting the observed one, among the replicates whose two
fits both succeeded. Simulating at the Poisson value itself keeps the
bootstrap valid at the boundary. `B = 0` requests the asymptotic
distribution (with unit fixed effects it returns the statistic without a
p-value), and a positive `B` sets the number of bootstrap replicates.
Every replicate refits both models, so the bootstrap takes about `B`
times as long as the original fit; `cores` spreads the replicates over
worker processes.

## References

Andrews, D. W. K. (2001). Testing when a parameter is on the boundary of
the maintained hypothesis. *Econometrica*, 69(3), 683-734. Cameron, A.
C. and Trivedi, P. K. (1990). Regression-based tests for overdispersion
in the Poisson model. *Journal of Econometrics*, 46(3), 347-364. Dean,
C. and Lawless, J. F. (1989). Tests for detecting overdispersion in
Poisson regression models. *Journal of the American Statistical
Association*, 84(406), 467-472. Self, S. G. and Liang, K.-Y. (1987).
Asymptotic properties of maximum likelihood estimators and likelihood
ratio tests under nonstandard conditions. *Journal of the American
Statistical Association*, 82(398), 605-610.

## See also

[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md),
[`dispersion_profile()`](https://bagozzib.github.io/underdisp/reference/dispersion_profile.md),
[`zi_test()`](https://bagozzib.github.io/underdisp/reference/zi_test.md)

## Examples

``` r
set.seed(2); n <- 400; x <- rnorm(n)
d <- data.frame(y = rgammacount(n, exp(1 + 0.4 * x), alpha = 2), x = x)
dispersion_test(count_reg(y ~ x, d, family = "gammacount"), alternative = "under")
#> 
#> Likelihood-ratio test of equidispersion (Gamma-count vs Poisson; asymptotic chi-square)
#> 
#> data: y ~ x
#> LR = 52.1022  (logLik: fitted family = -647.24, Poisson = -673.30),  p-value = 2.634e-13
#> alternative hypothesis: underdispersion
#> estimate: alpha = 1.9711  (Poisson value 1)
dispersion_test(count_reg(y ~ x, d, family = "poisson"), method = "auxiliary")
#> 
#> Cameron-Trivedi auxiliary regression test of equidispersion (Var = mu + a mu^2)
#> 
#> data: y ~ x
#> t = -9.4533,  p-value = < 2.2e-16
#> alternative hypothesis: two.sided
#> estimate: a = -0.1291  (Poisson value 0)
```
