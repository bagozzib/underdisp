# Regression-adjusted test of equidispersion

Tests whether a fitted count model's dispersion departs from the
Poisson, conditional on the covariates. For a fit from a family with a
dispersion parameter (`method = "lr"`, the default), the test is the
likelihood-ratio test of that parameter against its Poisson value, with
the Poisson refit on the same design (fixed effects included); the
alternative can be two-sided or directional (see Details). For a Poisson
fit, `method = "auxiliary"` runs the Cameron–Trivedi (1990) auxiliary
regression, which needs no alternative family.

## Usage

``` r
dispersion_test(
  object,
  alternative = c("two.sided", "under", "over"),
  method = c("lr", "auxiliary")
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

## Value

An object of class `c("dispersion_test", "htest")` with the statistic,
the p-value, the estimated dispersion parameter and its Poisson value,
the two log-likelihoods (LR method), and the alternative.

## Details

The Poisson is nested in every family at one value of the dispersion
parameter: `alpha = 1` (CPB), `delta = 1` (GEC), `nu = 1` (both
COM-Poisson parameterizations), `lambda = 0` (generalized Poisson),
`alpha = 1` (gamma-count), `theta = 1` (double Poisson), and
`1/theta = 0` (negative binomial). When that value is interior to the
parameter space the likelihood-ratio statistic is \\\chi^2_1\\ under the
null; a directional alternative (`"under"` or `"over"`) uses the signed
root \\r = \pm\sqrt{LR}\\, signed toward underdispersion when the
estimate lies on the underdispersed side of the Poisson value, which is
standard normal under the null, so the one-sided p-value is \\1 -
\Phi(r)\\ for `"under"` and \\\Phi(r)\\ for `"over"`. When the Poisson
value is on the boundary (the CPB, whose `alpha` lives in (0, 1); the
negative binomial, whose overdispersion parameter is non-negative) the
null distribution of the statistic is the \\\tfrac12\chi^2_0 +
\tfrac12\chi^2_1\\ mixture of Self and Liang (1987), the test is
one-sided by construction, and `alternative` is fixed accordingly
(`"under"` for the CPB, `"over"` for the negative binomial). For the
concentrated fixed-effects fits the Poisson null carries the same unit
effects, profiled out in closed form (Poisson) or by a per-unit Newton
step (zero-truncated Poisson).

## References

Cameron, A. C. and Trivedi, P. K. (1990). Regression-based tests for
overdispersion in the Poisson model. *Journal of Econometrics*, 46(3),
347-364. Self, S. G. and Liang, K.-Y. (1987). Asymptotic properties of
maximum likelihood estimators and likelihood ratio tests under
nonstandard conditions. *Journal of the American Statistical
Association*, 82(398), 605-610.

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
#> Likelihood-ratio test of equidispersion (Gamma-count vs Poisson)
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
