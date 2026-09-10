# Matched count regressions: Poisson, negative binomial, COM-Poisson, generalized Poisson, gamma-count, double Poisson

Fits a count regression with a log link from any of the package's count
families, optionally zero-truncated and/or with unit fixed effects,
returning an object that
[`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md),
[`score()`](https://bagozzib.github.io/underdisp/reference/score.md),
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md),
[`dispersion_profile()`](https://bagozzib.github.io/underdisp/reference/dispersion_profile.md),
and the broom methods treat on the same footing as a
[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md) or
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md) fit.
The families share one interface, so a gamma-count and a COM-Poisson can
be compared to a CPB with identical degrees-of-freedom, log-likelihood,
proper-score, weight, and robust-standard-error accounting rather than
reconciled across packages.

## Usage

``` r
count_reg(
  formula,
  data,
  family = c("poisson", "negbin", "compois", "mpcmp", "genpois", "gammacount",
    "doublepois"),
  truncated = FALSE,
  fe = NULL,
  offset = NULL,
  weights = NULL,
  se = c("analytic", "robust", "cluster", "none"),
  cluster = NULL
)
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame.

- family:

  One of `"poisson"`, `"negbin"`, `"compois"`, `"mpcmp"`, `"genpois"`,
  `"gammacount"`, `"doublepois"` (see Families).

- truncated:

  Logical; if `TRUE`, fit the zero-truncated form (requires all
  `Y >= 1`), the count analogue of the CPB's zero-truncated default.

- fe:

  Optional column name(s) for fixed effects, entered as factor dummies
  (so the degrees of freedom count each absorbed intercept, matching
  [`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)).
  Pass a vector of two columns for two-way (e.g. unit and time) fixed
  effects.

- offset:

  Optional offset entered on the linear-predictor (log) scale – the
  log-mean for the mean-parameterized families, the log-rate
  (\\\log\lambda\\) for `"compois"`. A numeric vector or the name of a
  column in `data`, e.g. \\\log(\text{exposure})\\ so the model becomes
  a rate model.

- weights:

  Optional frequency weights: a numeric vector or the name of a column
  in `data`. A weight of \\w_i\\ is equivalent to \\w_i\\ copies of row
  \\i\\ in the log-likelihood, the information, and the score products,
  and [`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the weight
  total. Survey (probability) weights call for `se = "robust"` or
  `se = "cluster"`.

- se:

  Standard errors: `"analytic"` (inverse information, default),
  `"robust"` (heteroskedasticity-consistent sandwich), `"cluster"`
  (cluster-robust; needs `cluster`), or `"none"`.

- cluster:

  Optional cluster identifier (a column name in `data` or a vector
  aligned to its rows) for `se = "cluster"`; supplying it selects
  `se = "cluster"` unless `se` is given explicitly.

## Value

An object of class `"count_reg"`. The component `$theta` holds the
dispersion parameter on its natural scale (the negative-binomial size,
the COM-Poisson \\\nu\\, the generalized-Poisson \\\lambda\\, the
gamma-count \\\alpha\\, the double-Poisson \\\theta\\; `NA` for the
Poisson); `$fitted.values` and `$residuals` are on the mean scale (the
conditional mean E(Y \| Y \>= 1) for a zero-truncated fit), as
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) return them, and
`$mu` is the family's natural parameter `exp(offset + x'b)`.

## Families

- `"poisson"`: the equidispersed baseline.

- `"negbin"`: negative binomial (overdispersion only); `$theta` is the
  size.

- `"compois"`: the classical Conway–Maxwell–Poisson in its rate
  parameterization, \\\log\lambda = x'\beta\\. The rate \\\lambda\\ is
  not the mean (though `predict(type = "response")` and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) return the
  mean); `$theta` is the dispersion \\\nu\\ (\\\nu \> 1\\
  underdispersed, \\\nu = 1\\ Poisson, \\\nu \< 1\\ overdispersed).

- `"mpcmp"`: the COM-Poisson in Huang's (2017) mean parameterization,
  \\\log\mathrm{E}(Y) = x'\beta\\, with the same dispersion \\\nu\\; the
  coefficients are effects on the log mean and rate ratios are exact.

- `"genpois"`: the Consul–Jain generalized Poisson with constant
  dispersion \\\lambda \in (-1, 1)\\ and mean \\\mu = \exp(x'\beta)\\,
  so \\\mathrm{Var}/\mathrm{Mean} = 1/(1-\lambda)^2\\; \\\lambda \< 0\\
  is underdispersion, on a finite support that the pmf is renormalized
  over (see
  [genpois-distribution](https://bagozzib.github.io/underdisp/reference/genpois-distribution.md)).

- `"gammacount"`: Winkelmann's (1995) gamma-count renewal-process model,
  underdispersed when the waiting times between events are more regular
  than exponential (\\\alpha \> 1\\); \\\exp(x'\beta)\\ is the long-run
  event rate and the exact mean is reported (see
  [gammacount-distribution](https://bagozzib.github.io/underdisp/reference/gammacount-distribution.md)).

- `"doublepois"`: Efron's (1986) double Poisson with the exact
  normalizing constant, \\\theta \> 1\\ underdispersed (see
  [doublepois-distribution](https://bagozzib.github.io/underdisp/reference/doublepois-distribution.md)).

The Poisson is the equidispersed member of every family, so
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md)
tests each family's dispersion parameter against it, and
[`dispersion_profile()`](https://bagozzib.github.io/underdisp/reference/dispersion_profile.md)
compares the families' implied variance-to-mean curves against the data.

Runtime: the Poisson, negative binomial, gamma-count, and double Poisson
and generalized Poisson families fit in under a second at a few hundred
rows; the COM-Poisson families evaluate a normalizing sum per
observation (and the mean parameterization solves a root per
observation), so they take seconds at a few hundred rows and minutes at
several thousand; `se = "robust"` and `"cluster"` add numerical scores
at the same cost per parameter.

## References

Huang, A. (2017). Mean-parametrized Conway–Maxwell–Poisson regression
models for dispersed counts. *Statistical Modelling*, 17(6), 359-380.
Winkelmann, R. (1995). Duration dependence and dispersion in count-data
models. *Journal of Business & Economic Statistics*, 13(4), 467-474.
Efron, B. (1986). Double exponential families and their use in
generalized linear regression. *Journal of the American Statistical
Association*, 81(395), 709-721. Consul, P. C. and Famoye, F. (1992).
Generalized Poisson regression model. *Communications in Statistics –
Theory and Methods*, 21(1), 89-109.

## See also

[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md),
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md),
[`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md),
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md),
[`dispersion_profile()`](https://bagozzib.github.io/underdisp/reference/dispersion_profile.md),
[`hurdle_count()`](https://bagozzib.github.io/underdisp/reference/hurdle_count.md),
[`zi_count()`](https://bagozzib.github.io/underdisp/reference/zi_count.md)

## Examples

``` r
set.seed(1); n <- 400; x <- rnorm(n)
y <- rgammacount(n, exp(1 + 0.5 * x), alpha = 2)     # regular event timing
d <- data.frame(y = y, x = x)
m <- count_reg(y ~ x, data = d, family = "gammacount")
m
#> Gamma-count count regression  
#> Call:  count_reg(formula = y ~ x, data = d, family = "gammacount")
#> 
#> Coefficients:
#> (Intercept)           x 
#>      0.9847      0.4955 
#> 
#> alpha (gamma-count dispersion; > 1 underdispersed): 1.9406  [Var/Mean ~ 0.515]
compare_models(poisson = count_reg(y ~ x, d, family = "poisson"), gammacount = m)
#>            df    logLik      AIC      BIC logscore       rps
#> gammacount  3 -637.0878 1280.176 1292.150 1.592719 0.6721447
#> poisson     2 -661.4718 1326.944 1334.927 1.653679 0.6905665
```
