# Compare fitted underdispersed-count models

Compares fitted models from this package — `cpb`, `cpb_fe`,
`hurdle_cpb`, `zi_cpb`, `gec`, `gec_fe`, `hurdle_gec`, `zi_gec`,
`count_reg`, `hurdle_count`, and `zi_count`, in any combination of fixed
effects and robust/clustered standard errors — on information criteria
and, where a predicted distribution is available, on proper scores.
Unlike
[`zi_test()`](https://bagozzib.github.io/underdisp/reference/zi_test.md)
this is not a hypothesis test: the models need not be nested, so it is
the appropriate tool for the non-nested hurdle-versus-mixture choice.
Passing an object that is not a fitted model from this package is an
error, and models fit on different data raise a warning.

## Usage

``` r
compare_models(...)
```

## Arguments

- ...:

  Two or more fitted models of the classes above, optionally named. With
  frequency weights, BIC uses the weight total
  ([`nobs()`](https://rdrr.io/r/stats/nobs.html)).

## Value

A data frame with `df`, `logLik`, `AIC`, `BIC`, and (where the predicted
distribution is available) `logscore` and `rps`, one row per model,
ordered by AIC, and the flags `truncated` (a zero-truncated likelihood,
conditional on `Y >= 1`) and `separated` (a zero-inflation or
participation equation the data do not bound, whose log-likelihood is a
supremum along a divergent ridge; see
[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md),
Separation). The function warns when the models mix truncated and
untruncated likelihoods, whose log-likelihoods are of different events,
or when a model is separated.

## Degrees of freedom

Each model's `df` is its own
([`logLik()`](https://rdrr.io/r/stats/logLik.html)'s). The hurdle and
zero-inflated fits do not count the parameters of a unit with no
positive count: with unit fixed effects in the participation or zero
equation (`part_fe`, `zero_fe`) or in the count equation (`fe`; dummies
written into a formula are covariates to the fit and are counted), such
a unit's dummy and its count effect sit at their limits, where every one
of its rows has probability one, so it adds nothing to the likelihood.
As fixest leaves such units out of a fixed-effects fit, these parameters
are left out of `df`, so that a hurdle and a zero-inflated model fit to
the same panel are charged alike. The fits record the units
(`part_zero_units` in the hurdles; `zero_units` and `zero_fe_zero_units`
in
[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md);
`fe_zero_units` and `zero_fe_zero_units` in
[`zi_count()`](https://bagozzib.github.io/underdisp/reference/zi_count.md);
`zero_fe_zero_units` in
[`zi_gec()`](https://bagozzib.github.io/underdisp/reference/zi_gec.md)),
and their print methods give the number. In new data
([`predict()`](https://rdrr.io/r/stats/predict.html),
[`score()`](https://bagozzib.github.io/underdisp/reference/score.md),
[`cv_score()`](https://bagozzib.github.io/underdisp/reference/cv_score.md))
the hurdles give such a unit's rows the participation probability at
that limit, zero, as the zero-inflated fits give them a count effect at
its limit. A unit whose every count is positive has its hurdle
participation dummy at the other limit, plus infinity; that dummy is
counted, as the zero-inflated fits count their zero-equation dummy for
the unit, and in new data the hurdles give the unit's rows the
participation probability one (the fits record these units in
`part_one_units`). The zero-inflated fits give a `zero_fe` unit's rows
the structural-zero probability at its limit in the same way: one for a
unit with no positive count, zero for a unit whose every count is
positive (`zero_fe_one_units`). CPB fits on the floor support
(`support = "floor"`) count these parameters, and in new data evaluate
these units at the coefficients the fit stopped at.

The rule is the two-part fits'. The single-equation fixed-effects fits
([`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md),
[`gec_fe()`](https://bagozzib.github.io/underdisp/reference/gec_fe.md),
`count_reg(fe = )`) count every unit's effect, a unit with no positive
count included; set beside a two-part fit to the same panel, such a fit
is charged two AIC points per unit with no positive count that the
two-part fit is not, and `compare_models()` says so in a message.

## See also

[`zi_test()`](https://bagozzib.github.io/underdisp/reference/zi_test.md),
[`compare_dispersion()`](https://bagozzib.github.io/underdisp/reference/compare_dispersion.md)

## Examples

``` r
# \donttest{
set.seed(1); n <- 400; x <- rnorm(n); z <- rnorm(n)
y <- rhurdle_cpb(n, exp(1.2 + 0.5 * x), 0.5, plogis(-0.2 + 0.8 * z))
d <- data.frame(y = y, x = x, z = z)
compare_models(hurdle = hurdle_cpb(y ~ x, data = d, participation = ~ z),
               zi = zi_cpb(y ~ x, data = d, zero = ~ z))
#>        df    logLik      AIC      BIC logscore       rps truncated separated
#> hurdle  5 -552.8141 1115.628 1135.586 1.382035 0.9930062     FALSE     FALSE
#> zi      5 -555.1695 1120.339 1140.296 1.387924 0.9981824     FALSE     FALSE
# }
```
