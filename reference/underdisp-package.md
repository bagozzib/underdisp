# underdisp: Diagnostics and Models for Underdispersed Count Data

Detect and model underdispersion (conditional variance below the
conditional mean) in count data.

## Screening and diagnostics

[`ud_screen()`](https://bagozzib.github.io/underdisp/reference/ud_screen.md)
(the zero-truncated-Poisson at-risk screen with calibrated or
parametric-bootstrap thresholds),
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md)
(regression-adjusted tests of equidispersion),
[`dispersion_profile()`](https://bagozzib.github.io/underdisp/reference/dispersion_profile.md)
(the conditional variance-to-mean curve against each family's implied
curve),
[`compare_dispersion()`](https://bagozzib.github.io/underdisp/reference/compare_dispersion.md),
[`zi_test()`](https://bagozzib.github.io/underdisp/reference/zi_test.md),
[`rootogram()`](https://bagozzib.github.io/underdisp/reference/rootogram.md),
[`pit_hist()`](https://bagozzib.github.io/underdisp/reference/pit_hist.md).

## Estimators

The hard-ceiling family
[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md) /
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)
and the free-dispersion family
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md) /
[`gec_fe()`](https://bagozzib.github.io/underdisp/reference/gec_fe.md),
with hurdle and zero-inflated forms
([`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md),
[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md),
[`hurdle_gec()`](https://bagozzib.github.io/underdisp/reference/hurdle_gec.md),
[`zi_gec()`](https://bagozzib.github.io/underdisp/reference/zi_gec.md));
the matched count families through one interface,
[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md),
[`hurdle_count()`](https://bagozzib.github.io/underdisp/reference/hurdle_count.md),
[`zi_count()`](https://bagozzib.github.io/underdisp/reference/zi_count.md)
(Poisson, negative binomial, COM-Poisson in the rate and the mean
parameterization, generalized Poisson, gamma-count, double Poisson), all
with offsets, frequency weights, fixed effects, and analytic, robust,
cluster, or bootstrap inference.

## Comparison and quantities of interest

[`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md),
[`score()`](https://bagozzib.github.io/underdisp/reference/score.md),
[`cv_score()`](https://bagozzib.github.io/underdisp/reference/cv_score.md);
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`implied_ceiling()`](https://bagozzib.github.io/underdisp/reference/implied_ceiling.md),
[`irr()`](https://bagozzib.github.io/underdisp/reference/irr.md),
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
(with the extensive/intensive decomposition for the two-part models),
[`confint()`](https://rdrr.io/r/stats/confint.html) for every class,
[`simulate()`](https://rdrr.io/r/stats/simulate.html) for DHARMa
diagnostics, and broom / texreg / modelsummary support.

## See also

Useful links:

- <https://CRAN.R-project.org/package=underdisp>

- <https://github.com/bagozzib/underdisp>

- <https://bagozzib.github.io/underdisp/>

- Report bugs at <https://github.com/bagozzib/underdisp/issues>

## Author

**Maintainer**: Benjamin E. Bagozzi <bagozzib@udel.edu>
([ORCID](https://orcid.org/0000-0002-6233-6453))
