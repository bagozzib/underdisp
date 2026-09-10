# Boundary-corrected test for zero-inflation (CPB vs ZI-CPB)

A likelihood-ratio test of whether a count needs a structural-zero
component. The (untruncated) CPB is nested in the zero-inflated CPB at a
structural-zero probability of zero. Because that value lies on the
boundary of the parameter space, the LR statistic follows a
\\\tfrac12\chi^2_0 + \tfrac12\chi^2_1\\ mixture (Self and Liang 1987),
which halves the naive \\\chi^2_1\\ p-value. The test uses an
intercept-only structural-zero probability, so it is a clean
single-parameter boundary test; a covariate-dependent structural-zero
model is better compared with information criteria and proper scores via
[`compare_dispersion()`](https://bagozzib.github.io/underdisp/reference/compare_dispersion.md).

## Usage

``` r
zi_test(object, object2 = NULL, data = NULL, ...)
```

## Arguments

- object:

  Either a model formula — then `data` is required and the nested pair
  is fit internally — or a fitted model (`cpb`, `cpb_fe`, or `zi_cpb`).

- object2:

  The second argument: the data frame when `object` is a formula (so
  `zi_test(y ~ x, mydata)` works), or the second fitted model when
  `object` is a fit. In the model interface exactly one of the two
  models must be a zero-inflated `zi_cpb` and the other its non-inflated
  nest (`cpb` or `cpb_fe`); the two may carry any combination of fixed
  effects and robust/clustered standard errors.

- data:

  A data frame; an alternative to passing it as `object2`.

- ...:

  Passed to
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md) and
  [`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md)
  in the formula interface.

## Value

An object of class `"zi_test"` with the two log-likelihoods, the LR
statistic, and the boundary-corrected p-value.

## Examples

``` r
# \donttest{
set.seed(1); x <- rnorm(500)
y <- rzicpb(500, lambda = exp(1.2 + 0.4 * x), alpha = 0.5, pi = 0.3)
zi_test(y ~ x, data = data.frame(y = y, x = x))
#> Warning: the fitted ceiling reaches max.support (500): alpha is bounded by the guard, not the data. The data may not be underdispersed (compare gec(), which lets the dispersion go either way), or raise max.support.
#> Boundary-corrected LR test for zero-inflation (CPB vs ZI-CPB)
#>   logLik: CPB = -1089.75, ZI-CPB = -855.76
#>   LR = 467.97,  p = 4.426e-104   (0.5 chi^2_0 + 0.5 chi^2_1 mixture)
#>   Reject: a structural-zero component improves fit.
# }
```
