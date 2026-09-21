# Compare fitted underdispersed-count models

Compares fitted models from this package — `cpb`, `cpb_fe`,
`hurdle_cpb`, or `zi_cpb`, in any combination of fixed effects and
robust/clustered standard errors — on information criteria and, where a
predicted distribution is available, on proper scores. Unlike
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

  Two or more fitted models (`cpb`, `cpb_fe`, `hurdle_cpb`, `zi_cpb`),
  optionally named.

## Value

A data frame with `df`, `logLik`, `AIC`, `BIC`, and (where the predicted
distribution is available) `logscore` and `rps`, one row per model,
ordered by AIC.

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
#>        df    logLik      AIC      BIC logscore       rps
#> hurdle  5 -547.3110 1104.622 1124.579 1.368277 0.9901023
#> zi      5 -550.3732 1110.746 1130.704 1.375933 0.9951678
# }
```
