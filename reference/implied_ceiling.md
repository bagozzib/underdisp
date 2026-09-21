# Implied ceiling with a profile-likelihood interval

Returns the observation- (or profile-) specific ceiling
lambda/(1-alpha), with bounds that carry the limits of the interval for
`alpha` (see
[`alpha_confint()`](https://bagozzib.github.io/underdisp/reference/alpha_confint.md);
calibrated when the fit went through
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md))
to the ceiling at the fitted rate. The rate is held at its estimate:
coefficient uncertainty in lambda is not propagated, so the bounds are
not a confidence interval for the ceiling; use
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
with `quantity = "ceiling"` for a fully bootstrapped contrast.

## Usage

``` r
implied_ceiling(object, ...)

# S3 method for class 'cpb'
implied_ceiling(object, newdata = NULL, level = 0.95, ...)

# S3 method for class 'cpb_fe'
implied_ceiling(object, newdata = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `"cpb"` object.

- ...:

  Further arguments passed to methods.

- newdata:

  Optional covariate profiles.

- level:

  Confidence level (default 0.95).

## Value

A data frame with `lambda`, `ceiling`, and `lower`/`upper` bounds.

## Examples

``` r
set.seed(6); x <- rnorm(300)
N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "none")
implied_ceiling(fit, newdata = data.frame(x = c(-1, 0, 1)))
#>     lambda   ceiling     lower     upper
#> 1 3.000664  5.808785  5.497061  6.495959
#> 2 4.457628  8.629225  8.166144  9.650054
#> 3 6.622017 12.819122 12.131193 14.335612
```
