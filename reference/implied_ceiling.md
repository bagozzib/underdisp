# Implied ceiling, and the ceiling at the limits of alpha's interval

Returns the observation- (or profile-) specific ceiling
lambda/(1-alpha), and the ceiling at the fitted rate with `alpha` at the
limits of its interval (see
[`alpha_confint()`](https://bagozzib.github.io/underdisp/reference/alpha_confint.md);
calibrated when the fit went through
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md)).
Those two columns hold the rate at its estimate, so they show how far
the dispersion's uncertainty alone moves the ceiling and are not a
confidence interval for it;
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
with `quantity = "ceiling"` propagates the coefficients' uncertainty
too. The ceiling is the CPB's continuous number of trials. No count
above its integer ceiling has positive probability on King's support
(the default), and none above its integer floor on the floor support
(`support = "floor"`); that largest feasible count is returned as
`max_count`. Whether held-out counts respect the bound is tested by
[`ceiling_exceedance()`](https://bagozzib.github.io/underdisp/reference/ceiling_exceedance.md).

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

A data frame with `lambda`, `ceiling`, `max_count` (the largest count
with positive probability), and `ceiling_at_alpha_lower` /
`ceiling_at_alpha_upper`, the ceiling at the fitted rate with `alpha` at
the lower and upper limits of its interval (for `cpb` fits, and hurdles
with a `cpb` intensity).

## Details

The ceiling is the fitted distribution's support end. With the rate and
`alpha` it is fixed by the CPB's variance-mean link, so it is a bound on
the counts only when the CPB is the data's process. On counts with a
common hard bound from another process (binomial counts out of six
trials, say), the fitted ceiling follows the mean and the dispersion
instead: in one simulation it ran from 3 where the mean was low to 12
where it was high, below the true bound in the first case and well above
it in the second.

## Examples

``` r
set.seed(6); x <- rnorm(300)
N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "none")
implied_ceiling(fit, newdata = data.frame(x = c(-1, 0, 1)))
#>     lambda   ceiling max_count ceiling_at_alpha_lower ceiling_at_alpha_upper
#> 1 2.968774  5.593221         6               5.011027               6.610784
#> 2 4.437736  8.360771         9               7.490505               9.881829
#> 3 6.633548 12.497719        13              11.196841              14.771403
```
