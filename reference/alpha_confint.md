# Profile-likelihood interval for the dispersion parameter alpha

The interval inverts the likelihood-ratio test of `alpha`, re-maximizing
the coefficients at each value. By default the cut is the chi-square one
(a first-order interval); it covers about 0.88 to 0.95 in simulations
from the CPB, with nearly all misses on the upper side, because the
estimate of `alpha` is biased toward zero
([`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md)
explains the mechanism). A fit that went through
[`calibrate_alpha()`](https://bagozzib.github.io/underdisp/reference/calibrate_alpha.md)
gets the interval calibrated by parametric bootstrap instead. Either
interval is model-based: it assumes the CPB and independent
observations.

## Usage

``` r
alpha_confint(object, level = 0.95)
```

## Arguments

- object:

  A `"cpb"` object.

- level:

  Confidence level (default 0.95).

## Value

A length-2 numeric vector (`lower`, `upper`) with attributes `alpha`
(the point estimate), `method` (first-order or calibrated) and
`boundary` (`TRUE` when the profile has not fallen to the cut by
`alpha = 0.005`, so that the lower limit is the parameter bound).

## Examples

``` r
set.seed(7); x <- rnorm(300)
N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "none")
alpha_confint(fit)
#>     lower     upper 
#> 0.4383949 0.5920518 
#> attr(,"alpha")
#>           
#> 0.4846687 
#> attr(,"method")
#> [1] "first-order profile likelihood"
#> attr(,"boundary")
#> [1] FALSE
```
