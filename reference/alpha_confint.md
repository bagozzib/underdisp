# Profile-likelihood interval for the dispersion parameter alpha

Profile-likelihood interval for the dispersion parameter alpha

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
(the point estimate) and `boundary` (`TRUE` if the lower bound is at the
feasibility boundary, i.e. strong underdispersion, where the interval is
one-sided).

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
#> attr(,"boundary")
#> boundary 
#>    FALSE 
```
