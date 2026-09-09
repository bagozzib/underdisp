# Implied ceiling with a profile-likelihood interval

Returns the observation- (or profile-) specific ceiling
lambda/(1-alpha), with an interval propagating the profile-likelihood
uncertainty in `alpha` at the fitted mean. (Coefficient uncertainty in
lambda is not propagated here; use
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
with `quantity = "ceiling"` for a fully bootstrapped contrast.)

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
#> 1 3.000724  5.808817  5.596138  6.496052
#> 2 4.457693  8.629227  8.313284  9.650141
#> 3 6.622079 12.819057 12.349712 14.335666
```
