# First difference for a CPB fit

The effect on a quantity of interest of moving one covariate `from` one
value `to` another, holding the other covariates at their means (or at
the values `at` gives them).

## Usage

``` r
first_difference(object, ...)

# S3 method for class 'cpb'
first_difference(
  object,
  variable,
  from,
  to,
  quantity = c("mean", "ceiling", "prob"),
  y = NULL,
  level = 0.95,
  at = NULL,
  ...
)
```

## Arguments

- object:

  A `"cpb"` object. With `se = "bootstrap"` at fit time the difference
  carries a bootstrap percentile interval; otherwise the point estimate
  is returned with `method = "none"`.

- ...:

  Further arguments passed to methods; unknown arguments error.

- variable:

  Name of a model-matrix column to vary.

- from, to:

  The two values of `variable` to contrast.

- quantity:

  `"mean"` (E(Y); for a zero-truncated fit the conditional mean E(Y \| Y
  \>= 1)), `"ceiling"` (lambda/(1-alpha)), or `"prob"` (P(Y = `y`)).

- y:

  The count value for `quantity = "prob"`.

- level:

  Confidence level (default 0.95).

- at:

  Optional named numeric vector fixing other model-matrix columns at
  given values instead of their means, e.g. `at = c(z = 1, regionB = 0)`
  (a factor's dummies are columns, so a reference level sets them to 0).

## Value

A `"ud_fd"` data frame – the package-wide first-difference contract
(columns `component`, `from`, `to`, `diff`, `lower`, `upper`, `method`)
– with one row for the requested quantity and a bootstrap percentile
interval on the difference (`method = "bootstrap (stored)"`). Every
`first_difference()` method in the package returns this same shape.

## Details

Every `first_difference()` method in the package takes its interval by
one rule. When the fit carries bootstrap draws (`se = "bootstrap"` at
fit time), the quantity is computed at every draw, jointly over the
coefficients and the dispersion parameter (and, for a fixed-effects fit,
the reference unit's effect), and the percentile interval is taken; the
draws of a fixed-effects fit are first recentered on the reported
estimate, because its replicates share the incidental-parameters bias.
Otherwise the interval is the delta method's over the fit's full
covariance, the dispersion parameter included (the matched count
families' analytic, robust or clustered covariance; a hurdle's two
margins, which are independent). A fit with neither gets the point
estimate with `method = "none"`.

## Examples

``` r
# \donttest{
set.seed(1); x <- rnorm(250)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1); y <- rbinom(250, N, 0.5)
fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 100)
first_difference(fit, "x", from = -1, to = 1, quantity = "mean")
#>  component  from    to  diff lower upper             method
#>       mean 3.069 8.394 5.325 5.006 5.615 bootstrap (stored)
# }
```
