# First difference for a CPB fit

The effect on a quantity of interest of moving one covariate `from` one
value `to` another, holding the other covariates at their means, with a
percentile interval from the model's bootstrap draws when the fit
carries them.

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

## Value

A `"ud_fd"` data frame – the package-wide first-difference contract
(columns `component`, `from`, `to`, `diff`, `lower`, `upper`, `method`)
– with one row for the requested quantity and a bootstrap percentile
interval on the difference (`method = "bootstrap (stored)"`). Every
`first_difference()` method in the package returns this same shape.

## Examples

``` r
# \donttest{
set.seed(1); x <- rnorm(400)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1); y <- rbinom(400, N, 0.5)
fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 200)
first_difference(fit, "x", from = -1, to = 1, quantity = "mean")
#>  component  from    to  diff lower upper             method
#>       mean 3.027 7.964 4.937 4.567 5.205 bootstrap (stored)
# }
```
