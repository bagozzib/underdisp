# Predictions from a CPB fit

Predictions from a CPB fit

## Usage

``` r
# S3 method for class 'cpb'
predict(
  object,
  newdata = NULL,
  type = c("response", "rate", "link", "ceiling", "prob"),
  at = NULL,
  offset = NULL,
  ...
)
```

## Arguments

- object:

  A `"cpb"` object.

- newdata:

  Optional data frame of new covariate profiles; if omitted, the fitted
  data are used.

- type:

  One of `"response"` (the exact mean of the fitted distribution; for a
  zero-truncated fit the conditional mean E(Y \| Y \>= 1)), `"rate"`
  (the CPB rate parameter lambda = exp(x'b), which equals the mean only
  when lambda/(1-alpha) is an integer), `"link"` (the linear predictor,
  log lambda), `"ceiling"` (the implied ceiling lambda/(1-alpha)), or
  `"prob"` (the probability that `Y` equals `at`).

- at:

  For `type = "prob"`, the count value(s) `y` whose probability is
  returned (length 1, or one per row of the prediction data).

- offset:

  Optional offset (log scale) for the `newdata` branch: a numeric vector
  or a column name in `newdata`. For `newdata = NULL` the fit's own
  offset is used.

- ...:

  Unused.

## Value

A numeric vector.

## Examples

``` r
set.seed(1); x <- rnorm(300)
N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "none")
predict(fit, newdata = data.frame(x = 0), type = "ceiling")
#>          
#> 9.543256 
predict(fit, newdata = data.frame(x = c(-1, 1)), type = "prob", at = 5)
#> [1] 0.08812504 0.05952484
```
