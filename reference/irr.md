# Incidence rate ratios for a CPB fit

Incidence rate ratios for a CPB fit

## Usage

``` r
irr(object, ...)

# S3 method for class 'cpb'
irr(object, level = 0.95, ...)
```

## Arguments

- object:

  A `"cpb"` object. With `se = "bootstrap"` at fit time the ratios carry
  bootstrap percentile intervals; otherwise point estimates are returned
  with `method = "none"`, as for every other `irr()` method.

- ...:

  Further arguments passed to methods.

- level:

  Confidence level (default 0.95).

## Value

A `"ud_irr"` data frame – the package-wide rate-ratio contract (columns
`term`, `equation`, `ratio`, `estimate`, `lower`, `upper`, `method`)
shared by every `irr()` method; here with bootstrap percentile intervals
from the stored draws.

## Examples

``` r
# \donttest{
set.seed(8); x <- rnorm(300)
N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 100)
irr(fit)
#>         term equation ratio estimate lower upper             method
#>  (Intercept)    count   IRR    4.486 4.265 4.655 bootstrap (stored)
#>            x    count   IRR    1.538 1.475 1.607 bootstrap (stored)
# }
```
