# Incidence rate ratios for a CPB fit

The ratios are \\\exp(\beta)\\, ratios of the CPB's rate parameter
\\\lambda\\. The CPB's mean equals \\\lambda\\ only where
\\\lambda/(1-\alpha)\\ is an integer, and a zero-truncated fit (the
[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md)
default) models \\E(Y \mid Y \ge 1)\\, so a ratio of expected counts can
differ from them;
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
reports the exact mean contrast.

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
set.seed(8); x <- rnorm(200)
N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(200, N, 0.5)
fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 100)
irr(fit)
#>         term equation ratio estimate lower upper             method
#>  (Intercept)    count   IRR    4.447 4.277 4.619 bootstrap (stored)
#>            x    count   IRR    1.504 1.441 1.558 bootstrap (stored)
# }
```
