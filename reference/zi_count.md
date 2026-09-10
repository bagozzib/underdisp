# Zero-inflated count regression for the matched count families

Fits a structural-zero mixture whose count component comes from any of
the package's count families (see
[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md),
Families), the count analogue of
[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md).
The count component uses a log link (on the mean for the
mean-parameterized families, on the rate \\\lambda\\ for `"compois"`)
and the inflation probability a link set by `link`. Returned as a
`"zi_count"` object that
[`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md)
and [`score()`](https://bagozzib.github.io/underdisp/reference/score.md)
accept.

## Usage

``` r
zi_count(
  formula,
  data,
  family = c("poisson", "negbin", "compois", "mpcmp", "genpois", "gammacount",
    "doublepois"),
  zero = NULL,
  fe = NULL,
  zero_fe = NULL,
  link = c("logit", "probit", "cloglog"),
  offset = NULL,
  weights = NULL,
  se = c("analytic", "robust", "cluster", "none"),
  cluster = NULL
)
```

## Arguments

- formula:

  Count formula (`y ~ x`).

- data:

  A data frame.

- family:

  The count-component family; see
  [`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md).

- zero:

  Optional one-sided formula for the inflation (structural-zero) model;
  defaults to the count right-hand side.

- fe:

  Optional fixed-effects column for the count equation (factor dummies).

- zero_fe:

  Optional fixed-effects column for the inflation equation (factor
  dummies); zero-equation fixed effects are opt-in.

- link:

  Link for the inflation probability: `"logit"` (default), `"probit"`,
  or `"cloglog"`. Positive inflation coefficients raise the probability
  of a *structural zero*.

- offset:

  Optional offset for the count component, on the linear-predictor (log)
  scale: a numeric vector or the name of a column in `data`.

- weights:

  Optional frequency weights (a numeric vector or a column name); see
  [`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md).

- se, cluster:

  Standard-error type and optional cluster; see
  [`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md).

## Value

An object of class `"zi_count"`; `$fitted.values` is the marginal mean
`(1 - pi) E(Y | count component)` that
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns, and
`$mu` the count component's natural parameter.

## See also

[`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md),
[`hurdle_count()`](https://bagozzib.github.io/underdisp/reference/hurdle_count.md),
[`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md)

## Examples

``` r
set.seed(2); n <- 300; x <- rnorm(n); z <- rnorm(n)
y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.8 * z)) == 1, 0L, rpois(n, exp(1 + 0.3 * x)))
zi_count(y ~ x, data.frame(y = y, x = x, z = z), family = "poisson",
         zero = ~ z, se = "none")
#> Zero-inflated Poisson regression
#> Call:  zi_count(formula = y ~ x, data = data.frame(y = y, x = x, z = z),     family = "poisson", zero = ~z, se = "none")
#> 
#> Count coefficients:
#> (Intercept)           x 
#>      1.0605      0.2989 
#> 
#> Zero-inflation (logit link) -- positive coefficients raise P(structural zero), i.e. lower the
#> chance of a positive count (the opposite direction from a hurdle participation model):
#> (Intercept)           z 
#>     -0.4266      0.5040 
```
