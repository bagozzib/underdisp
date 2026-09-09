# Hurdle continuous parameter binomial regression

Fits a participation (hurdle) model joined to a zero-truncated CPB for
the positive counts. Because the hurdle log-likelihood factorizes, the
two parts are fit separately: a logistic regression of participation on
all units, and a zero-truncated CPB
([`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md)) on
the positive counts. The result is the natural model for a bounded,
underdispersed participation process — most units at zero, participants
carrying a tight count.

## Usage

``` r
hurdle_cpb(
  formula,
  data,
  participation = NULL,
  fe = NULL,
  part_fe = NULL,
  link = c("logit", "probit", "cloglog"),
  offset = NULL,
  cluster = NULL,
  se = c("none", "bootstrap"),
  B = 500,
  weights = NULL,
  cores = 1L,
  ...
)
```

## Arguments

- formula:

  Intensity model formula (`y ~ x`).

- data:

  A data frame.

- participation:

  Optional one-sided formula (`~ z`) for the participation (hurdle)
  model; defaults to the intensity model's right-hand side.

- fe:

  Optional column name for unit fixed effects in the intensity model,
  absorbed by a concentrated likelihood
  ([`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)).
  This is how the within-unit underdispersion is recovered.

- part_fe:

  Optional column name for fixed effects in the participation model
  (added as factors). Units with no within-unit variation in
  participation are uninformative for a fixed-effects logit.

- link:

  Link for the participation model: `"logit"` (default), `"probit"`, or
  `"cloglog"`. Positive participation coefficients raise the probability
  of a positive count.

- offset:

  Optional offset on the log-mean scale for the intensity (an exposure):
  a numeric vector or the name of a column in `data`.

- cluster:

  Optional cluster identifier (a column name in `data` or a vector) for
  cluster-robust bootstrap inference on the intensity coefficients; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).
  Cluster-robust errors for the participation logit are available
  directly via `sandwich::vcovCL(fit$participation, cluster = ...)`.

- se:

  Inference for the intensity coefficients: `"none"` or `"bootstrap"`.

- B:

  Bootstrap replicates when `se = "bootstrap"`.

- weights:

  Optional frequency weights (a numeric vector or a column name),
  applied to both margins; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- cores:

  Worker processes for the intensity bootstrap; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- ...:

  Passed to
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md) or
  [`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md).

## Value

An object of class `"hurdle_cpb"` with elements `participation` (a
`glm`), `intensity` (a `cpb` or `cpb_fe`), and bookkeeping.

## Examples

``` r
set.seed(1)
n <- 800; x <- rnorm(n); z <- rnorm(n)
y <- rhurdle_cpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5,
                 p = plogis(-0.2 + 0.8 * z))
fit <- hurdle_cpb(y ~ x, data = data.frame(y = y, x = x, z = z),
                  participation = ~ z)
fit
#> Hurdle Continuous Parameter Binomial
#> Call:  hurdle_cpb(formula = y ~ x, data = data.frame(y = y, x = x, z = z),     participation = ~z)
#> Units: 800 (392 participate, 49%)
#> 
#> Participation (logit link) -- positive coefficients raise P(Y > 0), i.e. participation:
#> (Intercept)           z 
#>     -0.0437      0.7652 
#> 
#> Intensity (zero-truncated CPB) coefficients:
#> (Intercept)           x 
#>      1.2831      0.5031 
#> Intensity alpha (shape): 0.5106 
```
