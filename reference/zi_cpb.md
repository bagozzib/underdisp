# Zero-inflated continuous parameter binomial regression

Fits the zero-inflated CPB. The model mixes a structural-zero process (a
logistic model for the probability `pi` that a unit is a structural
zero) with an untruncated CPB for the count, so a zero can arise either
structurally or as a sampling zero from the CPB. Unlike
[`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md),
the likelihood does not factorize, so the parameters are estimated by
direct joint maximum likelihood, seeded from separate fits (a
zero-truncated CPB on the positives and a logit of the zero indicator)
and cross-checked against an expectation-maximization climb of the same
observed-data likelihood, keeping the better optimum. An EM-only path
(`method = "em"`) is retained for comparison.

## Usage

``` r
zi_cpb(
  formula,
  data,
  zero = NULL,
  fe = NULL,
  zero_fe = NULL,
  method = c("ml", "em"),
  se = c("none", "bootstrap"),
  B = 500,
  cluster = NULL,
  max.support = 500,
  maxit = 200,
  tol = 1e-06,
  offset = NULL,
  weights = NULL,
  cores = 1L
)
```

## Arguments

- formula:

  Intensity (count) model formula (`y ~ x`).

- data:

  A data frame.

- zero:

  Optional one-sided formula (`~ z`) for the zero-inflation model;
  defaults to the intensity model's right-hand side.

- fe:

  Optional column name for unit fixed effects in the intensity (count)
  component. Under `method = "ml"` the unit effects enter as plug-in
  offsets from a first-stage fit; under `method = "em"` a classification
  (hard-assignment) step handles the structural zeros.

- zero_fe:

  Optional column name for fixed effects in the zero-inflation equation,
  entered as factor dummies (opt-in).

- method:

  Estimation method: `"ml"` (default) is robust direct joint maximum
  likelihood seeded from separate fits; `"em"` is expectation-
  maximization, retained for comparison but prone to a degenerate
  structural-zero-probability collapse under heavy zero-inflation.

- se:

  Coefficient inference: `"none"` (default) or `"bootstrap"`.

- B:

  Bootstrap resamples when `se = "bootstrap"`.

- cluster:

  Optional cluster for the bootstrap (a column name or vector); under
  fixed effects the units are resampled by default. The mixture EM makes
  the bootstrap costly, so keep `B` modest.

- max.support:

  Maximum support for the CPB pmf.

- maxit:

  Optimizer iteration budget (scaled internally for the joint
  maximization; also the EM iteration cap under `method = "em"`).

- tol:

  Relative convergence tolerance on the observed-data log-likelihood.

- offset:

  Optional exposure offset (log scale) for the count component: a
  numeric vector or the name of a column in `data`, making the count a
  rate model. Supported for `method = "ml"` without `fe` (the
  fixed-effects path's plug-in unit effects have no offset handling yet,
  and errors loudly). The bootstrap carries the offset through every
  replicate.

- weights:

  Optional frequency weights (a numeric vector or a column name); see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).
  Supported for `method = "ml"`.

- cores:

  Worker processes for the bootstrap; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

## Value

An object of class `"zi_cpb"`.

## Details

The inflation equation uses a logit link. For a probit or cloglog
inflation link, use
[`zi_count()`](https://bagozzib.github.io/underdisp/reference/zi_count.md)
(Poisson/NB/COM-Poisson count), whose `link` argument covers the binary
stage.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 800; x <- rnorm(n); z <- rnorm(n)
y <- rzicpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5, pi = plogis(-0.5 + 0.8 * z))
zi_cpb(y ~ x, data = data.frame(y = y, x = x, z = z), zero = ~ z)
#> Zero-Inflated Continuous Parameter Binomial
#> Call:  zi_cpb(formula = y ~ x, data = data.frame(y = y, x = x, z = z),     zero = ~z)
#> N: 800   EM iterations: 236   logLik: -1239.2
#> 
#> Intensity (CPB) coefficients:
#> (Intercept)           x 
#>      1.2967      0.5011 
#> Intensity alpha (shape): 0.5021 
#> 
#> Zero-inflation (logit link) -- positive coefficients raise P(structural zero), i.e. lower the
#> chance of a positive count (the opposite direction from a hurdle participation model):
#> (Intercept)           z 
#>     -0.3674      0.7065 
#> Mean structural-zero probability: 0.418 
# }
```
