# Hurdle count regression for the matched count families

Fits a participation model joined to a zero-truncated intensity from any
of the package's count families (see
[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md),
Families), the count analogue of
[`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md).
Returned as a `"hurdle_count"` object that
[`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md)
and [`score()`](https://bagozzib.github.io/underdisp/reference/score.md)
accept, so a hurdle-CPB and a hurdle gamma-count can be compared on one
footing.

## Usage

``` r
hurdle_count(
  formula,
  data,
  family = c("poisson", "negbin", "compois", "mpcmp", "genpois", "gammacount",
    "doublepois", "hyperpois", "ewp3", "ewp2"),
  participation = NULL,
  fe = NULL,
  part_fe = NULL,
  link = c("logit", "probit", "cloglog"),
  offset = NULL,
  weights = NULL,
  se = c("analytic", "robust", "cluster", "none"),
  cluster = NULL
)
```

## Arguments

- formula:

  Intensity formula (`y ~ x`).

- data:

  A data frame.

- family:

  The intensity family; see
  [`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md).

- participation:

  Optional one-sided formula for the participation model; defaults to
  the intensity right-hand side.

- fe, part_fe:

  Optional fixed-effects column names for the intensity and
  participation models (entered as factor dummies). The participation
  slopes carry the dummy logit's incidental-parameters bias in short
  panels (see
  [`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md)).
  The participation dummy of a unit that never participates is not
  counted in the degrees of freedom (see
  [`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md),
  Degrees of freedom). In new data such a unit, and one whose every
  count is positive, get the participation probability at their limits
  (zero and one; with one `part_fe` column), and a unit either equation
  has not seen (the intensity never sees a unit with no positive count)
  enters that equation at its reference unit, as in
  [`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md);
  the participation equation is checked for separation
  (`part_separation`).

- link:

  Link for the participation model: `"logit"` (default), `"probit"`, or
  `"cloglog"`. Positive participation coefficients raise the probability
  of a *positive count* (participation).

- offset:

  Optional offset for the intensity, on the linear-predictor (log)
  scale: a numeric vector or the name of a column in `data`.

- weights:

  Optional frequency weights (a numeric vector or a column name),
  applied to both margins; see
  [`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md).

- se, cluster:

  Standard-error type and optional cluster; see
  [`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md).
  They apply to both margins: the participation equation's standard
  errors are model-based, robust or cluster-robust with the intensity's.

## Value

An object of class `"hurdle_count"`.

## See also

[`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md),
[`zi_count()`](https://bagozzib.github.io/underdisp/reference/zi_count.md),
[`compare_models()`](https://bagozzib.github.io/underdisp/reference/compare_models.md)

## Examples

``` r
set.seed(1); n <- 300; x <- rnorm(n); z <- rnorm(n)
y <- ifelse(rbinom(n, 1, plogis(0.4 + 0.8 * z)) == 1, rpois(n, exp(1 + 0.3 * x)) + 1L, 0L)
hurdle_count(y ~ x, data.frame(y = y, x = x, z = z), family = "poisson",
             participation = ~ z, se = "none")
#> Hurdle Poisson regression
#> Units: 300 (166 participate, 55%)
#> 
#> Participation (logit link) -- positive coefficients raise P(Y > 0), i.e. participation:
#> (Intercept)           z 
#>      0.2493      0.7361 
#> 
#> Intensity (zero-truncated Poisson):
#> (Intercept)           x 
#>      1.3394      0.2183 
```
