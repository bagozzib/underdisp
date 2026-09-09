# Compare count-model dispersion across the model family

Fits the Poisson and negative binomial defaults alongside the two
underdispersed workhorses—the soft-tail Conway–Maxwell–Poisson (fit
natively, no external dependency) and the hard-ceiling CPB—and returns a
comparison on both fit (log-likelihood, AIC, BIC) and calibration (mean
logarithmic score and ranked probability score, lower is better, plus
the fitted share of zeros against the observed share). Whether the tail
is better described by an accelerated decay (COM-Poisson) or a binding
ceiling (CPB) is a testable question, not an assumption; the scoring
rules are the calibration counterpart to the information criteria. The
CPB's ceiling-exceedance share (observations above the implied ceiling)
is reported as a falsification check on the hard bound.

## Usage

``` r
compare_dispersion(
  formula,
  data,
  max.support = 500,
  hurdle = FALSE,
  zi = FALSE
)
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame.

- max.support:

  Passed to
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- hurdle:

  If `TRUE` and the data contain zeros, also fit and score a
  [`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md).

- zi:

  If `TRUE` and the data contain zeros, also fit and score a
  [`zi_cpb()`](https://bagozzib.github.io/underdisp/reference/zi_cpb.md).
  The mixture EM is slow on large panels, so it is a separate opt-in
  from `hurdle`.

## Value

A list with `table` (a data frame of `df`, `logLik`, `AIC`, `BIC`,
`logscore`, `rps`, and `zero_fit` per model), `obs_zero` (the observed
zero share), `nu` (the COM-Poisson dispersion, `>1` = underdispersion,
`NA` if the COM-Poisson fit failed), `alpha` (the CPB shape parameter),
`ceiling_exceedance` (share of observations above the CPB ceiling), and
`cpb_ok` (`FALSE` if the single-equation CPB could not satisfy its
feasibility constraint on the data, e.g. under heavy zero-inflation with
a wide count range — itself a signal that a hurdle or zero-inflated
model is needed).

## Details

With `hurdle = TRUE` and zeros present, a hurdle-CPB is added, so a
zero-inflated underdispersed process can be compared to the
single-equation models on the same footing.

## See also

[`ud_screen()`](https://bagozzib.github.io/underdisp/reference/ud_screen.md),
[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md),
[`score()`](https://bagozzib.github.io/underdisp/reference/score.md)

## Examples

``` r
set.seed(1); x <- rnorm(120)
N <- pmax(round(exp(1.6 + 0.4 * x) / 0.5), 1); y <- rbinom(120, N, 0.5)
compare_dispersion(y ~ x, data = data.frame(y = y, x = x))$table
#>             df    logLik      AIC      BIC logscore       rps    zero_fit
#> Poisson      2 -238.4130 480.8261 486.4011 1.986775 0.9540835 0.015806078
#> NegBinomial  3 -238.4141 482.8282 491.1906 1.986784 0.9540871 0.015806727
#> COM-Poisson  3 -228.4099 462.8198 471.1822 1.903416 0.9341754 0.003138777
#> CPB          3 -230.4031 466.8062 475.1687 1.920026 0.9375246 0.005861094
#> GEC          3 -230.4234 466.8469 475.2093 1.920195 0.9374796 0.005886297
```
