# Calibrate the interval for the CPB dispersion parameter

An opt-in step that stores, in the fit, the parametric-bootstrap
distribution of the signed root of the profile likelihood ratio for
`alpha`. Once it is stored,
[`alpha_confint()`](https://bagozzib.github.io/underdisp/reference/alpha_confint.md),
[`confint.cpb()`](https://bagozzib.github.io/underdisp/reference/confint.cpb.md),
[`implied_ceiling()`](https://bagozzib.github.io/underdisp/reference/implied_ceiling.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) report the
calibrated interval instead of the first-order one. It is the remedy for
floor-support fits (`support = "floor"`). On King's support, the
default, the first-order interval covers close to its level, and the
calibration is an optional refinement.

## Usage

``` r
calibrate_alpha(object, B = 199, cores = 1L, seed = NULL, ...)

# S3 method for class 'cpb'
calibrate_alpha(object, B = 199, cores = 1L, seed = NULL, ...)

# S3 method for class 'hurdle_cpb'
calibrate_alpha(object, B = 199, cores = 1L, seed = NULL, ...)

# Default S3 method
calibrate_alpha(object, B = 199, cores = 1L, seed = NULL, ...)
```

## Arguments

- object:

  A pooled
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md) fit,
  or a
  [`hurdle_cpb()`](https://bagozzib.github.io/underdisp/reference/hurdle_cpb.md)
  fit whose intensity is a pooled CPB.

- B:

  Number of bootstrap replicates (default 199).

- cores:

  Worker processes for the refits.

- seed:

  Optional seed for the simulated responses; the caller's random number
  state is restored afterwards, as in
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html).

- ...:

  Unused.

## Value

`object`, with the component `alpha.cal`: a list holding the signed
roots `r` (`NA` for a replicate that failed to fit), `B`, `B_ok` and
`seed`.

## Details

The first-order interval cuts the profile log-likelihood
`qchisq(level, 1) / 2` below its maximum on both sides. On King's
support that cut is close to calibrated, because the probabilities
change continuously with the parameters. The simulations from the CPB
cover nine designs (n from 200 to 1,000, true `alpha` 0.3, 0.5 and 0.7).
There the first-order 95% interval covers between 0.940 and 0.956, 0.947
on average, and more of its misses fall above the upper limit than below
the lower one. With `B = 39` replicates, the calibrated interval covers
0.948 on average in the same simulations. It is about 12% longer, and
its misses are more nearly balanced.

On the floor support the cut is not calibrated. The support
`0, ..., floor(lambda / (1 - alpha))` moves with the parameters, the
log-ceiling coefficients behave like endpoint parameters, and
`alpha`-hat inherits a bias toward zero from them. Across the same nine
designs the first-order 95% interval covers between 0.90 and 0.94, and
most of its misses have the true `alpha` above the upper limit.

The calibration simulates `B` responses from the fitted model (the
fitted rates, `alpha`-hat, the same truncation and offset), refits each
one, and records the signed root
`r* = sign(alpha* - alpha-hat) sqrt(2 (l* - l*_p(alpha-hat)))`. The
interval is then `{alpha : q_lo <= r(alpha) <= q_hi}`, with `q_lo` and
`q_hi` the order statistics of the `r*` at ranks `(B + 1)(1 - level)/2`
from each end (the 5th smallest and 5th largest of 199 at the 95%
level), so the upper limit sits `q_lo^2 / 2` and the lower limit
`q_hi^2 / 2` below the maximum of the profile. On the floor support, in
the same simulations, the calibrated 95% interval covers 0.92 to 0.96,
with its misses balanced between the two sides. On the floor support,
the first-order interval falls short only when the design has many
distinct covariate patterns (the endpoint effect needs many distinct
ceilings): with an intercept only, or a covariate taking ten or fewer
values, it covers at or above its level, and the calibration, which
stays close to the nominal level there too (0.95 with a binary
covariate), is not needed.

The responses are all drawn first, on the calling process, and the
refits use no random numbers, so the result depends on `seed` but not on
`cores`. With `B = 199` the two cuts carry simulation error (a standard
deviation of about 0.2 in `q_lo` and `q_hi`), which moves a limit by a
fraction of its distance from the estimate; `B = 999` reduces it and is
needed for a 99% interval. Each replicate costs about a third of the
original fit.

The interval is model-based. It assumes the CPB and independent
observations, and it is not cluster-robust: when the fit was given
`cluster`, a message says so. Under misspecification `alpha` has no
fixed target (on data with a soft upper tail its estimate rises with the
sample size), and no interval computed under the fitted CPB repairs
that.

Calibration is refused, with the reason, where it has nothing to work
with: a fit whose ceiling reaches `max.support` (the guard, not the
data, bounds `alpha`), weights that are not whole numbers, and a fit
saved by version 0.1.0. Whole-number weights are treated as frequency
weights: a row of weight `w` contributes `w` independent simulated
responses.

## See also

[`alpha_confint()`](https://bagozzib.github.io/underdisp/reference/alpha_confint.md),
[`dispersion_test()`](https://bagozzib.github.io/underdisp/reference/dispersion_test.md)

## Examples

``` r
# \donttest{
set.seed(7); x <- rnorm(100)
d <- data.frame(y = rcpb(100, exp(1.5 + 0.4 * x), 0.5, truncated = TRUE), x = x)
fit <- cpb(y ~ x, d, se = "none")
alpha_confint(fit)                       # first-order
#>     lower     upper 
#> 0.3866204 0.6670207 
#> attr(,"alpha")
#>           
#> 0.4955812 
#> attr(,"method")
#> [1] "first-order profile likelihood"
#> attr(,"boundary")
#> [1] FALSE
## 19 replicates keep the example short and can serve a 90% interval; use 199 or more
fit <- calibrate_alpha(fit, B = 19, seed = 1)
alpha_confint(fit, level = 0.90)         # calibrated
#>     lower     upper 
#> 0.4185265 0.6310002 
#> attr(,"alpha")
#>           
#> 0.4955812 
#> attr(,"method")
#> [1] "profile likelihood, signed root calibrated by parametric bootstrap (19 of 19 replicates)"
#> attr(,"boundary")
#> [1] FALSE
# }
```
