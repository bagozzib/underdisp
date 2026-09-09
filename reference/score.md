# Proper scoring rules for a fitted count model

Computes the mean logarithmic score and the ranked probability score
(RPS) of a fitted model's predicted distribution against observed
counts. Lower is better for both, and both are proper, so they are the
natural way to compare the calibration of competing count models.

## Usage

``` r
score(fit, newdata = NULL, kmax = NULL)
```

## Arguments

- fit:

  A fitted `underdisp` count model (`cpb`, `cpb_fe`, `hurdle_cpb`,
  `zi_cpb`, `count_reg`, `hurdle_count`, `zi_count`, `gec`, `gec_fe`,
  `hurdle_gec`, or `zi_gec`).

- newdata:

  Optional data frame of held-out observations (including the response).
  When supplied, the fitted model's predicted distribution is evaluated
  on these rows; fixed-effect units unseen in training fall back to the
  mean fixed effect. Offsets are assumed absent on `newdata`.

- kmax:

  Highest count to evaluate; defaults to the maximum observed count in
  the fitting data, raised to the maximum held-out count when `newdata`
  contains larger values. Predicted probabilities are floored at 1e-12
  in the log score, so an observation outside a hard-ceiling model's
  support contributes 27.6 to it.

## Value

A named numeric vector `c(logscore, rps)`.

## Details

By default the scores are computed **in sample** (against the data the
model was fit to). Supply `newdata` to score a fitted model on
**held-out** observations, or use
[`cv_score()`](https://bagozzib.github.io/underdisp/reference/cv_score.md)
for a cross-validated score; an in-sample log score equals `-logLik/n`
and does not penalize model complexity, so for comparing models of
different size the held-out or cross-validated score is the honest
criterion.

## See also

[`cv_score()`](https://bagozzib.github.io/underdisp/reference/cv_score.md)

## Examples

``` r
set.seed(1)
y <- rcpb(400, lambda = 3, alpha = 0.5)
fit <- cpb(y ~ 1, data = data.frame(y = y), truncated = FALSE, se = "none")
score(fit)
#>  logscore       rps 
#> 1.6095598 0.6678615 
```
