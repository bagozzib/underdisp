# Cross-validated proper scores for a count model

K-fold cross-validated log score and RPS: the data are split into `k`
folds, the model is refit on each training set via `fitfun`, and its
predicted distribution is scored on the held-out fold. Because the score
is genuinely out of sample it penalizes over-parameterization, so it is
the honest criterion for comparing models that differ in the number of
parameters (for example a hurdle with per-unit participation fixed
effects against a zero-inflated model). Per-observation held-out scores
are returned so two models fit on the *same* folds can be compared with
a paired, cluster-robust test.

## Usage

``` r
cv_score(fitfun, data, k = 5, kmax = NULL, folds = NULL, cores = 1L)
```

## Arguments

- fitfun:

  A function of one argument (a training data frame) returning a fitted
  `underdisp` model, e.g.
  `function(d) hurdle_cpb(y ~ x, d, participation = ~ z, fe = "unit", se = "none")`.

- data:

  The full data frame.

- k:

  Number of folds (default 5).

- kmax:

  Highest count to evaluate; if `NULL`, taken from a fit on the full
  data (one extra fit).

- folds:

  Optional integer vector of length `nrow(data)` giving a fixed fold
  assignment (so two models can be scored on identical folds). If
  `NULL`, folds are drawn at random.

- cores:

  Worker processes for the fold refits (default 1); see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

## Value

An object of class `"cv_score"`: a list with the mean held-out
`logscore` and `rps`, the per-observation vectors `logscore_i`/`rps_i`,
and the `folds` used.

## See also

[`score()`](https://bagozzib.github.io/underdisp/reference/score.md)

## Examples

``` r
# \donttest{
set.seed(1)
d <- data.frame(y = rcpb(500, lambda = 3, alpha = 0.5))
cv <- cv_score(function(tr) cpb(y ~ 1, tr, truncated = FALSE, se = "none"), d, k = 5)
cv$logscore
#> [1] 1.613497
# }
```
