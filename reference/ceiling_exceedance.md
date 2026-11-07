# Held-out ceiling exceedance of a CPB fit

Whether counts respect the CPB's ceiling can be checked, but not in
sample: its likelihood is minus infinity for any count above the bound,
so every feasible fit respects the bound on the data it was fit to. The
check is one-sided: it detects a ceiling that is too low, where held-out
counts reach it, and has no power against one that is too high (a fitted
ceiling above a true bound is never exceeded), nor where counts stay
well below the ceiling (see
[`implied_ceiling()`](https://bagozzib.github.io/underdisp/reference/implied_ceiling.md)).
`ceiling_exceedance()` takes held-out observations and reports how many
of their counts exceed the largest count the fitted model allows at
their covariates (`max_count` of
[`implied_ceiling()`](https://bagozzib.github.io/underdisp/reference/implied_ceiling.md));
for a two-part model the bound is the count component's, and a
fixed-effects unit seen in training enters at its own effect (an unseen
one at the mean effect), as in
[`score()`](https://bagozzib.github.io/underdisp/reference/score.md). A
breached ceiling also zeroes the held-out likelihood, which
[`score()`](https://bagozzib.github.io/underdisp/reference/score.md)
records with its 1e-12 floor.

## Usage

``` r
ceiling_exceedance(fit, newdata, weights = NULL, offset = NULL)
```

## Arguments

- fit:

  A `"cpb"`, `"cpb_fe"`, `"hurdle_cpb"`, or `"zi_cpb"` object.

- newdata:

  Held-out observations, including the response.

- weights:

  Optional frequency weights for the rows of `newdata`: a column name of
  `newdata` or a numeric vector with one value per row.

- offset:

  Optional offset (log scale) of the count for the `newdata` rows; by
  default a fit made with an offset reads it from the same column of
  `newdata` (see
  [`predict.cpb()`](https://bagozzib.github.io/underdisp/reference/predict.cpb.md)).

## Value

A one-row data frame with `n` (held-out rows, weighted), `exceed` (rows
above the bound, weighted), and `share`.

## See also

[`implied_ceiling()`](https://bagozzib.github.io/underdisp/reference/implied_ceiling.md),
[`score()`](https://bagozzib.github.io/underdisp/reference/score.md)

## Examples

``` r
set.seed(7); x <- rnorm(400); y <- rcpb(400, exp(1.4 + 0.3 * x), 0.5)
d <- data.frame(y = y, x = x)
fit <- cpb(y ~ x, d[1:300, ], truncated = FALSE, se = "none")
ceiling_exceedance(fit, d[301:400, ])
#>     n exceed share
#> 1 100      0     0
```
