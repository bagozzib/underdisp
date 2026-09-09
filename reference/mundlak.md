# Mundlak (correlated random effects) device

Augments a data frame with the unit-level means of the time-varying
numeric covariates in `formula`, and returns the augmented formula and
data. Fitting any of the package's estimators on the result implements
the Mundlak / correlated-random-effects specification: the coefficients
on the original covariates recover the within-unit
(fixed-effects-consistent) effects, while the coefficients on the unit
means capture (and test) the correlation between the covariates and the
unit effect. It is a lighter, between-variation-preserving alternative
to full fixed effects, and it is model-agnostic — the same device feeds
[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md),
[`gec()`](https://bagozzib.github.io/underdisp/reference/gec.md), or
[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md).

## Usage

``` r
mundlak(formula, data, unit, suffix = "_mean")
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame.

- unit:

  Column name identifying the panel unit.

- suffix:

  Suffix for the added unit-mean columns (default `"_mean"`).

## Value

A list with `formula` (the augmented formula), `data` (the augmented
data frame), and `added` (the names of the unit-mean columns).

## See also

[`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md),
[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md)

## Examples

``` r
set.seed(1)
d <- data.frame(y = rpois(200, 3), x = rnorm(200), unit = factor(rep(1:20, each = 10)))
m <- mundlak(y ~ x, d, unit = "unit")
fit <- count_reg(m$formula, m$data, family = "poisson")
```
