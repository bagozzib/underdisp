# Average marginal effects

The average, over the rows of `newdata`, of the change in a predicted
quantity when one variable changes: its derivative (a numeric variable,
with `from` and `to` left out) or the contrast between two of its values
(`from` and `to`, for a variable of any type). Every row keeps its own
values of the other variables, and the variable enters through the
model's formula, so interactions and transformations of it (`x:z`,
`I(x^2)`, `log(x)`) are handled.
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md)
instead moves one model-matrix column at a single reference profile.

## Usage

``` r
ame(
  object,
  variable,
  newdata,
  from = NULL,
  to = NULL,
  type = "response",
  eps = NULL,
  level = 0.95,
  offset = NULL,
  weights = NULL,
  B = 0,
  cores = 1L
)
```

## Arguments

- object:

  A fitted model of any class in the package.

- variable:

  The name of a variable of the model, a column of `newdata`.

- newdata:

  The rows to average over, usually the data the model was fit to.

- from, to:

  Two values of `variable` for a contrast. Leave both out for the
  derivative of a numeric variable.

- type:

  The quantity, as the `type` of the class's
  [`predict()`](https://rdrr.io/r/stats/predict.html) method. The
  default `"response"` is the expected count (for a zero-truncated fit,
  such as
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md)'s
  default, the conditional mean E(Y \| Y \>= 1)); the two-part models
  also give their equation-specific quantities.

- eps:

  The step of the central difference for the derivative; by default 1e-4
  times the range of `variable` in `newdata`.

- level:

  Confidence level.

- offset:

  Passed to [`predict()`](https://rdrr.io/r/stats/predict.html).

- weights:

  Optional weights for the average over the rows: a numeric vector or
  the name of a column of `newdata`, such as the frequency weights of a
  fit to tabulated data (the average then equals the one over the
  expanded rows). By default the rows count equally.

- B:

  Refit-bootstrap replicates for a percentile interval; `0` (default)
  uses the fit's own inference (see Details).

- cores:

  Worker processes for the refits; see
  [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

## Value

A data frame of class `"ud_ame"` with `term`, `contrast` (`"dY/dx"`, or
`"to - from"`), `type`, `estimate`, `lower`, `upper`, `method` (the
interval's source) and `n`, the rows averaged (a row with a missing
value is left out).

## Details

The interval follows the rule of
[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md).
With stored bootstrap draws, the average is computed at every draw,
jointly over the coefficients, the dispersion parameter and a
fixed-effects fit's reference unit (the draws of a fixed-effects fit
recentered on its estimate), and the percentile interval is taken.
Otherwise it uses the delta method over the fit's covariance, and
otherwise there is none. A row at a fixed-effects unit's own effect, or
at its own zero-equation dummy, depends on a parameter the draws do not
carry, so an average over such rows, the sample's average partial
effect, has no interval of this kind (`B` below gives it one). Without
the unit column in `newdata` every row enters at the reference unit, the
mean of the identified units' effects on the log scale (a typical,
geometric-mean unit), whose draws the fit carries: that is the average
effect for such a unit over the rows' covariates, not the sample's
average partial effect, and the two can differ substantially when the
unit effects are dispersed.

A fixed-effects count family
([`count_reg()`](https://bagozzib.github.io/underdisp/reference/count_reg.md)
with `fe`) carries its unit dummies in its covariance, so its average
partial effect has a delta-method interval. The slope of a fit on the
floor support (`support = "floor"`) is refused, since its mean jumps at
integer ceilings; a contrast is available.

`B > 0` gives the average a refit-bootstrap interval instead, for any
class: the model is refit on `B` resampled copies of `newdata`, which
must then be the data the model was fit to, resampling whole clusters of
the fit's `cluster`, else whole units of its first fixed-effects column,
else rows (a unit drawn twice becomes two units). Each replicate
averages over its own rows at its own unit effects, and the percentile
interval is taken. Resampling the units treats them as a sample, so this
is an interval for a fixed-effects fit's average partial effect in the
population of units they represent, which carries the variability of
which units are drawn: with few units, or a few with large effects, it
is wide and asymmetric. Each refit is the fit's own call with its data
replaced and its standard errors off (on King's support, every replicate
of a jackknife-corrected fit is corrected without re-running the
validity gate, as in
[`cpb_fe()`](https://bagozzib.github.io/underdisp/reference/cpb_fe.md)'s
own bootstrap, and none is when the gate refused the fit's correction);
the call's other arguments are evaluated where `ame()` is called, as
[`update()`](https://rdrr.io/r/stats/update.html) does, and an offset or
weights given to the fit as one value per row follow the resampled rows.
The run takes about `B` times one fit (see `cores`).

## See also

[`first_difference()`](https://bagozzib.github.io/underdisp/reference/first_difference.md),
[`predict.cpb()`](https://bagozzib.github.io/underdisp/reference/predict.cpb.md)
(`se.fit`, `interval`).

## Examples

``` r
# \donttest{
set.seed(1); n <- 400; x <- rnorm(n); g <- factor(sample(c("a", "b"), n, TRUE))
d <- data.frame(y = rcpb(n, exp(1 + 0.4 * x + 0.3 * (g == "b")), 0.5), x = x, g = g)
fit <- cpb(y ~ x * g, data = d, truncated = FALSE, se = "bootstrap", B = 100)
ame(fit, "x", d)
#>  term contrast     type estimate lower upper             method   n
#>     x    dY/dx response    1.426 1.299 1.526 bootstrap (stored) 400
ame(fit, "g", d, from = "a", to = "b")
#>  term contrast     type estimate lower upper             method   n
#>     g    b - a response    1.108   0.9 1.307 bootstrap (stored) 400

## a fixed-effects fit's average partial effect over its own units
pd <- do.call(rbind, lapply(1:30, function(i) {
  x <- rnorm(8); data.frame(unit = i, x = x, y = rcpb(8, exp(rnorm(1, 0, 0.5) + 0.4 * x), 0.5))
}))
ffit <- cpb_fe(y ~ x, data = pd, fe = "unit")
ame(ffit, "x", pd, B = 100)
#>  term contrast     type estimate lower upper                   method   n
#>     x    dY/dx response    0.407 0.308 0.554 bootstrap (refit, B=100) 240
# }
```
