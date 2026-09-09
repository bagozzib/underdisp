# Screen a count outcome for underdispersion

Applies the diagnostic sequence developed in Bagozzi (2026): a
*marginal* verdict from the conditional Pearson statistic and a
through-origin score regression (which gives the direction of
dispersion), a negative-binomial-versus-Poisson test for the
overdispersion call, and—most importantly—an *at-risk* verdict on the
positive counts benchmarked against a **zero-truncated Poisson**. The
last step is what separates genuine underdispersion from the artifact of
conditioning on \\Y\>0\\. When the data are not zero-dominated it also
fits the CPB and reports a ceiling-exceedance diagnostic, and it fits a
generalized-Poisson soft-tail comparator.

## Usage

``` r
ud_screen(
  formula,
  data,
  run_cpb = TRUE,
  cpb_max_n = 3000,
  run_gp = TRUE,
  run_comp = TRUE,
  comp_max_par = 30,
  comp_max_n = 5000,
  ztp_threshold = c("calibrated", "bootstrap"),
  ztp_boot_B = 199L,
  cores = 1L,
  digits = 3
)
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame.

- run_cpb:

  Logical; fit the CPB when the data are not zero-dominated (default
  `TRUE`).

- cpb_max_n:

  Skip the CPB fit above this sample size (default 3000).

- run_gp:

  Logical; fit the generalized-Poisson comparator (default `TRUE`).

- run_comp:

  Logical; fit the native COM-Poisson comparator (default `TRUE`).
  Skipped when the mean model carries more than `comp_max_par`
  parameters (the COM-Poisson has no concentrated fixed-effects path, so
  dummy-heavy screens would be slow) or when `n` exceeds `comp_max_n`.

- comp_max_par, comp_max_n:

  Parameter and sample-size gates for the COM-Poisson comparator
  (defaults 30 and 5000).

- ztp_threshold:

  How to set the at-risk test's underdispersion cutoff. `"calibrated"`
  (default) uses the simulation-calibrated rule \\1 -
  2.27/\sqrt{n\_+}\\, whose constant is an estimated standard deviation
  fitted to one calibration grid (no fixed effects, and that grid's rate
  profile); `"bootstrap"` calibrates the cutoff on the data at hand by a
  parametric bootstrap of the fitted zero-truncated Poisson null
  (simulate `ztp_boot_B` at-risk panels at the fitted rates, refit, and
  take the empirical 5\\ anti-conservative in two separable regimes:
  when the mean model's parameter count is a nontrivial share of the
  positive observations (dummy-heavy fixed-effects screens, where the
  null's center drifts with the parameter share), and when the fitted
  rates concentrate at small values (roughly \\\lambda\\ below 2, and
  the more severely the smaller the rates), where the null's spread
  exceeds the fitted constant even without fixed effects (simulated size
  roughly 0.07–0.10 against the nominal 0.05). The bootstrap absorbs
  both departures by construction and is the recommended choice in
  either regime.

- ztp_boot_B:

  Number of parametric-bootstrap replicates (default 199).

- cores:

  Worker processes for the parametric-bootstrap replicates (default 1);
  see [`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md).

- digits:

  Printing precision.

## Value

An object of class `"ud_screen"` with `verdict_marginal`,
`verdict_atrisk`, the conditional and at-risk (ZTP-benchmarked) Pearson
statistics, the NB-vs-Poisson LR test, a log-likelihood comparison,
(when fit) the CPB alpha and ceiling-exceedance share, and the
over-conditioning guard state: `atrisk_skipped` (`TRUE` when the mean
model nearly saturates the positive counts, so the at-risk statistic is
not computed and the printout says why), `sat_ratio` (the fitted
parameter share of the positives), and `overconditioned` (`TRUE` when
that share reaches 0.10, the region where the calibrated threshold is
anti-conservative; the printout then flags the verdict as diagnostic
rather than probative and recommends `ztp_threshold = "bootstrap"`).

## References

King, G. (1989). Variance specification in event count models. *AJPS*
33(3):762-784.

Bagozzi, B. E. (2026). Revisiting underdispersion in political science.
Companion manuscript.

## See also

[`cpb()`](https://bagozzib.github.io/underdisp/reference/cpb.md),
[`compare_dispersion()`](https://bagozzib.github.io/underdisp/reference/compare_dispersion.md)

## Examples

``` r
set.seed(1); x <- rnorm(250)
N <- pmax(round(exp(1.4 + 0.4 * x) / 0.6), 1); y <- rbinom(250, N, 0.6)
ud_screen(y ~ x, data = data.frame(y = y, x = x))
#> 
#> === Underdispersion screen ===
#> Formula: y ~ x 
#> N=250  mean=4.328  max=14  %zero=0.0  (unconditional var/mean=1.07)
#> 
#> MARGINAL verdict: UNDERDISPERSED 
#>    Pearson=0.370  prop.slope=-0.627 (p=<2e-16)
#>    NB vs Poisson LR = 0 (p= 0.5 ; sig => overdispersion)
#> 
#> AT-RISK (y>0) verdict:UNDERDISPERSED  [n_pos=250]
#>    ZTP-Pearson = 0.429  (underdispersed if < 0.856, the calibrated 5% threshold)
#> 
#> Model comparison (log-lik):
#>  Poisson       NB       GP     COMP      CPB 
#> -450.900 -450.901 -450.900 -404.542 -403.106 
#> 
#> COM-Poisson (full data): nu = 2.91  (> 1 = underdispersed, soft tail)
#> 
#> CPB (full data): alpha = 0.433; ceiling-exceedance = 0.0% of obs
#> 
```
