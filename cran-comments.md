## underdisp 0.1.1

This is a maintenance and feature release of `underdisp` (first published on CRAN
2026-08-20). It adds four count families, two dispersion diagnostics, frequency
weights, parallel bootstraps, and an opt-in calibration of the interval for the
CPB dispersion parameter, and it fixes the issues found in two external review
rounds (see NEWS.md).

## Test environments

* local Windows 11, R 4.6.1 (`R CMD check --as-cran`): see the results below.
* GitHub Actions (r-lib/actions): ubuntu-latest (R release, devel, oldrel-1),
  macOS-latest (release), windows-latest (release), with NOT_CRAN = true so the
  full test suite runs.
* win-builder R-devel and R-release.

## R CMD check results

Local (Windows 11, R 4.6.1, `R CMD check --as-cran`): 0 errors | 0 warnings |
1 note. The note is "Skipping checking math rendering: package 'V8' unavailable",
a property of the local check environment.

Overall check time is under eight minutes on that machine. The 0.1.0 submission
was auto-rejected in 2026-08 for exceeding the ten-minute budget, so for this
release the heavy test blocks, examples and vignette chunks were measured
individually and reduced: the CRAN-visible test run is about one minute, the
examples about forty seconds (including `--run-donttest`), and the vignette
rebuild about two minutes.

The CRAN-mode test run passes with 0 failures and 0 errors; the full validation
batteries run with NOT_CRAN = true in the package's continuous integration, which
is green on all five platforms above.

## Notes for the reviewer

* The package compiles a small amount of C++ (Rcpp) for the likelihoods; no
  system requirements beyond a C++11 compiler.
* Long-running test batteries are gated behind `testthat::skip_on_cran()` so the
  check stays inside the ten-minute budget; they run in the package's continuous
  integration.
* Tests and examples that compare against `glmmTMB`, `gamlss.dist`, `rmutil`,
  `COMPoissonReg`, and the `Ecdat`/`wooldridge` data sets are guarded by
  `requireNamespace()` and skip when those Suggests are unavailable.
* Words flagged as possibly misspelled in the Description (Katz, equidispersion,
  rootograms) are standard statistical terms.
