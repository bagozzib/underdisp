## underdisp 0.1.1

This is a maintenance and feature release of `underdisp` (first published on CRAN
2026-08-20). It adds four count families, two dispersion diagnostics, frequency
weights, and parallel bootstraps, and it fixes the issues found in an external
review round (see NEWS.md).

## Test environments

* local Windows 11, R 4.6.1 (R CMD check --as-cran): see the results below.
* GitHub Actions (r-lib/actions): ubuntu-latest (R release, devel, oldrel-1),
  macOS-latest (release), windows-latest (release), with NOT_CRAN = true so the
  full test suite runs.
* win-builder R-devel and R-release (submitted 2026-09-10).

## R CMD check results

Local (Windows 11, R 4.6.1, `R CMD check --as-cran`, 2026-09-10):
0 errors | 0 warnings | 1 note. The note is "Files 'README.md' or 'NEWS.md'
cannot be checked without 'pandoc' being installed", a property of the local
check environment. The CRAN-mode test run passes (183 expectations, 0
failures); the full validation batteries run with NOT_CRAN=true in the
package CI on macOS (arm64), Windows, and Linux.

## Notes for the reviewer

* The package compiles a small amount of C++ (Rcpp) for the likelihoods; no
  system requirements beyond a C++11 compiler.
* Long-running test batteries are gated behind `testthat::skip_on_cran()` so
  the check stays inside the ten-minute budget; they run in the package's
  continuous integration.
* Tests and examples that compare against `glmmTMB`, `gamlss.dist`, `rmutil`,
  `COMPoissonReg`, and the `Ecdat`/`wooldridge` data sets are guarded by
  `requireNamespace()` and skip when those Suggests are unavailable.
* Words flagged as possibly misspelled in the Description (Katz, equidispersion,
  rootograms) are standard statistical terms.
