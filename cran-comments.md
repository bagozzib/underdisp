## underdisp 0.1.2

This release fixes the undefined behaviour reported for 0.1.1 at
https://www.stats.ox.ac.uk/pub/bdr/M1-SAN/underdisp, which is why it follows
0.1.1 so closely.

The fitted ceiling lambda / (1 - alpha) leaves the range of int as the
dispersion parameter approaches its equidispersed limit, so converting it was
undefined. Ceilings are now capped at max.support + 1 before the conversion.
Callers already treated anything above max.support as out of range, so fitted
values are unchanged. The other conversions of the same kind in src/ get the
same treatment.

## Test environments

* local Windows 11, R 4.6.1 (R CMD check --as-cran)
* r-hub containers clang-ubsan and clang-asan
* GitHub Actions: ubuntu-latest (R release, devel, oldrel-1), macOS-latest,
  windows-latest, with NOT_CRAN = true so the full test suite runs
* win-builder R-devel

## R CMD check results

On win-builder R-devel: 0 errors | 0 warnings | 1 note, the incoming
feasibility note. 0.1.1 was published on 2026-09-28 and the sanitizer report
arrived the same day; this submission answers it. Check time was 361 seconds,
unchanged from 0.1.1.

Two further notes appear only on my own machine: 'V8' is unavailable there, so
math rendering is skipped, and implied_ceiling's example crosses the five-second
elapsed threshold under load (4.6 seconds of CPU time). Neither appears on
win-builder.

The clang-ubsan and clang-asan containers are clean on this version, and the
same check reproduces the reported errors on 0.1.1.

## Notes for the reviewer

* The package compiles a small amount of C++ (Rcpp); no system requirements
  beyond a C++11 compiler.
* Long-running test batteries are gated behind testthat::skip_on_cran() so the
  check stays inside the ten-minute budget; they run in continuous integration.
* Tests and examples comparing against glmmTMB, gamlss.dist, rmutil,
  COMPoissonReg and the Ecdat/wooldridge data sets are guarded by
  requireNamespace() and skip when those Suggests are unavailable.
