## underdisp 0.2.0

This is a feature release. It changes one default: the continuous parameter
binomial (CPB) now follows King's (1989) definition, on the support
0, ..., ceiling(lambda / (1 - alpha)). The definition used by earlier versions
remains available through `support = "floor"`, which reproduces their
estimates. The release also adds the hyper-Poisson and exponentially weighted
Poisson families, average marginal effects (`ame()`), confidence intervals for
predictions, a held-out check of the CPB's ceiling (`ceiling_exceedance()`),
and separation checks for the binary equations of the hurdle and zero-inflated
models. NEWS.md lists the changes.

## Test environments

* local Windows 11, R 4.6.1: R CMD check --as-cran, and again with
  _R_CHECK_DEPENDS_ONLY_=true
* win-builder: R-devel, R 4.6.2, R 4.5.3
* R-hub: macos-arm64 and windows (R-devel); the Linux containers clang-asan,
  clang-ubsan, valgrind, nold, lto, gcc16, donttest, mkl, atlas and rchk
* GitHub Actions: macOS and Windows (R release), Windows (R devel), Ubuntu
  (R devel, release and oldrel-1), with NOT_CRAN = true so the full test
  suite runs
* GitHub Actions: gcc and clang AddressSanitizer and
  UndefinedBehaviorSanitizer (the r-hub gcc-asan, clang-asan and clang-ubsan
  containers)

## R CMD check results

0 errors | 0 warnings | 0 notes on win-builder R-devel, R 4.6.2 and R 4.5.3.

An earlier R 4.5.3 run listed possibly misspelled words in DESCRIPTION: CPB,
Katz, Underdispersed, equidispersion, rootograms, underdispersion. They are
technical terms and a surname, spelled as intended.

Check time on win-builder R-devel was 537 seconds.

No sanitizer or valgrind issues. rchk reports nothing in the package's own
code; its only reports are in Rcpp's headers.

## Reverse dependencies

None.

## Notes

* The package compiles C++ code through Rcpp and has no system requirements
  beyond a C++ compiler.
* Long-running tests are skipped on CRAN with testthat::skip_on_cran(), so the
  check stays within the ten-minute limit; they run in continuous integration.
* Tests that compare against suggested packages (DGLMExtPois, ewp,
  gamlss.dist, glmmTMB, rmutil, pscl, sandwich, broom, texreg, DHARMa) skip
  when those packages are not installed.
