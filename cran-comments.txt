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

* win-builder: R-devel, R-release, R-oldrelease
* macOS builder (mac.R-project.org)
* R-hub containers: gcc-asan, valgrind, nold, rchk
* GitHub Actions: macOS and Windows (R release), Ubuntu (R devel, release and
  oldrel-1), with NOT_CRAN = true so the full test suite runs
* GitHub Actions: clang AddressSanitizer and UndefinedBehaviorSanitizer
  (r-hub clang-asan and clang-ubsan containers)

## R CMD check results

[TO FILL IN from the final runs: 0 errors | 0 warnings | N notes, with one
line per note.]

[TO FILL IN: the result of the sanitizer, valgrind and rchk runs, e.g. "No
sanitizer, valgrind or rchk issues."]

## Reverse dependencies

None.

## Notes

* The package compiles C++ code through Rcpp and has no system requirements
  beyond a C++ compiler.
* Long-running tests are skipped on CRAN with testthat::skip_on_cran(), so the
  check stays within the time limit; they run in continuous integration.
* Tests that compare against suggested packages (DGLMExtPois, ewp,
  gamlss.dist, glmmTMB, rmutil, pscl, sandwich, broom, texreg, DHARMa) skip
  when those packages are not installed.
