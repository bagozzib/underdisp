## underdisp 0.2.0

DRAFT: the check results below are to be filled in from the final checks when
the release is submitted.

This is a feature release with one changed default. The continuous parameter
binomial now follows King's (1989) definition, on the support
0, ..., ceiling(lambda / (1 - alpha)); the definition of earlier versions stays
available through `support = "floor"`, which reproduces their estimates. The
release also adds the hyper-Poisson and exponentially weighted Poisson
families, average marginal effects (`ame()`), held-out ceiling checks
(`ceiling_exceedance()`), and separation checks for the zero-inflation and
participation equations. NEWS.md lists the changes.

## Test environments

* GitHub Actions: ubuntu-latest (R release, devel, oldrel-1), macOS-latest and
  windows-latest (R release), and windows-latest (R devel), with
  NOT_CRAN = true so the full test suite runs
* GitHub Actions: the clang-ubsan and clang-asan sanitizer containers used for
  0.1.2

## R CMD check results

0 errors | 0 warnings | notes: to be filled in.

## Reverse dependencies

None.

## Notes for the reviewer

* The package compiles a small amount of C++ (Rcpp); no system requirements
  beyond a C++11 compiler.
* Long-running test batteries are gated behind testthat::skip_on_cran() so the
  check stays inside the ten-minute budget; they run in continuous integration.
* Tests and examples comparing against glmmTMB, gamlss.dist, rmutil,
  COMPoissonReg, DGLMExtPois, ewp and the Ecdat/wooldridge data sets are
  guarded by requireNamespace() and skip when those Suggests are unavailable.
