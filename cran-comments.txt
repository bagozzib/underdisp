## underdisp 0.1.2

This release exists to correct the undefined behaviour reported for 0.1.1 on
https://www.stats.ox.ac.uk/pub/bdr/M1-SAN/underdisp (mail of 2026-09-28, to be
corrected before 2026-10-19).

## The issue and the fix

The sanitizer reported "outside the range of representable values of type 'int'"
at `cpb_pmf.h:30`, `cpb_fe.cpp:106` (twice) and `cpb_fe.cpp:164`. The value being
converted is a fitted ceiling, `lambda / (1 - alpha)` for the CPB and
`mu / (1 - delta)` for the GEC. Both diverge as the dispersion parameter
approaches its equidispersed limit, and the limit is reachable exactly, because
the parameter is carried on the logit scale and `plogis()` returns 1 once its
argument passes about 40. The ceiling was then `+Inf` and the conversion
undefined.

Every caller already rejected a ceiling above `max.support`, so the fix caps the
value at `max.support + 1` before converting, through one rule in the new
`src/support_cap.h`. I audited every `double` to `int` conversion in `src/` and
routed all six of this kind through it, including two the sanitizer's runs did
not reach (`cpb_fe.cpp:67` and `gec_fe.cpp:104`), so the same defect cannot
resurface from a path your tests happen to exercise later. `tests/testthat/
test-support-cap.R` pins the behaviour at and beyond the limit.

## Test environments

* local Windows 11, R 4.6.1 (`R CMD check --as-cran`).
* GitHub Actions (r-lib/actions): ubuntu-latest (R release, devel, oldrel-1),
  macOS-latest (release), windows-latest (release), with NOT_CRAN = true so the
  full test suite runs.
* win-builder R-devel and R-release.

## R CMD check results

0 errors | 0 warnings | 1 note. The note is the usual incoming-feasibility one,
naming the maintainer and flagging Katz, equidispersion and rootograms in the
Description; these are standard statistical terms.

Check time is unchanged and well inside the ten-minute budget (the 0.1.1 incoming
check ran 363 seconds).

## Notes for the reviewer

* No user-visible interface changed, and no fitted value changed: I compared the
  0.1.1 and 0.1.2 builds on an identical battery of fits and every cross-sectional
  and fixed-effects CPB result, and the calibrated interval, agree bit for bit.
* Long-running test batteries remain gated behind `testthat::skip_on_cran()`.
