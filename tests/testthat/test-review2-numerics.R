## Regression tests for the numerics items of the second panel (referee I; see
## verification/referee_panels_2026-09/panel2/response_p2.md).

test_that("the inner search finds a unit's supremum at ceilings in the hundreds (I-R1)", {
  set.seed(3); Tn <- 10L; x <- rnorm(Tn); y <- rpois(Tn, 400)
  X <- matrix(x, ncol = 1); off <- rep(0, Tn); w <- rep(1, Tn); ms <- max(500L, 10L * max(y))
  r <- underdisp:::cpb_fe_nll_cpp(c(0.5, qlogis(0.5)), X, as.integer(y), off, w, c(0L, Tn), 1L, ms, FALSE, 30L, numeric(0))
  a_cpp <- r[[2]][1]
  unit_ll <- function(a) sum(dcpb(y, exp(a + 0.5 * x), 0.5, log = TRUE))
  o <- optimize(unit_ll, c(a_cpp - 0.02, a_cpp + 0.02), maximum = TRUE, tol = 1e-12)
  expect_gte(-r[[1]], o$objective - 1e-6)
})

test_that("count-family covariances are finite at large fitted means (I-R2)", {
  set.seed(4); n <- 300; x <- rnorm(n); d <- data.frame(y = rpois(n, exp(7 + 0.3 * x)), x = x)
  f <- count_reg(y ~ x, d, family = "gammacount")
  expect_true(all(is.finite(f$se.beta)))
  g <- count_reg(y ~ x, d, family = "gammacount", se = "robust")
  expect_true(all(is.finite(g$se.beta)))
})

test_that("the COM-Poisson mean parameterization is right at strong underdispersion (I-R4, I-R7)", {
  expect_equal(sum(dcompois(0:400, lambda = 1, nu = 5, mu = 200)), 1, tolerance = 1e-8)
  set.seed(1); expect_equal(mean(rcompois(2000, lambda = 1, nu = 5, mu = 200)), 200, tolerance = 0.01)
  expect_true(is.na(dcompois(1, NA, 2))); expect_true(is.na(pcompois(1, Inf, 2)))
})

test_that("the GEC support guard refuses instead of renormalizing (I-R5)", {
  set.seed(3); n <- 200; x <- rnorm(n); d <- data.frame(y = rpois(n, exp(3 + 0.5 * x)), x = x)
  g <- gec(y ~ x, d, se = "none", max.support = 500)
  expect_false(g$support_binding)
  expect_warning(p <- predict(g, newdata = data.frame(x = c(0, 6))), "max.support")
  expect_true(is.finite(p[1]) && is.na(p[2]))
  expect_warning(dd <- dgec(300, 1000, 1), "max.support"); expect_true(is.na(dd))
  expect_equal(dgec(300, 1000, 1, max.support = 5000), dpois(300, 1000), tolerance = 1e-10)
  expect_error(gec(y ~ x, d, se = "none", max.support = 60), "max.support")
})

test_that("a matched-family fit that never leaves the infeasible plateau is not converged (I-R3)", {
  ## the family's mean wall is 1e8: counts near 3e8 leave no finite point for the optimizer
  set.seed(5); n <- 60; x <- rnorm(n); d <- data.frame(y = rpois(n, 3e8), x = x)
  expect_warning(f <- count_reg(y ~ x, d, family = "poisson", se = "none"), "did not converge")
  expect_false(f$converged); expect_equal(f$loglik, -Inf)
})
