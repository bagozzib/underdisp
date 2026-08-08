test_that("ud_screen distinguishes underdispersion from overdispersion", {
  set.seed(2)
  n <- 500; x <- rnorm(n)

  # underdispersed (binomial thinning, var/mean approx 0.6)
  N <- pmax(round(exp(1.4 + 0.4 * x) / 0.6), 1)
  yu <- rbinom(n, N, 0.6)
  su <- ud_screen(yu ~ x, data = data.frame(yu = yu, x = x), run_cpb = FALSE)
  expect_s3_class(su, "ud_screen")
  expect_match(su$verdict_marginal, "UNDER")
  expect_true(su$dd$pearson < 1)

  # overdispersed (negative binomial)
  yo <- MASS::rnegbin(n, mu = exp(1.4 + 0.4 * x), theta = 1)
  so <- ud_screen(yo ~ x, data = data.frame(yo = yo, x = x), run_cpb = FALSE)
  expect_match(so$verdict_marginal, "OVER")
  expect_true(so$dd$pearson > 1)
})
