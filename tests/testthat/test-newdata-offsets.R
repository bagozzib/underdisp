## A fit made with an offset carries it to new rows: predict(), score() and ceiling_exceedance() read the fit's
## offset column from newdata (or take `offset =`), refuse rather than assume exposure 1, and first differences
## hold the offset at its mean in every class.

set.seed(5); n <- 400; x <- rnorm(n); z <- rnorm(n); lexpo <- log(runif(n, 1, 4))
d_off <- data.frame(y = rcpb(n, exp(0.8 + 0.3 * x + lexpo), 0.5), x = x, z = z, lexpo = lexpo)

test_that("new rows take the fit's offset column, and a fit with an offset refuses rows without one", {
  f <- cpb(y ~ x, d_off[1:300, ], truncated = FALSE, offset = "lexpo", se = "none")
  expect_equal(predict(f, newdata = d_off[1:300, ], type = "rate"), f$rate)
  expect_identical(ceiling_exceedance(f, d_off[301:400, ])$exceed, 0)          # the model's own draws respect the bound
  expect_error(predict(f, newdata = d_off[1:5, c("y", "x")]), "fitted with an offset")
  expect_equal(predict(f, newdata = d_off[1:5, c("y", "x")], offset = d_off$lexpo[1:5], type = "rate"), f$rate[1:5])
  fv <- cpb(y ~ x, d_off, truncated = FALSE, offset = d_off$lexpo, se = "none")  # a vector offset has no column
  expect_error(predict(fv, newdata = d_off[1:5, ]), "fitted with an offset")
  sc <- score(f, newdata = d_off[301:400, ])
  expect_true(is.finite(sc[["logscore"]]))
  g <- gec(y ~ x, d_off, offset = "lexpo", se = "none")
  expect_equal(predict(g, newdata = d_off, type = "link"), g$linear.predictors)
})

test_that("zi_cpb carries its offset to new rows and to the bound check", {
  zo <- suppressWarnings(zi_cpb(y ~ x, d_off, zero = ~ 1, offset = "lexpo"))
  expect_equal(predict(zo, newdata = d_off, type = "rate"), zo$lambda_full)
  expect_identical(ceiling_exceedance(zo, d_off)$exceed, 0)
  expect_equal(zo$int_fe_ref, mean(d_off$lexpo))
})

test_that("two-part first differences hold the offset at its mean", {
  d <- d_off; d$y0 <- ifelse(runif(n) < 0.3, 0L, d$y)
  zc <- zi_count(y0 ~ x, d, family = "poisson", zero = ~ z, offset = "lexpo")
  fd <- first_difference(zc, "x", 0, 1)
  expect_equal(fd$from[fd$component == "count"], unname(exp(mean(d$lexpo) + zc$coefficients[1])))
  hc <- hurdle_count(y0 ~ x, d, family = "poisson", participation = ~ z, offset = "lexpo")
  expect_equal(hc$int_fe_ref, mean(d$lexpo))
  hg <- hurdle_gec(y0 ~ x, d, participation = ~ z, offset = "lexpo")
  fh <- first_difference(hg, "x", 0, 0)
  expect_equal(fh$from[fh$component == "intensity"],
               predict(hg, newdata = data.frame(x = 0, z = mean(d$z), lexpo = mean(d$lexpo)), type = "intensity"),
               tolerance = 1e-6)
})
