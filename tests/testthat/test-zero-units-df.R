## A unit with no positive count adds nothing to a two-part fit's likelihood, and none of its parameters counts in
## the degrees of freedom (?compare_models, Degrees of freedom): on a panel with such units each class's df equals its
## df on the panel without them, whichever equation carries the unit dummies. The panel's first unit (the reference
## level of the dummies) is one of them.
zero_df_panel <- function(seed = 11, U = 12, TT = 10, zu = c("u01", "u05", "u09")) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(U), function(i) {
    x <- rnorm(TT); a <- rnorm(1, 0.3, 0.4)
    part <- rbinom(TT, 1, plogis(0.4 + 0.5 * x))
    data.frame(unit = sprintf("u%02d", i), t = seq_len(TT), x = x,
               y = part * (1L + rpois(TT, exp(0.5 + a + 0.3 * x))))
  }))
  d$y[d$unit %in% zu] <- 0L
  list(full = d, sub = d[!d$unit %in% zu, ], zu = zu)
}

test_that("the hurdles do not count the participation dummies of units that never participate", {
  p <- zero_df_panel()
  h <- hurdle_cpb(y ~ x, p$full, participation = ~ x, fe = "unit", part_fe = "unit")
  hs <- hurdle_cpb(y ~ x, p$sub, participation = ~ x, fe = "unit", part_fe = "unit")
  expect_identical(h$df, hs$df)
  expect_setequal(h$part_zero_units, p$zu)
  expect_equal(as.numeric(logLik(h)), as.numeric(logLik(hs)), tolerance = 1e-6)
  expect_identical(attr(logLik(h), "df"), h$df)
  expect_output(print(h), "3 units have no positive count. Their participation dummies are at the limit")
  cm <- suppressMessages(suppressWarnings(compare_models(h = h, z = zi_cpb(y ~ x, p$full, zero = ~ x, fe = "unit",
                                                                         zero_fe = "unit"))))
  expect_equal(cm["h", "df"], h$df)
  ## the floor support counts them, as before 0.2.0
  hf <- hurdle_cpb(y ~ x, p$full, participation = ~ x, fe = "unit", part_fe = "unit", support = "floor")
  expect_identical(hf$df, length(coef(hf$participation)) + hf$intensity$df)
  expect_length(hf$part_zero_units, 0L)
  ## a hurdle and its count-family twin share the participation equation, so they are charged alike
  hc <- hurdle_count(y ~ x, p$full, family = "poisson", participation = ~ x, fe = "unit", part_fe = "unit",
                     se = "none")
  hcs <- hurdle_count(y ~ x, p$sub, family = "poisson", participation = ~ x, fe = "unit", part_fe = "unit",
                      se = "none")
  expect_identical(hc$df, hcs$df)
  expect_identical(hc$df - hc$intensity$df, h$df - h$intensity$df)
  expect_setequal(hc$part_zero_units, p$zu)
})

test_that("hurdle_gec, zi_gec, zi_count and zi_cpb(zero_fe) do not count them either", {
  skip_on_cran()
  p <- zero_df_panel()
  g <- hurdle_gec(y ~ x, p$full, participation = ~ x, part_fe = "unit")
  gs <- hurdle_gec(y ~ x, p$sub, participation = ~ x, part_fe = "unit")
  expect_identical(g$df, gs$df); expect_setequal(g$part_zero_units, p$zu)
  zg <- suppressWarnings(zi_gec(y ~ x, p$full, zero = ~ 1, zero_fe = "unit"))
  zgs <- suppressWarnings(zi_gec(y ~ x, p$sub, zero = ~ 1, zero_fe = "unit"))
  expect_identical(zg$df, zgs$df); expect_setequal(zg$zero_fe_zero_units, p$zu)
  zc <- suppressWarnings(zi_count(y ~ x, p$full, family = "poisson", zero = ~ 1, fe = "unit", zero_fe = "unit",
                                  se = "none"))
  zcs <- suppressWarnings(zi_count(y ~ x, p$sub, family = "poisson", zero = ~ 1, fe = "unit", zero_fe = "unit",
                                   se = "none"))
  expect_identical(zc$df, zcs$df)
  expect_setequal(zc$fe_zero_units, p$zu); expect_setequal(zc$zero_fe_zero_units, p$zu)
  expect_output(print(zc), "Their count effects and zero-equation dummies are at the limit")
  zc1 <- suppressWarnings(zi_count(y ~ x, p$full, family = "poisson", zero = ~ 1, fe = "unit", se = "none"))
  zc1s <- suppressWarnings(zi_count(y ~ x, p$sub, family = "poisson", zero = ~ 1, fe = "unit", se = "none"))
  expect_identical(zc1$df, zc1s$df)
  ## zi_cpb with zero-equation dummies only: on King's support not counted, on the floor support counted
  z <- suppressWarnings(zi_cpb(y ~ x, p$full, zero = ~ 1, zero_fe = "unit"))
  zs <- suppressWarnings(zi_cpb(y ~ x, p$sub, zero = ~ 1, zero_fe = "unit"))
  expect_identical(z$df, zs$df); expect_setequal(z$zero_fe_zero_units, p$zu)
  zf <- suppressWarnings(zi_cpb(y ~ x, p$full, zero = ~ 1, zero_fe = "unit", support = "floor"))
  zfs <- suppressWarnings(zi_cpb(y ~ x, p$sub, zero = ~ 1, zero_fe = "unit", support = "floor"))
  expect_identical(zf$df, zfs$df + length(p$zu)); expect_length(zf$zero_fe_zero_units, 0L)
  ## with count fixed effects zi_cpb leaves the units out, and the matched hurdle is charged alike
  zfe <- suppressMessages(suppressWarnings(zi_cpb(y ~ x, p$full, zero = ~ x, fe = "unit", zero_fe = "unit")))
  zfes <- suppressMessages(suppressWarnings(zi_cpb(y ~ x, p$sub, zero = ~ x, fe = "unit", zero_fe = "unit")))
  expect_identical(zfe$df, zfes$df)
})
