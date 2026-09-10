## Diagnose the macOS-only failure of gec_fe(): every start is reported infeasible.
## Prints the toolchain, the outcome of the help-page example, and, at each of
## gec_fe()'s two dispersion starts, which units the concentrated objective marks
## infeasible, next to the pooled Katz log-likelihood of the same observations
## (src/gec_ll.cpp, which works on macOS) for comparison.
suppressPackageStartupMessages(library(underdisp))
cat(R.version.string, "|", Sys.info()[["machine"]], "| underdisp", as.character(utils::packageVersion("underdisp")), "\n")

## the help-page example
set.seed(5)
d <- do.call(rbind, lapply(1:30, function(i) {
  x <- rnorm(10); N <- pmax(round(exp(1 + rnorm(1, 0, 0.4) + 0.3 * x) / 0.5), 1)
  data.frame(unit = i, x = x, y = rbinom(10, N, 0.5))
}))
fit <- tryCatch(gec_fe(y ~ x, data = d, fe = "unit"), error = function(e) e)
if (inherits(fit, "error")) cat("gec_fe: ERROR:", conditionMessage(fit), "\n") else
  cat(sprintf("gec_fe: OK  delta = %.5f  beta = %.5f  loglik = %.5f\n", fit$delta, coef(fit)[["x"]], fit$loglik))
cf <- tryCatch(cpb_fe(y ~ x, data = d, fe = "unit", truncated = FALSE), error = function(e) e)
if (inherits(cf, "error")) cat("cpb_fe: ERROR:", conditionMessage(cf), "\n") else
  cat(sprintf("cpb_fe: OK  alpha = %.5f  loglik = %.5f\n", cf$alpha, cf$loglik))
pg <- tryCatch(gec(y ~ x + factor(unit), data = d), error = function(e) e)
if (inherits(pg, "error")) cat("pooled gec with dummies: ERROR:", conditionMessage(pg), "\n") else
  cat(sprintf("pooled gec with dummies: OK  delta = %.5f  loglik = %.5f\n", pg$delta, pg$loglik))

## gec_fe()'s internals, as in R/gec.R
d <- d[order(d$unit), ]
Y <- as.integer(d$y); X <- matrix(d$x, ncol = 1, dimnames = list(NULL, "x"))
w <- rep(1, length(Y)); off <- rep(0, length(Y))
s <- underdisp:::.ud_colscale(X, w); Xs <- sweep(X, 2, s, "/")
uf <- factor(d$unit); nu <- nlevels(uf)
ustart <- as.integer(c(0, cumsum(tabulate(as.integer(uf), nu))))
ms <- as.integer(max(500L, 10L * max(Y)))
bs <- stats::glm.fit(cbind(1, Xs), Y, family = stats::poisson())$coefficients[-1]
cat(sprintf("\ncolumn scale %.6f  Poisson slope start %.6f  max.support %d\n", s, bs, ms))

nll <- underdisp:::gec_fe_nll_cpp
for (ld in c(-0.5, 0.2)) {
  par <- c(bs, ld)
  r <- nll(par, Xs, Y, off, w, ustart, nu, ms, 30L, numeric(0))
  cat(sprintf("\nlog delta = %.2f (delta = %.4f): concentrated nll = %s\n", ld, exp(ld), format(r[[1]], digits = 12)))
  bad <- 0L
  for (u in seq_len(nu)) {
    idx <- (ustart[u] + 1L):ustart[u + 1L]
    one <- function(it, awarm = numeric(0))
      nll(par, Xs[idx, , drop = FALSE], Y[idx], off[idx], w[idx], c(0L, length(idx)), 1L, ms, it, awarm)
    ru <- one(30L)
    if (ru[[1]] < 1e9) next
    bad <- bad + 1L
    if (bad > 3L) next
    o <- as.numeric(Xs[idx, 1] * bs)
    apois <- log(sum(Y[idx]) / sum(exp(o)))
    cat(sprintf("  unit %d infeasible (cold, inner_it 30); y = %s; Poisson intercept %.4f\n",
                u, paste(Y[idx], collapse = " "), apois))
    for (it in c(1L, 5L)) cat(sprintf("    cold, inner_it %d: nll %s\n", it, format(one(it)[[1]], digits = 12)))
    rw <- one(30L, awarm = apois)
    cat(sprintf("    warm start at the Poisson intercept: nll %s, intercept %s\n",
                format(rw[[1]], digits = 12), format(rw[[2]], digits = 10)))
    grid <- seq(apois - 1.5, apois + 1.5, by = 0.25)
    lp <- vapply(grid, function(a) sum(underdisp:::gec_lp0_cpp(c(a, 1, ld), cbind(1, o), Y[idx], off[idx], ms)[, 1]), 0)
    cat("    pooled Katz log-likelihood of the unit across intercepts:\n")
    print(round(rbind(intercept = grid, loglik = lp), 4))
  }
  cat(sprintf("  units infeasible at this start: %d of %d\n", bad, nu))
}
