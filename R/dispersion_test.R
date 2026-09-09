## ---------------------------------------------------------------------------
## Regression-adjusted tests of equidispersion.
##
## Two tests, one interface. (1) The likelihood-ratio test of a fitted family's
## dispersion parameter against its Poisson value, on the fitted model's own
## design: the Poisson is the equidispersed member of every family in the
## package (CPB alpha = 1, GEC delta = 1, COM-Poisson nu = 1, generalized
## Poisson lambda = 0, gamma-count alpha = 1, double Poisson theta = 1, negative
## binomial 1/theta = 0). Where the Poisson value is interior to the parameter
## space the LR statistic is chi-square with one degree of freedom under the
## null and a directional alternative uses the signed root of the statistic,
## r = sign(theta_hat - theta_0) * sqrt(LR) ~ N(0, 1). Where the Poisson value
## is on the boundary (the CPB's alpha = 1, the negative binomial's 1/theta = 0)
## the null distribution is the 0.5 chi^2_0 + 0.5 chi^2_1 mixture of Self and
## Liang (1987), whose tail probability equals the one-sided signed-root
## p-value; the test is then necessarily one-sided. (2) The auxiliary-regression
## test of Cameron and Trivedi (1990) for a Poisson fit: regress
## ((y - mu)^2 - y) / mu on mu through the origin; the t statistic of the slope
## tests Var = mu against Var = mu + a * mu^2, with a < 0 underdispersion.
## ---------------------------------------------------------------------------

#' Regression-adjusted test of equidispersion
#'
#' Tests whether a fitted count model's dispersion departs from the Poisson,
#' conditional on the covariates. For a fit from a family with a dispersion
#' parameter (`method = "lr"`, the default), the test is the likelihood-ratio
#' test of that parameter against its Poisson value, with the Poisson refit on
#' the same design (fixed effects included); the alternative can be two-sided
#' or directional (see Details). For a Poisson fit, `method = "auxiliary"` runs
#' the Cameron--Trivedi (1990) auxiliary regression, which needs no
#' alternative family.
#'
#' @details
#' The Poisson is nested in every family at one value of the dispersion
#' parameter: `alpha = 1` (CPB), `delta = 1` (GEC), `nu = 1` (both COM-Poisson
#' parameterizations), `lambda = 0` (generalized Poisson), `alpha = 1`
#' (gamma-count), `theta = 1` (double Poisson), and `1/theta = 0` (negative
#' binomial). When that value is interior to the parameter space the
#' likelihood-ratio statistic is \eqn{\chi^2_1} under the null; a directional
#' alternative (`"under"` or `"over"`) uses the signed root
#' \eqn{r = \pm\sqrt{LR}}, signed toward underdispersion when the estimate lies
#' on the underdispersed side of the Poisson value, which is standard normal
#' under the null, so the one-sided p-value is \eqn{1 - \Phi(r)} for
#' `"under"` and \eqn{\Phi(r)} for `"over"`. When the Poisson value is on the
#' boundary (the CPB, whose `alpha` lives in (0, 1); the negative binomial,
#' whose overdispersion parameter is non-negative) the null distribution of the
#' statistic is the \eqn{\tfrac12\chi^2_0 + \tfrac12\chi^2_1} mixture of Self
#' and Liang (1987), the test is one-sided by construction, and `alternative`
#' is fixed accordingly (`"under"` for the CPB, `"over"` for the negative
#' binomial). For the concentrated fixed-effects fits the Poisson null carries
#' the same unit effects, profiled out in closed form (Poisson) or by a
#' per-unit Newton step (zero-truncated Poisson).
#'
#' @param object A fitted `cpb`, `cpb_fe`, `gec`, `gec_fe`, or `count_reg`
#'   model (for `method = "auxiliary"`, a `count_reg` fit with
#'   `family = "poisson"`).
#' @param alternative `"two.sided"` (default where the null is interior),
#'   `"under"` (underdispersion), or `"over"` (overdispersion).
#' @param method `"lr"` (likelihood ratio against the Poisson refit; default)
#'   or `"auxiliary"` (the Cameron--Trivedi auxiliary regression for a
#'   Poisson fit).
#' @return An object of class `c("dispersion_test", "htest")` with the
#'   statistic, the p-value, the estimated dispersion parameter and its Poisson
#'   value, the two log-likelihoods (LR method), and the alternative.
#' @references Cameron, A. C. and Trivedi, P. K. (1990). Regression-based
#'   tests for overdispersion in the Poisson model. \emph{Journal of
#'   Econometrics}, 46(3), 347-364. Self, S. G. and Liang, K.-Y. (1987).
#'   Asymptotic properties of maximum likelihood estimators and likelihood
#'   ratio tests under nonstandard conditions. \emph{Journal of the American
#'   Statistical Association}, 82(398), 605-610.
#' @examples
#' set.seed(2); n <- 400; x <- rnorm(n)
#' d <- data.frame(y = rgammacount(n, exp(1 + 0.4 * x), alpha = 2), x = x)
#' dispersion_test(count_reg(y ~ x, d, family = "gammacount"), alternative = "under")
#' dispersion_test(count_reg(y ~ x, d, family = "poisson"), method = "auxiliary")
#' @seealso [count_reg()], [dispersion_profile()], [zi_test()]
#' @export
dispersion_test <- function(object, alternative = c("two.sided", "under", "over"),
                            method = c("lr", "auxiliary")) {
  method <- match.arg(method); alt_missing <- missing(alternative)
  alternative <- match.arg(alternative)
  if (method == "auxiliary") return(.disp_auxiliary(object, alternative))
  .disp_lr(object, alternative, alt_missing)
}

## Poisson (or zero-truncated Poisson) log-likelihood maximized on the fit's own
## design, weights, and offset
.disp_poisson_refit <- function(object) {
  if (inherits(object, "count_reg"))
    return(.count_fit(object$X, object$Y, .count_fam("poisson"), isTRUE(object$truncated),
                      offset = object$offset, w = object$weights)$loglik)
  if (inherits(object, c("cpb_fe", "gec_fe"))) return(.disp_fe_null(object))
  .cpb_baseline_loglik(object$X, object$Y,
                       if (is.null(object$offset)) rep(0, length(object$Y)) else object$offset,
                       isTRUE(object$truncated), w = object$weights)
}

## Poisson null for the concentrated fixed-effects fits, with the unit effects
## profiled out: in closed form for the Poisson (a_u = log(sum_t w y) - log(sum_t
## w exp(eta))) and by a vectorized Newton iteration on all units at once for the
## zero-truncated Poisson (score sum_t w (y - m), m = lambda / (1 - e^-lambda))
.disp_fe_null <- function(object) {
  X <- object$X; Y <- object$Y; u <- object$unit
  if (is.null(X) || is.null(u)) stop("The fixed-effects fit does not carry its design; refit with the current package version.")
  w <- .ud_w1(object$weights, length(Y))
  off <- if (is.null(object$offset)) rep(0, length(Y)) else object$offset
  trunc <- isTRUE(object$truncated)
  s <- .ud_colscale(X); Xs <- sweep(X, 2, s, "/")
  sy <- as.numeric(rowsum(w * Y, u)[, 1])
  ll_conc <- function(b) {
    eta <- off + as.numeric(Xs %*% b)
    a <- log(pmax(sy, 1e-300)) - log(as.numeric(rowsum(w * exp(eta), u)[, 1]))
    a[sy <= 0] <- -30
    if (!trunc) return(sum(w * stats::dpois(Y, exp(a[u] + eta), log = TRUE)))
    ## a unit whose counts are all one has its zero-truncated Poisson supremum
    ## at lambda -> 0 (log-likelihood 0): it is held at the floor and skipped
    ones <- as.numeric(rowsum(w * (Y > 1), u)[, 1]) <= 0
    a[ones] <- -8 - max(eta)
    for (it in seq_len(60L)) {
      lam <- exp(a[u] + eta); em <- exp(-lam); m <- lam / (1 - em)
      g <- as.numeric(rowsum(w * (Y - m), u)[, 1])
      h <- as.numeric(rowsum(w * ((1 - em * (1 + lam)) / (1 - em)^2) * lam, u)[, 1])
      step <- pmin(pmax(g / pmax(h, 1e-10), -2), 2)
      step[ones | !is.finite(step)] <- 0
      a <- pmax(a + step, -8 - max(eta))
      if (max(abs(step)) < 1e-9) break
    }
    lam <- exp(a[u] + eta)
    sum(w * (stats::dpois(Y, lam, log = TRUE) - log1p(-exp(-lam))))
  }
  b0 <- tryCatch({ v <- stats::glm.fit(cbind(1, Xs), Y, weights = w, offset = off, family = stats::poisson())$coefficients[-1]
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, ncol(X)))
  o <- stats::optim(b0, function(b) { v <- ll_conc(b); if (is.finite(v)) -v else 1e10 },
                    method = "BFGS", control = list(maxit = 1000, reltol = 1e-12))
  -o$value
}

.disp_lr <- function(object, alternative, alt_missing) {
  boundary <- FALSE; forced <- NULL
  if (inherits(object, c("cpb", "cpb_fe"))) {        # cpb or cpb_fe
    fam_lab <- "CPB"; par_name <- "alpha"; est <- object$alpha; null <- 1
    boundary <- TRUE; forced <- "under"; under_dir <- est < null
    ll1 <- object$loglik
    ll0 <- if (!inherits(object, "cpb_fe") && !is.null(object$loglik.null) && is.finite(object$loglik.null))
             object$loglik.null else .disp_poisson_refit(object)
  } else if (inherits(object, c("gec", "gec_fe"))) { # gec or gec_fe
    fam_lab <- "GEC"; par_name <- "delta"; est <- object$delta; null <- 1; under_dir <- est < null
    ll1 <- object$loglik; ll0 <- .disp_poisson_refit(object)
  } else if (inherits(object, "count_reg")) {
    fam <- .count_fam(object$family)
    if (!fam$nshape) stop("The fitted family (Poisson) has no dispersion parameter; use method = \"auxiliary\" ",
                          "or fit a family with a dispersion parameter.")
    fam_lab <- fam$label; est <- object$theta; null <- fam$null_shape
    par_name <- sub("^(log|atanh)\\((.*)\\)$", "\\2", fam$shape_name)
    ## which side of the Poisson value is underdispersion, per family
    under_dir <- switch(object$family,
      negbin = FALSE, compois = est > null, mpcmp = est > null, genpois = est < null,
      gammacount = est > null, doublepois = est > null)
    if (object$family == "negbin") { boundary <- TRUE; forced <- "over"; est <- 1 / est; null <- 0; par_name <- "1/theta" }
    ll1 <- object$loglik; ll0 <- .disp_poisson_refit(object)
  } else stop("dispersion_test() takes a cpb, cpb_fe, gec, gec_fe, or count_reg fit.")
  if (boundary) {
    if (!alt_missing && alternative != forced)
      warning("the Poisson value is on the boundary of this family's parameter space; the test is one-sided (alternative = \"",
              forced, "\").")
    alternative <- forced
  } else if (alt_missing) alternative <- "two.sided"
  LR <- max(2 * (ll1 - ll0), 0)
  r <- sqrt(LR) * if (isTRUE(under_dir)) 1 else -1     # signed root: positive = underdispersed direction
  p <- if (boundary) 0.5 * stats::pchisq(LR, 1, lower.tail = FALSE)   # Self-Liang mixture
       else switch(alternative,
         two.sided = stats::pchisq(LR, 1, lower.tail = FALSE),
         under = stats::pnorm(r, lower.tail = FALSE),
         over  = stats::pnorm(r))
  structure(list(statistic = c(LR = LR), parameter = c(df = 1), p.value = p,
                 estimate = stats::setNames(est, par_name), null.value = stats::setNames(null, par_name),
                 alternative = switch(alternative, two.sided = "two.sided", under = "underdispersion", over = "overdispersion"),
                 method = paste0("Likelihood-ratio test of equidispersion (", fam_lab, " vs Poisson",
                                 if (boundary) "; boundary null, Self-Liang mixture" else "", ")"),
                 data.name = paste(deparse(object$call$formula), collapse = ""),
                 loglik = c(alternative = ll1, poisson = ll0), boundary = boundary),
            class = c("dispersion_test", "htest"))
}

.disp_auxiliary <- function(object, alternative) {
  if (!inherits(object, "count_reg") || !identical(object$family, "poisson"))
    stop("method = \"auxiliary\" applies to a count_reg(family = \"poisson\") fit.")
  if (isTRUE(object$truncated)) stop("the auxiliary regression is defined for the untruncated Poisson.")
  y <- object$Y; mu <- object$fitted.values; w <- .ud_w1(object$weights, length(y))
  z <- ((y - mu)^2 - y) / mu                      # E[(y-mu)^2 - y] = a * mu^2 under Var = mu + a mu^2
  fit <- stats::lm(z ~ mu + 0, weights = w)
  cf <- summary(fit)$coefficients
  a <- cf[1, 1]; tstat <- cf[1, 3]
  p <- switch(alternative,
    two.sided = 2 * stats::pnorm(-abs(tstat)),
    under = stats::pnorm(tstat),                  # a < 0
    over  = stats::pnorm(tstat, lower.tail = FALSE))
  structure(list(statistic = c(t = tstat), parameter = c(df = fit$df.residual), p.value = p,
                 estimate = c(a = a), null.value = c(a = 0),
                 alternative = switch(alternative, two.sided = "two.sided", under = "underdispersion", over = "overdispersion"),
                 method = "Cameron-Trivedi auxiliary regression test of equidispersion (Var = mu + a mu^2)",
                 data.name = paste(deparse(object$call$formula), collapse = ""), boundary = FALSE),
            class = c("dispersion_test", "htest"))
}

#' @export
print.dispersion_test <- function(x, digits = 4, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("data: ", x$data.name, "\n", sep = "")
  cat(sprintf("%s = %.4f", names(x$statistic), x$statistic),
      if (!is.null(x$loglik)) sprintf("  (logLik: fitted family = %.2f, Poisson = %.2f)", x$loglik[1], x$loglik[2]),
      sprintf(",  p-value = %s\n", format.pval(x$p.value, digits = digits)), sep = "")
  cat("alternative hypothesis: ", x$alternative, "\n", sep = "")
  cat(sprintf("estimate: %s = %.4f  (Poisson value %s)\n", names(x$estimate), x$estimate, format(x$null.value)))
  invisible(x)
}
