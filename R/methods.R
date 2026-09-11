## ---------------------------------------------------------------------------
## S3 methods for class "cpb"
## ---------------------------------------------------------------------------

#' @export
print.cpb <- function(x, ...) {
  cat("Continuous Parameter Binomial regression",
      if (x$truncated) "(zero-truncated)" else "(untruncated)", "\n")
  cat("Call:  ", deparse(x$call), "\n\n", sep = "")
  cat("Coefficients:\n"); print(round(x$coefficients, 4))
  cat("\nalpha (shape parameter):", round(x$alpha, 4),
      "  median implied bound:", round(stats::median(x$ceiling), 2), "\n")
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  if (isTRUE(x$support_binding)) cat("Note: the fitted ceiling reaches max.support; alpha is bounded by the guard, not the data.\n")
  invisible(x)
}

#' @method coef cpb
#' @export
coef.cpb <- function(object, ...) object$coefficients

#' @method vcov cpb
#' @export
vcov.cpb <- function(object, ...) {
  if (is.null(object$vcov)) return(NULL)                 # no bootstrap requested
  object$vcov
}

#' @method logLik cpb
#' @export
logLik.cpb <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = nobs(object), class = "logLik")

#' @method nobs cpb
#' @export
nobs.cpb <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted

#' @method fitted cpb
#' @export
fitted.cpb <- function(object, ...) object$fitted.values

#' @method residuals cpb
#' @export
residuals.cpb <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  r <- object$Y - object$fitted.values
  if (type == "pearson") {                          # exact variance of the fitted pmf
    v <- .cpb_moments(.cpb_rate(object), object$alpha, isTRUE(object$truncated))$var
    r <- r / sqrt(pmax(v, 1e-12))
  }
  r
}

#' Summarize a CPB fit
#'
#' @param object A `"cpb"` object.
#' @param ... Unused.
#' @return An object of class `"summary.cpb"` with the coefficient table, the
#'   dispersion parameter and its profile-likelihood interval, the implied ceiling,
#'   fit statistics, and the likelihood-ratio test against a (zero-truncated) Poisson.
#' @method summary cpb
#' @export
summary.cpb <- function(object, ...) {
  z <- object$coefficients / object$se.beta
  ctab <- cbind(Estimate = object$coefficients, `Std. Error` = object$se.beta,
                `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z)))
  aci <- .cpb_alpha_profile_ci(object)
  ## the likelihood-ratio statistic against the (zero-truncated) Poisson limit and
  ## its asymptotic boundary-mixture p-value, by the rule dispersion_test()
  ## states: a value the support guard pushed below zero is the boundary value 0,
  ## and the p-value is flagged when the estimated mean parameters push the
  ## asymptotic test's first-order size past the calibration tolerance
  lr  <- .disp_boundary_lr(if (is.na(object$loglik.null)) NA_real_ else -2 * (object$loglik.null - object$loglik),
                           isTRUE(object$support_binding))
  des <- .disp_design(object); size1 <- .disp_first_order_size(des$p, des$n, "under")
  if (is.finite(lr$p) && size1 > .disp_size_tol)
    lr$note <- c(lr$note, sprintf(paste0("the p-value is asymptotic; with %d mean parameters on %s observations its ",
                                         "first-order size at the 5%% level is %.3f, and dispersion_test() calibrates ",
                                         "it by parametric bootstrap"), as.integer(des$p), format(des$n), size1))
  out <- list(call = object$call, truncated = object$truncated, coefficients = ctab,
              alpha = object$alpha, alpha.ci = aci, ceiling = object$ceiling,
              loglik = object$loglik, aic = -2 * object$loglik + 2 * object$df,
              LR = lr$LR, LR.p = lr$p, LR.note = lr$note,
              support_binding = isTRUE(object$support_binding), converged = object$converged,
              n = object$n, nobs = nobs(object), se.type = object$se.type,
              nboot_ok = if (!is.null(object$boot)) attr(object$boot, "nboot_ok") else NA)
  class(out) <- "summary.cpb"
  out
}

#' @export
print.summary.cpb <- function(x, ...) {
  cat("\nContinuous Parameter Binomial regression",
      if (x$truncated) "(zero-truncated)\n" else "(untruncated)\n")
  cat("N =", x$n, if (!is.null(x$nobs) && x$nobs != x$n) paste0(" (weight total ", format(x$nobs), ")"),
      "   inference:", x$se.type, "\n\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
  cat("\nalpha =", round(x$alpha, 4),
      sprintf("  (profile %d%% CI: %.3f to %.3f%s)\n", 95L, x$alpha.ci["lower"],
              x$alpha.ci["upper"], if (x$alpha.ci["boundary"] == 1) ", at feasibility boundary" else ""))
  cat("Implied ceiling lambda/(1-alpha): median", round(stats::median(x$ceiling), 2),
      "  range", paste(round(range(x$ceiling), 2), collapse = " to "), "\n")
  cat("logLik =", round(x$loglik, 2), "   AIC =", round(x$aic, 2), "\n")
  if (!is.na(x$LR))
    cat(sprintf("LR vs %sPoisson (H0: alpha = 1): %.2f, p %s\n",
                if (x$truncated) "ZT-" else "", x$LR, format.pval(x$LR.p)))
  for (nt in x$LR.note) cat("Note: ", nt, ".\n", sep = "")
  if (isTRUE(x$support_binding))
    cat("Note: the fitted ceiling reaches max.support; alpha is bounded by the guard, not the data.\n")
  if (!is.na(x$nboot_ok))
    cat("(", x$nboot_ok, " bootstrap resamples converged)\n", sep = "")
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  invisible(x)
}

#' Confidence intervals for a CPB fit
#'
#' Coefficient intervals use the cold-multistart bootstrap percentile method
#' (validated to nominal coverage); the interval for `alpha` uses the
#' profile-likelihood method, which is reliable except under strong underdispersion,
#' where `alpha` sits at the feasibility boundary and the interval is one-sided.
#'
#' @param object A `"cpb"` object fit with `se = "bootstrap"`.
#' @param parm Optional subset of parameters (coefficient names and/or `"alpha"`).
#' @param level Confidence level (default 0.95).
#' @param ... Unused.
#' @return A matrix of lower/upper bounds.
#' @method confint cpb
#' @export
confint.cpb <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$boot))
    stop("Bootstrap required; refit with se = \"bootstrap\".")
  a  <- (1 - level) / 2
  ok <- object$boot[complete.cases(object$boot), , drop = FALSE]
  cib <- t(apply(ok[, 1:object$p, drop = FALSE], 2, quantile, c(a, 1 - a)))
  rownames(cib) <- names(object$coefficients)
  aci <- .cpb_alpha_profile_ci(object, level = level)
  ci  <- rbind(cib, alpha = aci[c("lower", "upper")])
  colnames(ci) <- c(paste0(format(100 * a), "%"), paste0(format(100 * (1 - a)), "%"))
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
