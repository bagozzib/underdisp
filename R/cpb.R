## ---------------------------------------------------------------------------
## Continuous Parameter Binomial (CPB) regression: estimator + internal helpers
## ---------------------------------------------------------------------------

## Factor levels of a model frame (for safe prediction on newdata), avoiding the
## non-exported stats:::.getXlevels.
.cpb_xlevels <- function(mf) {
  fac <- vapply(mf, is.factor, logical(1))
  if (!any(fac)) return(NULL)
  lapply(mf[fac], levels)
}

## Negative log-likelihood via the compiled C++ backend (see src/cpb_ll.cpp).
## `offset` is a length-n vector added to the log-mean (an exposure); 0 gives the
## plain model. With frequency weights `w` the weighted form is used.
.cpb_nll <- function(params, X, Y, offset, max.support, truncated, w = NULL) {
  if (is.null(w)) cpb_nll_cpp(params, X, as.integer(Y), offset, as.integer(max.support), truncated)
  else cpb_wnll_cpp(params, X, as.integer(Y), w, offset, as.integer(max.support), truncated)
}

## Single Nelder-Mead fit from one starting value on the (column-scaled) design.
.cpb_fit_one <- function(Xs, Y, start, offset, max.support, truncated, maxit, reltol, w = NULL) {
  optim(start, .cpb_nll, X = Xs, Y = Y, offset = offset, max.support = max.support, truncated = truncated,
        w = w, method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol))
}

## Multi-start fit. The CPB surface is piecewise smooth: every time a parameter
## move shifts an observation's support floor, the normalizing constant jumps,
## so a simplex can settle on a local tooth. The fit therefore multi-starts over
## the dispersion parameter, polishes the incumbent by Nelder-Mead restarts, and
## then restarts from deterministic perturbations of the incumbent whose
## intercept is raised to keep every count below its ceiling, keeping the best.
## The design columns are scaled to unit standard deviation for the optimizer
## and the coefficients mapped back, so the fit is invariant to the covariates'
## units. Returns the lowest-objective fit with `par` on the original scale.
.cpb_fit_ms <- function(X, Y, offset, max.support, truncated, alpha.starts, maxit, reltol, w = NULL,
                        restarts = 4L) {
  s <- .ud_colscale(X); Xs <- sweep(X, 2, s, "/"); p <- ncol(X)
  wv <- .ud_w1(w, nrow(X))
  fn <- function(par) .cpb_nll(par, Xs, Y, offset, max.support, truncated, w)
  bs <- tryCatch({ v <- glm.fit(Xs, Y, weights = wv, offset = offset, family = poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  cand <- lapply(alpha.starts, function(a0)
    tryCatch(.cpb_fit_one(Xs, Y, c(bs, qlogis(a0)), offset, max.support, truncated, maxit, reltol, w),
             error = function(e) NULL))
  cand <- cand[!vapply(cand, is.null, logical(1))]
  if (!length(cand)) return(NULL)
  best <- cand[[which.min(vapply(cand, function(f) f$value, numeric(1)))]]
  best <- .ud_nm_polish(best, fn, maxit, reltol)
  if (restarts > 0L && best$value < 1e9) {
    icol <- match("(Intercept)", colnames(X)); pos <- Y > 0
    pats <- list(rep(0.15, p + 1L), rep(-0.15, p + 1L), 0.15 * (-1)^seq_len(p + 1L), -0.15 * (-1)^seq_len(p + 1L))
    for (dlt in pats[seq_len(min(restarts, 4L))]) {
      st <- best$par + dlt
      if (!is.na(icol) && any(pos)) {                    # feasibility repair: raise the intercept
        a <- plogis(st[p + 1L]); eta <- offset + as.numeric(Xs %*% st[seq_len(p)])
        need <- max(log(Y[pos] * (1 - a)) - eta[pos]) + 1e-6
        if (need > 0) st[icol] <- st[icol] + need
      }
      o2 <- tryCatch(.cpb_fit_one(Xs, Y, st, offset, max.support, truncated, maxit, reltol, w),
                     error = function(e) NULL)
      if (!is.null(o2) && o2$value < best$value - 1e-8) best <- .ud_nm_polish(o2, fn, maxit, reltol)
    }
  }
  best$par[seq_len(p)] <- best$par[seq_len(p)] / s
  best
}

## Baseline log-likelihood: (zero-truncated) Poisson MLE, the correct H0: alpha = 1.
.cpb_baseline_loglik <- function(X, Y, offset, truncated, w = NULL) {
  wv <- .ud_w1(w, nrow(X)); s <- .ud_colscale(X); Xs <- sweep(X, 2, s, "/")
  nll <- if (truncated)
    function(b) { l <- as.numeric(exp(offset + Xs %*% b))
      v <- sum(wv * (dpois(Y, l, log = TRUE) - log(1 - exp(-l)))); if (!is.finite(v)) 1e10 else -v }
  else
    function(b) { l <- as.numeric(exp(offset + Xs %*% b))
      v <- sum(wv * dpois(Y, l, log = TRUE)); if (!is.finite(v)) 1e10 else -v }
  bs <- tryCatch({ v <- glm.fit(Xs, Y, weights = wv, offset = offset, family = poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, ncol(X)))
  o <- tryCatch(optim(bs, nll, method = "BFGS"), error = function(e) NULL)
  if (is.null(o)) NA_real_ else -o$value
}

#' Fit a Continuous Parameter Binomial (CPB) regression
#'
#' Fits the underdispersed continuous parameter binomial model of King (1989), in
#' which the conditional variance is a fraction of the conditional mean,
#' \eqn{\mathrm{Var}(Y\mid x) = \alpha\,\mathrm{E}(Y\mid x)} with
#' \eqn{0 < \alpha < 1}, and each observation has an endogenous ceiling
#' \eqn{\lambda_i/(1-\alpha)}. A zero-truncated variant (the default) conditions on
#' \eqn{Y \ge 1}, appropriate when underdispersion lives among the positive counts of
#' an otherwise zero-inflated outcome.
#'
#' The mean is modelled log-linearly, \eqn{\lambda_i = \exp(x_i'\beta)}. Because the
#' support \eqn{0, \ldots, \lfloor \lambda_i/(1-\alpha) \rfloor} depends on the
#' parameters, the log-likelihood is discontinuous: whenever a parameter move
#' carries an observation's ceiling across an integer its normalizing constant
#' jumps, by about \eqn{(1-\alpha)^{k}} at ceiling \eqn{k}. The maximizer is
#' found by a deterministic multistart, BFGS from the Poisson solution at each
#' value of `alpha.start` polished by Nelder-Mead with feasibility-repaired
#' restarts, on covariates scaled to unit standard deviation (the coefficients
#' are mapped back, so the fit does not depend on the covariates' units). Two
#' parameterizations of the same design (say, treatment and sum contrasts) can
#' settle on different teeth of this surface, with log-likelihoods differing by
#' the order of the jumps, which is far inside the sampling variability of the
#' estimates; the fixed-effects estimator [cpb_fe()] solves the corresponding
#' problem in the unit intercepts exactly. The numerical Hessian is unreliable on
#' such a surface, so inference uses a cold-multistart bootstrap for the
#' coefficients (validated to nominal coverage) and a profile-likelihood interval
#' for \eqn{\alpha} (see [confint.cpb()]).
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param truncated Logical; if `TRUE` (default) fit the zero-truncated CPB
#'   (requires all `Y >= 1`); if `FALSE` fit the untruncated CPB on `Y >= 0`.
#'   **Note the sibling default differs:** [cpb_fe()] defaults to
#'   `truncated = FALSE`, because its typical call fits a whole panel including
#'   zeros, while `cpb()`'s typical call fits the positive counts of a
#'   zero-inflated outcome. State `truncated` explicitly when moving a
#'   specification between the two.
#' @param se Inference method: `"none"` (default; fast, no standard errors) or
#'   `"bootstrap"` (cold-multistart pairs/cluster bootstrap; slower).
#' @param B Number of bootstrap resamples (default 500).
#' @param cluster Optional cluster identifier for cluster-robust inference: a
#'   column name in `data` or a vector aligned to its rows. When supplied, the
#'   bootstrap resamples whole clusters (a block bootstrap), giving
#'   cluster-robust standard errors and intervals; when `NULL` (default) it
#'   resamples observations, giving heteroskedasticity-robust errors. This is
#'   the appropriate route to robust inference here because the CPB's
#'   parameter-dependent support makes a Hessian-based sandwich unreliable.
#' @param offset Optional offset on the log-mean scale (an exposure): a numeric
#'   vector or the name of a column in `data`, making the model a rate model.
#' @param weights Optional frequency weights: a numeric vector or the name of a
#'   column in `data`. A weight of \eqn{w_i} is equivalent to \eqn{w_i} copies of
#'   row \eqn{i}; the bootstrap resamples rows together with their weights, and
#'   `nobs()` returns the weight total.
#' @param cores Number of worker processes for the bootstrap (default 1). With
#'   `cores > 1` the replicates run on a socket cluster whose random streams are
#'   seeded from the calling session, so `set.seed()` before the call
#'   reproduces a parallel run (parallel and serial draws differ).
#' @param alpha.start Starting value for the dispersion parameter (default 0.5); the
#'   fit also multi-starts over a spread of alpha values.
#' @param max.support Guard on the maximum evaluated support: parameter values
#'   whose implied ceiling exceeds it are treated as infeasible (they lie in the
#'   `alpha` near 1 region where the support explodes). The default `NULL`
#'   sets it to `max(500, 10 * max(y))`. A fit whose largest ceiling reaches the
#'   guard is flagged with a warning and in `$support_binding`, because `alpha`
#'   is then bounded by the guard rather than by the data.
#' @param maxit,reltol Optimizer controls passed to [stats::optim()].
#'
#' @return An object of class `"cpb"`: a list with `coefficients`, `alpha`, bootstrap
#'   `se.beta`/`vcov`/`ci.beta`, `loglik`, `loglik.null`, `fitted.values`, `ceiling`
#'   (the observation-specific implied ceiling), the model frame pieces, and the
#'   original `call`.
#'
#' @references King, G. (1989). Variance specification in event count models.
#'   \emph{American Journal of Political Science}, 33(3), 762-784.
#'
#' @examples
#' set.seed(1)
#' n <- 300; x <- rnorm(n)
#' N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
#' y <- rbinom(n, N, 0.5)                 # underdispersed (var/mean approx 0.5)
#' d <- data.frame(y = y, x = x)
#' fit <- cpb(y ~ x, data = d[d$y > 0, ], se = "none")
#' summary(fit)
#' @seealso [ud_screen()], [confint.cpb()], [predict.cpb()], [implied_ceiling()],
#'   [dispersion_test()]
#' @export
cpb <- function(formula, data, truncated = TRUE, se = c("none", "bootstrap"),
                B = 500, cluster = NULL, offset = NULL, weights = NULL, cores = 1L,
                alpha.start = 0.5, max.support = NULL, maxit = 20000, reltol = 1e-8) {
  se <- match.arg(se)
  cl <- match.call()
  mf <- model.frame(formula, data, na.action = na.omit)
  Y  <- model.response(mf); X <- model.matrix(formula, mf)
  n  <- length(Y); p <- ncol(X); rows <- .ud_kept_rows(mf, data)
  if (!is.numeric(Y)) stop("Response must be a numeric count; got ", class(Y)[1L], ".")
  if (n == 0L) stop("No complete cases on the model variables.")
  off <- rep_len(0, n)                                          # log-scale exposure offset
  if (!is.null(offset)) {
    ov <- if (is.character(offset) && length(offset) == 1L) data[[offset]] else offset
    if (is.null(ov)) stop("'offset' must name a column of 'data' or be a numeric vector.")
    if (length(ov) != nrow(data)) stop("'offset' must have one value per row of 'data'.")
    off <- as.numeric(ov[rows])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  }
  w <- .ud_weights(weights, data, rows)

  ## optional clustering variable, aligned to the (possibly na.omit-reduced) model frame
  clid <- NULL
  if (!is.null(cluster)) {
    if (se != "bootstrap") {
      warning("'cluster' only affects the bootstrap; it is ignored with se = \"none\" (use se = \"bootstrap\").")
    } else {
      cv <- if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
      if (is.null(cv)) stop("'cluster' must be a column name in 'data' or a vector.")
      if (length(cv) != nrow(data)) stop("'cluster' must have one value per row of 'data'.")
      clid <- cv[rows]
      if (anyNA(clid)) stop("'cluster' has missing values on the estimation rows.")
    }
  }
  if (any(Y < 0) || any(Y != floor(Y)))
    stop("Response must be non-negative integer counts.")
  if (truncated && any(Y < 1))
    stop("truncated = TRUE requires all Y >= 1; use truncated = FALSE for zero-inclusive data.")
  .ud_rank_check(X)
  if (is.null(max.support)) max.support <- max(500L, 10L * max(Y))

  astarts <- unique(c(alpha.start, 0.15, 0.4, 0.65, 0.9))
  fit <- .cpb_fit_ms(X, Y, off, max.support, truncated, astarts, maxit, reltol, w = w)
  if (is.null(fit)) stop("All optimization starts failed.")
  if (fit$value >= 1e9)
    stop("no feasible fit: every start left some count above its implied ceiling (raise max.support, ",
         "currently ", max.support, ", or check the data).")
  if (fit$convergence != 0)
    warning("Optimization did not converge (code ", fit$convergence, ").")

  beta  <- fit$par[1:p]; names(beta) <- colnames(X)
  alpha <- plogis(fit$par[p + 1]); loglik <- -fit$value
  lam   <- as.numeric(exp(off + X %*% beta))
  fitted <- .cpb_mean(lam, alpha, truncated)                # exact mean of the fitted pmf
  binding <- max(lam / (1 - alpha)) >= 0.95 * max.support
  if (binding)
    warning("the fitted ceiling reaches max.support (", max.support, "): alpha is bounded by the guard, not the data. ",
            "The data may not be underdispersed (compare gec(), which lets the dispersion go either way), ",
            "or raise max.support.")
  loglik0 <- .cpb_baseline_loglik(X, Y, off, truncated, w = w)

  se.beta <- setNames(rep(NA_real_, p), colnames(X)); vcv <- NULL
  boot <- NULL; ci.beta <- NULL
  if (se == "bootstrap") {
    if (!is.null(clid)) .ud_cluster_guard(clid)
    if (!is.null(clid)) .ud_cluster_guard(clid)
    grp <- if (!is.null(clid)) split(seq_len(n), clid) else NULL   # NULL = ordinary bootstrap
    one <- function(b) {
      idx <- if (is.null(grp)) sample.int(n, n, replace = TRUE)     # resample rows, or
             else unlist(grp[sample.int(length(grp), length(grp), replace = TRUE)],
                         use.names = FALSE)                          # resample whole clusters
      fb <- .cpb_fit_ms(X[idx, , drop = FALSE], Y[idx], off[idx], max.support, truncated,
                        c(0.2, 0.5, 0.8), maxit, 1e-6, w = if (is.null(w)) NULL else w[idx], restarts = 0L)
      if (!is.null(fb) && fb$convergence == 0 && fb$value < 1e9) {
        bp <- fb$par; bp[p + 1] <- plogis(bp[p + 1]); bp
      } else rep(NA_real_, p + 1L)
    }
    boot <- do.call(rbind, .ud_lapply(seq_len(B), one, cores))
    ok <- boot[complete.cases(boot), , drop = FALSE]
    if (nrow(ok) < 2) warning("Too few bootstrap resamples converged for stable inference.")
    se.beta <- setNames(apply(ok[, 1:p, drop = FALSE], 2, sd), colnames(X))
    vcv <- cov(ok[, 1:p, drop = FALSE]); dimnames(vcv) <- list(colnames(X), colnames(X))
    ci.beta <- t(apply(ok[, 1:p, drop = FALSE], 2, quantile, c(.025, .975)))
    rownames(ci.beta) <- colnames(X)
    attr(boot, "nboot_ok") <- nrow(ok)
  }

  structure(list(
    coefficients = beta, alpha = alpha, se.beta = se.beta, vcov = vcv,
    ci.beta = ci.beta, boot = boot, loglik = loglik, loglik.null = loglik0,
    fitted.values = fitted, rate = lam, linear.predictors = as.numeric(off + X %*% beta), offset = off,
    weights = w, nobs_weighted = if (is.null(w)) n else sum(w),
    ceiling = lam / (1 - alpha), residuals = Y - fitted,
    n = n, df = p + 1, p = p, truncated = truncated, max.support = max.support, support_binding = binding,
    se.type = se, clustered = !is.null(clid),
    n_clusters = if (!is.null(clid)) length(unique(clid)) else NA_integer_,
    converged = (fit$convergence == 0),
    formula = formula, terms = attr(mf, "terms"), levels = .cpb_xlevels(mf),
    contrasts = attr(X, "contrasts"), X = X, Y = Y, call = cl),
    class = "cpb")
}

## Profile-likelihood interval for alpha: invert the LR test by fixing alpha and
## re-optimizing beta (on the column-scaled design). Well calibrated except at
## strong underdispersion, where alpha sits at the feasibility boundary and the
## interval becomes one-sided (see paper).
## The profile is traced by continuation: each alpha starts from the solution
## at the previous alpha, and the start is first repaired to feasibility (a
## smaller alpha shrinks every ceiling lambda_i/(1-alpha), so the previous
## coefficients can leave observations above their ceiling; raising the
## intercept restores feasibility) and then re-maximized over beta with
## restarts. The lower limit is bracketed on a descending grid and located by
## uniroot inside the bracket. The interval is one-sided (`boundary`) only when
## the profile has not dropped to the cut by alpha = 0.005.
.cpb_alpha_profile_ci <- function(object, level = 0.95) {
  X <- object$X; Y <- object$Y; ms <- object$max.support; tr <- object$truncated
  off <- if (!is.null(object$offset)) object$offset else rep_len(0, length(Y))
  w <- object$weights
  s <- .ud_colscale(X); Xs <- sweep(X, 2, s, "/")
  llmax <- object$loglik; cut <- llmax - qchisq(level, 1) / 2
  icol <- match("(Intercept)", colnames(X))
  pos <- Y > 0
  repair <- function(b, a) {                     # shift the intercept so every count is feasible
    if (is.na(icol) || !any(pos)) return(b)
    eta <- off + as.numeric(Xs %*% b)
    need <- max(log(Y[pos] * (1 - a)) - eta[pos]) + 1e-6
    if (need > 0) b[icol] <- b[icol] + need * s[icol]
    b
  }
  prof <- function(a, b0) {
    b0 <- repair(b0, a)
    fn <- function(bb) .cpb_nll(c(bb, qlogis(a)), Xs, Y, off, ms, tr, w)
    o <- optim(b0, fn, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-8))
    o <- .ud_nm_polish(o, fn, 3000, 1e-8, restarts = 2L)
    list(ll = -o$value, b = o$par)
  }
  ah <- object$alpha; bh <- object$coefficients * s
  ## downward: descending grid with continuation until the profile crosses the cut
  lo <- NA_real_; a_prev <- ah; b_prev <- bh; ll_prev <- llmax
  grid <- rev(seq(0.005, ah - 1e-4, by = max(0.005, ah / 40)))
  for (a in grid) {
    r <- prof(a, b_prev)
    if (r$ll < cut) {
      f <- function(aa) prof(aa, b_prev)$ll - cut
      lo <- tryCatch(uniroot(f, c(a, a_prev), tol = 1e-4)$root, error = function(e) a)
      break
    }
    a_prev <- a; b_prev <- r$b; ll_prev <- r$ll
  }
  ## upward: the previous solution stays feasible as alpha rises
  hi <- tryCatch(uniroot(function(aa) prof(aa, bh)$ll - cut, c(ah + 1e-4, 0.999), tol = 1e-4)$root,
                 error = function(e) NA_real_)
  boundary <- is.na(lo)
  c(lower = if (is.na(lo)) 0.005 else unname(lo),
    upper = if (is.na(hi)) 0.999 else unname(hi),
    boundary = as.numeric(boundary))
}
