## ---------------------------------------------------------------------------
## Internal utilities shared by the estimators: estimation-row bookkeeping,
## frequency weights, the parallel replicate driver used by every bootstrap.
## ---------------------------------------------------------------------------

## indices (into the data frame passed to model.frame) of the rows that
## survived na.omit -- independent of row names
.ud_kept_rows <- function(mf, data) {
  na <- attr(mf, "na.action")
  idx <- seq_len(nrow(data))
  if (is.null(na)) idx else idx[-as.integer(na)]
}

## resolve a `weights` argument (column name or numeric vector, aligned to
## `data`) to the estimation rows `rows`; NULL stays NULL
.ud_weights <- function(weights, data, rows) {
  if (is.null(weights)) return(NULL)
  w <- if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(data)) stop("'weights' must name a column of 'data' or be a numeric vector.")
    data[[weights]]
  } else weights
  if (!is.numeric(w)) stop("'weights' must be numeric (frequency weights).")
  if (length(w) != nrow(data)) stop("'weights' must have one value per row of 'data' (", nrow(data), ").")
  w <- as.numeric(w[rows])
  if (anyNA(w)) stop("'weights' has missing values on the estimation rows.")
  if (any(!is.finite(w)) || any(w < 0)) stop("'weights' must be finite and non-negative.")
  w
}

## weight vector for the likelihood: ones when no weights were supplied
.ud_w1 <- function(w, n) if (is.null(w)) rep(1, n) else w

## Parallel replicate driver for the bootstraps and the parametric-bootstrap
## screen threshold. `cores = 1` runs lapply(); `cores > 1` runs a PSOCK cluster
## (portable across Windows, macOS, and Linux) that loads this package from the
## calling session's library path and draws its random numbers from L'Ecuyer
## streams seeded from the calling session's RNG, so `set.seed()` before the
## call makes a parallel run reproducible (the parallel and serial streams
## differ, so a parallel result equals a serial one only in distribution).
.ud_lapply <- function(X, FUN, cores = 1L) {
  cores <- as.integer(cores)[1L]
  if (is.na(cores) || cores < 1L) cores <- 1L
  if (cores == 1L || length(X) < 2L) return(lapply(X, FUN))
  cores <- min(cores, length(X))
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(lp) { .libPaths(lp); loadNamespace("underdisp"); NULL }, .libPaths())
  parallel::clusterSetRNGStream(cl, iseed = sample.int(.Machine$integer.max, 1L))
  parallel::parLapply(cl, X, FUN)
}

## Column scales for the optimizer: each non-constant design column is divided
## by its standard deviation before optimization and the coefficients mapped
## back, so the fit (and the numerical Hessians differenced in the scaled
## coordinates) do not depend on the units the covariates are measured in.
## Constant columns (the intercept, a column of zeros) keep scale 1.
.ud_colscale <- function(X) {
  s <- apply(X, 2, function(v) {
    r <- stats::sd(v)
    if (!is.finite(r) || r <= 1e-10 * max(1, max(abs(v)))) 1 else r
  })
  s[colnames(X) %in% "(Intercept)"] <- 1
  unname(s)
}

## rank check for a design matrix: refuse an aliased column (an empty factor
## level, an exact linear combination) with a message naming it, rather than
## letting a flat likelihood direction reach the optimizer and the bootstrap
.ud_rank_check <- function(X, what = "design matrix") {
  if (ncol(X) == 0L) return(invisible(TRUE))
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    bad <- colnames(X)[qx$pivot[(qx$rank + 1L):ncol(X)]]
    stop("the ", what, " is rank deficient: column(s) ", paste(bad, collapse = ", "),
         " are collinear with the others (an empty factor level or an exact linear combination); ",
         "drop them or use droplevels().", call. = FALSE)
  }
  invisible(TRUE)
}

## within-unit rank check for the concentrated fixed-effects fits: a covariate
## constant within units is absorbed by the unit intercepts
.ud_rank_check_fe <- function(X, unit) {
  if (ncol(X) == 0L) return(invisible(TRUE))
  Xd <- X - apply(X, 2, function(v) stats::ave(v, unit))
  wvar <- apply(Xd, 2, function(v) max(abs(v)))
  const <- colnames(X)[wvar <= 1e-10 * pmax(1, apply(abs(X), 2, max))]
  if (length(const))
    stop("covariate(s) ", paste(const, collapse = ", "), " do not vary within units and are absorbed by ",
         "the unit fixed effects; remove them from the formula.", call. = FALSE)
  qx <- qr(Xd)
  if (qx$rank < ncol(X)) {
    bad <- colnames(X)[qx$pivot[(qx$rank + 1L):ncol(X)]]
    stop("the within-unit design is rank deficient: column(s) ", paste(bad, collapse = ", "),
         " are collinear with the others after removing the unit effects.", call. = FALSE)
  }
  invisible(TRUE)
}

## Nelder-Mead restart until the improvement stops paying (a single run can
## stall on a collapsed simplex short of the optimum)
.ud_nm_polish <- function(fit, fn, maxit, reltol, restarts = 3L, tol = 1e-6) {
  for (rs in seq_len(restarts)) {
    op2 <- tryCatch(stats::optim(fit$par, fn, method = "Nelder-Mead",
                                 control = list(maxit = maxit, reltol = reltol)), error = function(e) NULL)
    if (is.null(op2)) break
    if (op2$value < fit$value) fit <- op2
    if (is.null(op2) || fit$value > op2$value - tol) break
  }
  fit
}

## guard for a cluster bootstrap: no resampling variation with one cluster
.ud_cluster_guard <- function(clid) {
  ng <- length(unique(clid))
  if (ng < 2L) stop("the cluster bootstrap needs at least two clusters; 'cluster' is constant on the estimation rows.", call. = FALSE)
  if (ng < 30L) warning(sprintf("only %d clusters; the cluster bootstrap may be unreliable with so few.", ng), call. = FALSE)
  invisible(ng)
}

## first differences hold every other model-matrix column at its reference
## value; a variable that also enters an interaction cannot be moved through one
## column alone
.fd_check_terms <- function(variable, cols) {
  inter <- cols[grepl(":", cols, fixed = TRUE)]
  hit <- inter[vapply(strsplit(inter, ":", fixed = TRUE), function(p) variable %in% p, logical(1))]
  if (length(hit))
    stop("'", variable, "' enters the interaction term(s) ", paste(hit, collapse = ", "),
         "; first_difference() moves one model-matrix column and cannot recompute the product. ",
         "Use predict(newdata = ) with the two full covariate profiles instead.", call. = FALSE)
  invisible(TRUE)
}

## Outer optimizer for the concentrated fixed-effects likelihoods (cpb_fe, gec_fe).
## `nll_raw(par, awarm)` returns list(nll, unit intercepts, edge rows) and
## `grad_raw(par, a, edge)` the gradient of the concentrated objective at those
## intercepts (the envelope theorem with the breakpoint correction; see
## src/cpb_fe.cpp). The closure caches the last evaluation for the gradient.
## Each dispersion start runs BFGS from the Poisson slopes (Nelder-Mead if BFGS
## fails); the incumbent is polished by Nelder-Mead. The returned `$a` are the
## intercepts that attain `$value`.
.fe_outer <- function(nll_raw, grad_raw, bs, dstarts, maxit, reltol) {
  ## every evaluation searches each unit from cold: a warm start carried from a
  ## distant trial point can leave a unit on the wrong tooth, and the resulting
  ## path dependence defeats the gradient methods
  env <- new.env(); env$last <- NULL
  fn <- function(par) {
    r <- nll_raw(par, numeric(0)); v <- r[[1L]]
    if (!is.finite(v)) v <- 1e10                       # an overflowed trial point is infeasible
    if (v < 1e9) env$last <- list(par = as.numeric(par), a = r[[2L]], edge = r[[3L]])
    v
  }
  gr <- function(par) {
    if (is.null(env$last) || !identical(as.numeric(par), env$last$par)) fn(par)
    if (is.null(env$last) || !identical(as.numeric(par), env$last$par)) return(rep(0, length(par)))
    g <- grad_raw(par, env$last$a, env$last$edge)
    if (anyNA(g) || !all(is.finite(g))) rep(0, length(par)) else g
  }
  run <- function(start) {
    o <- tryCatch(stats::optim(start, fn, gr, method = "BFGS", control = list(maxit = maxit, reltol = reltol)),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$value) || o$value >= 1e9)
      o <- tryCatch(stats::optim(start, fn, method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol)),
                    error = function(e) NULL)
    o
  }
  cand <- lapply(dstarts, function(d0) run(c(bs, d0)))
  cand <- cand[!vapply(cand, is.null, logical(1))]
  if (!length(cand)) return(NULL)
  best <- cand[[which.min(vapply(cand, function(f) f$value, numeric(1)))]]
  if (best$value >= 1e9) return(best)
  best <- .ud_nm_polish(best, fn, maxit, reltol, restarts = 1L)
  rc <- nll_raw(best$par, numeric(0)); best$value <- rc[[1L]]; best$a <- rc[[2L]]
  best
}
