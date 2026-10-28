## ---------------------------------------------------------------------------
## The participation equation of the hurdles (hurdle_cpb, hurdle_gec,
## hurdle_count): a binary glm of P(Y > 0) on every row.
##
## Its covariance is the glm's model-based one, or the sandwich with the meat
## of the count families' (.count_vcov): per-row scores (hurdle_count with
## se = "robust"), or scores summed within clusters with the G/(G - 1) factor
## (a fit given `cluster`), frequency weights entering as duplicated rows. The
## hurdles store it as `part_vcov`, with `part_vcov_type`, and every table,
## interval and quantity of interest reads it from there. A clustered covariance
## asked for without a cluster is NA, as count_reg()'s is (hurdle_count's
## intensity fit says so).
##
## Its separation is checked by the zero-inflated fits' rule (.zi_separation;
## see ?zi_cpb, Separation), which on a binary equation is the profile of the
## glm's own likelihood: a covariate coefficient whose profile 100 standard
## deviations of its covariate further from zero is less than
## qchisq(0.95, 1) / 2 below its value at the estimate is flagged. The intercept
## and the part_fe dummies are not checked (a unit that always or never
## participates has its dummy at infinity, as in any fixed-effects logit).
## ---------------------------------------------------------------------------

.part_fit_vcov <- function(pfit, type = "model", cluster = NULL) {
  V <- tryCatch(stats::vcov(pfit), error = function(e) NULL)
  if (is.null(V) || type == "model") return(V)
  if (type == "cluster" && is.null(cluster)) { V[] <- NA_real_; return(V) }
  ok <- !is.na(stats::coef(pfit)); X <- stats::model.matrix(pfit)[, ok, drop = FALSE]
  mu <- pfit$fitted.values; w <- pfit$prior.weights
  S <- X * ((pfit$y - mu) * pfit$family$mu.eta(pfit$linear.predictors) / pmax(mu * (1 - mu), 1e-300))
  meat <- if (type == "cluster") {
    G <- length(unique(cluster))
    (G / (G - 1)) * crossprod(as.matrix(rowsum(S * w, group = cluster)))
  } else crossprod(S * sqrt(w))
  bread <- V[ok, ok, drop = FALSE]
  V[ok, ok] <- bread %*% meat %*% bread
  V
}

## the participation equation's covariance of a hurdle (an object from before 0.2.0: the glm's)
.part_vcov <- function(object) {
  if (!is.null(object$part_vcov)) return(object$part_vcov)
  tryCatch(stats::vcov(object$participation), error = function(e) NULL)
}

## the separation check of the participation glm: its successes enter .zi_separation() as the "zeros" of a
## mixture whose count component has no zero mass, which leaves exactly the binary likelihood
.part_separation <- function(pfit, part_fe = NULL, group = NULL) {
  g <- stats::coef(pfit); ok <- !is.na(g)
  Z <- stats::model.matrix(pfit)[, ok, drop = FALSE]
  sep <- tryCatch(.zi_separation(Z, g[ok], ifelse(pfit$y > 0.5, 0L, 1L), rep(0, length(pfit$y)),
                                 w = pfit$prior.weights, group = group, fe_cols = .zi_fe_cols(Z, part_fe),
                                 link = pfit$family$link),
                  error = function(e) NULL)
  if (is.null(sep)) return(NULL)
  sep[c("diverging", "profile", "n_extreme", "groups_extreme")]
}

.part_warn_separation <- function(sep, fun) {
  if (is.null(sep) || !length(sep$diverging)) return(invisible(NULL))
  pl <- length(sep$diverging) > 1L
  warning(fun, ": the participation equation is separated or nearly so. The profile likelihood of the coefficient",
          if (pl) "s", " on ", paste(sep$diverging, collapse = ", "), " falls by less than qchisq(0.95, 1) / 2 over ",
          "100 standard deviations of ", if (pl) "their covariates" else "its covariate", " further from zero, so ",
          "the data put no useful bound on ", if (pl) "them" else "it", ", and ", if (pl) "their" else "its",
          " standard error, odds ratio and effect on the participation probability are not meaningful (",
          sep$n_extreme, " observations have a participation probability numerically 0 or 1).", call. = FALSE)
}

## The participation equation of a fitted hurdle: its covariance (model-based, robust, or clustered on `cl`) and
## separation check, as fields of the fit
.part_fields <- function(pfit, part_fe, data, type, cl, fun) {
  sep <- .part_separation(pfit, part_fe, if (!is.null(part_fe)) data[[part_fe]])
  .part_warn_separation(sep, fun)
  list(part_vcov = .part_fit_vcov(pfit, type, cl), part_vcov_type = type,
       part_n_clusters = if (type == "cluster") length(unique(cl)) else NA_integer_, part_separation = sep)
}

## the dummy columns of fixed-effects variables `fe` among the names `nm`
.is_fe_col <- function(nm, fe) {
  out <- rep(FALSE, length(nm))
  for (v in fe) out <- out | startsWith(nm, paste0("factor(", v, ")"))
  out
}

## The participation table of a hurdle's summary: the participation covariance's standard errors, the part_fe
## dummies left out, and the separated terms without a standard error
.part_block <- function(object) {
  b <- object$part_beta; V <- .part_vcov(object)
  se <- if (is.null(V)) NULL else sqrt(pmax(diag(V), 0))[names(b)]
  keep <- !.is_fe_col(names(b), object$part_fe) & !is.na(b)
  .sum_mask(.sum_block(b[keep], if (is.null(se)) NULL else se[keep]), object$part_separation$diverging)
}
.part_note <- function(object) {
  ty <- object$part_vcov_type
  c(if (identical(ty, "cluster") && isTRUE(object$part_n_clusters > 0))
      sprintf("Participation standard errors are cluster-robust (%d clusters).", object$part_n_clusters)
    else if (identical(ty, "robust")) "Participation standard errors are heteroskedasticity-robust.",
    .sep_note(object$part_separation, "participation"))
}

## a summary table's separated rows lose their standard error, z and p
.sum_mask <- function(bl, terms) {
  if (is.null(bl) || !length(terms) || ncol(bl) < 4L) return(bl)
  bl[rownames(bl) %in% terms, 2:4] <- NA
  bl
}
.sep_note <- function(sep, equation) {
  if (is.null(sep) || !length(sep$diverging)) return(NULL)
  paste0("Note: the data put no useful bound on ", paste(sep$diverging, collapse = ", "), " in the ", equation,
         " equation (separation); ", if (length(sep$diverging) > 1L) "their standard errors are" else
           "its standard error is", " not shown.")
}

## the separated terms of a fit (zero-inflation or participation equation), for irr(), zi_test() and compare_models()
.ud_separated <- function(object) c(object$separation$diverging, object$part_separation$diverging)
