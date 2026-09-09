// cpb_fe.cpp -- CPB regression with high-dimensional unit fixed effects via a
// CONCENTRATED (profiled) likelihood. Each unit's intercept is solved by a 1-D inner
// maximization, so the outer optimizer sees only (beta, alpha). No unit dummies.
// Data must be SORTED BY UNIT; ustart[u]..ustart[u+1]-1 are unit u's rows.
//
// The inner problem. Observation t has rate lambda_t = exp(a + o_t) and ceiling
// N_t = lambda_t / (1 - alpha), and its pmf is renormalized over 0..floor(N_t).
// As a increases through a breakpoint a = log(k (1-alpha)) - o_t (N_t = k) the
// support gains the point k, the normalizer grows, and the log pmf of the
// observed y_t < k drops by about (1-alpha)^k. The unit log-likelihood is
// therefore a saw-tooth in a: smooth between breakpoints, with a downward jump
// at every breakpoint, and its smooth part is unimodal. Its supremum over a is
// attained either at the interior critical point of the tooth that contains the
// peak of the smooth part or as the left limit at a breakpoint, never at a
// tooth's left end. The search below (i) locates the peak of the smooth part
// (a coarse grid around the unit's Poisson intercept; a warm start may be
// supplied, but the package evaluates every unit from cold, since a warm start
// carried from a distant trial point can leave a unit on the wrong tooth),
// (ii) evaluates the
// left limit at every breakpoint within a window of that peak, and (iii)
// refines by golden section inside the tooth holding the peak and inside the
// tooth ending at the best breakpoint. Solving the inner problem to its
// supremum is what makes the concentrated likelihood smooth in alpha (a change
// in alpha shifts every breakpoint of a unit rigidly, so the maximizing tooth
// slides) and continuous with kinks in beta.
//
// The gradient. By the envelope theorem the gradient of the concentrated
// log-likelihood is the partial gradient at the solved intercepts; where a
// unit's supremum is the left limit at a breakpoint of observation t*, the
// intercept tracks that breakpoint, a* = log(k (1-alpha)) - o_{t*}, so the
// total derivative adds (dl/da) (da*/dtheta) with da*/dbeta = -x_{t*} and
// da*/dalpha = -1/(1-alpha).
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// log P_CPB(Y=y | lambda, alpha), using a precomputed log-factorial table lf.
static inline double cpb_logpmf(int y, double lam, double alpha, double la, double l1a,
                                const std::vector<double>& lf, int max_support, bool truncated) {
  double ni = lam / (1.0 - alpha);
  int Ki = (int)std::floor(ni + 1e-9);
  if (y > Ki || Ki > max_support) return -1e300;
  double lgni1 = R::lgammafn(ni + 1.0);
  double g = lgni1, mx = -1e300;
  // first pass: max (the work buffer is reused across calls: no allocation). The
  // log terms are unimodal in k, so the sum stops once it is past the mode and
  // 36 log units below it (and past y), which bounds the work at large ceilings.
  static std::vector<double> lw; if ((int)lw.size() < Ki + 1) lw.resize(Ki + 1);
  int kend = Ki;
  for (int k = 0; k <= Ki; k++) {
    if (k > 0) g -= std::log(ni - k + 1.0);
    lw[k] = lgni1 - lf[k] - g + k * l1a + (ni - k) * la;
    if (lw[k] > mx) mx = lw[k];
    else if (k >= y && lw[k] < mx - 36.0) { kend = k; break; }
  }
  double s = 0.0; for (int k = 0; k <= kend; k++) s += std::exp(lw[k] - mx);
  double logD = mx + std::log(s);
  double lp = lw[y] - logD;
  if (truncated) {
    double e = std::exp(lw[0] - logD);
    if (e >= 1.0 - 1e-15) return -1e300;
    lp -= std::log1p(-e);
  }
  return lp;
}

// d log P_CPB(y | lambda, alpha) / d eta (eta = log lambda) and / d alpha.
// Writing psi_k = digamma(N - k + 1) and E[.] for the expectation under the
// (truncated-support) pmf, the terms common to every k cancel and
//   d/deta   = N (E[psi_k] - psi_y),
//   d/dalpha = -(y - E[k]) (1/(1-alpha) + 1/alpha) - N/(1-alpha) (psi_y - E[psi_k]);
// with zero truncation the same quantities at y = 0 enter through
//   -d log(1 - p0) = p0/(1-p0) d log p0.
static inline bool cpb_dlogpmf(int y, double lam, double alpha, double la, double l1a,
                               const std::vector<double>& lf, int max_support, bool truncated,
                               double& deta, double& dalpha) {
  double ni = lam / (1.0 - alpha);
  int Ki = (int)std::floor(ni + 1e-9);
  if (y > Ki || Ki > max_support) return false;
  double lgni1 = R::lgammafn(ni + 1.0);
  double g = lgni1, mx = -1e300;
  static std::vector<double> lw, ps; if ((int)lw.size() < Ki + 1) { lw.resize(Ki + 1); ps.resize(Ki + 1); }
  int kend = Ki;
  for (int k = 0; k <= Ki; k++) {
    if (k > 0) g -= std::log(ni - k + 1.0);
    lw[k] = lgni1 - lf[k] - g + k * l1a + (ni - k) * la;
    ps[k] = R::digamma(ni - k + 1.0);
    if (lw[k] > mx) mx = lw[k];
    else if (k >= y && lw[k] < mx - 36.0) { kend = k; break; }
  }
  double s = 0.0, epsi = 0.0, ek = 0.0;
  for (int k = 0; k <= kend; k++) { double w = std::exp(lw[k] - mx); s += w; epsi += w * ps[k]; ek += w * k; }
  epsi /= s; ek /= s;
  double c = 1.0 / (1.0 - alpha) + 1.0 / alpha, r = ni / (1.0 - alpha);
  deta   = ni * (epsi - ps[y]);
  dalpha = -((double)y - ek) * c - r * (ps[y] - epsi);
  if (truncated) {
    double p0 = std::exp(lw[0] - mx) / s;
    if (p0 >= 1.0 - 1e-15) return false;
    double f = p0 / (1.0 - p0);
    deta   += f * (ni * (epsi - ps[0]));
    dalpha += f * (-(0.0 - ek) * c - r * (ps[0] - epsi));
  }
  return true;
}

// Everything the inner search needs about one unit.
struct Unit {
  const IntegerVector& Y; const std::vector<double>& o; const NumericVector& w;
  int lo, hi; double alpha, la, l1a; const std::vector<double>& lf; int max_support; bool truncated;
  double ll(double a) const {
    double v = 0.0;
    for (int t = lo; t < hi; t++) {
      double lp = cpb_logpmf(Y[t], std::exp(a + o[t]), alpha, la, l1a, lf, max_support, truncated);
      if (lp <= -1e299 || !R_finite(lp)) return -1e300;
      v += w[t] * lp;
    }
    return v;
  }
  // breakpoints bounding the tooth that contains a (prev <= a < next) and the
  // observations that own them
  void tooth(double a, double& prev, double& next, int& prev_t, int& next_t) const {
    prev = -1e300; next = 1e300; prev_t = -1; next_t = -1;
    for (int t = lo; t < hi; t++) {
      double N = std::exp(a + o[t]) / (1.0 - alpha);
      int K = (int)std::floor(N + 1e-9);
      if (K >= 1) { double b = std::log((double)K * (1.0 - alpha)) - o[t]; if (b > prev) { prev = b; prev_t = t; } }
      double b2 = std::log((double)(K + 1) * (1.0 - alpha)) - o[t]; if (b2 < next) { next = b2; next_t = t; }
    }
  }
  void tooth(double a, double& prev, double& next) const { int pt, nt; tooth(a, prev, next, pt, nt); }
};

// window (log scale) around the peak within which breakpoints are enumerated
// (extended once on a side where the best breakpoint sits near the edge), and
// the offset below a breakpoint at which its left limit is taken
static const double BP_WINDOW = 0.15;
static const double BP_EPS    = 1e-8;
static const double GR = 0.6180339887498949;

// golden section for a maximum on [L, H]: at most `iters` iterations, stopping
// once the bracket is narrower than 1e-6 on the log scale (an intercept error
// of that size moves the objective by less than 1e-9); returns the value, a_at
// the location
static double golden_max(const Unit& U, double L, double H, int iters, double& a_at) {
  double c = H - GR * (H - L), d = L + GR * (H - L), fc = U.ll(c), fd = U.ll(d);
  for (int it = 0; it < iters && (H - L) > 1e-6; it++) {
    if (fc < fd) { L = c; c = d; fc = fd; d = L + GR * (H - L); fd = U.ll(d); }
    else         { H = d; d = c; fd = fc; c = H - GR * (H - L); fc = U.ll(c); }
  }
  a_at = 0.5 * (L + H);
  return U.ll(a_at);
}

// best left limit over the breakpoints inside (WL, WH]; best_t = its observation
static void scan_breakpoints(const Unit& U, double WL, double WH, double L0,
                             double& best_a, double& best_v, bool& best_bp, int& best_t) {
  // a breakpoint at ceiling k drops the objective by about (1-alpha)^k: only
  // ceilings whose jump exceeds 1e-7 are worth visiting
  double kcap = std::log(1e-7) / std::log1p(-U.alpha);
  int kmax_teeth = (kcap > (double)U.max_support) ? U.max_support : (int)std::floor(kcap);
  for (int t = U.lo; t < U.hi; t++) {
    double nlo = std::exp(WL + U.o[t]) / (1.0 - U.alpha), nhi = std::exp(WH + U.o[t]) / (1.0 - U.alpha);
    int klo = std::max((int)std::ceil(nlo), U.Y[t] + 1), khi = (int)std::floor(nhi);
    if (khi > kmax_teeth) khi = kmax_teeth;
    for (int k = klo; k <= khi; k++) {
      double b = std::log((double)k * (1.0 - U.alpha)) - U.o[t] - BP_EPS;
      if (b <= L0 || b <= WL) continue;
      double v = U.ll(b);
      if (v > best_v) { best_v = v; best_a = b; best_bp = true; best_t = t; }
    }
  }
}

// returns the unit's maximized log-likelihood; a_star = the intercept; edge_t =
// the row whose breakpoint the intercept sits at (-1 for an interior optimum)
static double maximize_unit(const Unit& U, double a0, double& a_star, int& edge_t, int inner_it) {
  double alo = -1e300; int ymax = 0; double sy = 0.0, se = 0.0;
  edge_t = -1;
  for (int t = U.lo; t < U.hi; t++) {
    sy += U.w[t] * (double)U.Y[t]; se += U.w[t] * std::exp(U.o[t]);
    if (U.Y[t] > 0) {
      double b = std::log(U.Y[t] * (1.0 - U.alpha)) - U.o[t];   // ceiling >= count
      if (b > alo) alo = b;
      if (U.Y[t] > ymax) ymax = U.Y[t];
    }
  }
  if (ymax == 0) { a_star = -30.0; return 0.0; }
  double apois = (sy > 0.0 && se > 0.0) ? std::log(sy / se) : alo;   // Poisson unit intercept
  double L0 = alo + 1e-6;
  bool warm = R_finite(a0), best_bp = false;
  double best_a = L0, best_v = -1e300, centre = L0; int best_t = -1;
  for (int pass = 0; pass < 2; pass++) {
    // (i) peak of the smooth part: warm start, or coarse grid + golden refinement
    if (warm) centre = std::max(a0, L0);
    else {
      double L = std::max(L0, apois - 2.0), H = std::max(L + 0.5, apois + 2.0);
      double bv = -1e300; centre = L;
      for (int ext = 0; ext < 3; ext++) {
        int G = (int)std::ceil((H - L) / 0.5);
        for (int g = 0; g <= G; g++) {
          double a = L + (H - L) * (double)g / (double)G;
          double v = U.ll(a);
          if (v > bv) { bv = v; centre = a; }
        }
        if (centre >= H - 1e-9) { L = H; H = H + 2.0; continue; }          // peak beyond the grid: extend
        if (centre <= L + 1e-9 && L > L0 + 1e-9) { H = L; L = std::max(L0, L - 2.0); continue; }
        break;
      }
      double ac; double vc = golden_max(U, std::max(L0, centre - 0.5), centre + 0.5, 12, ac);
      if (vc > bv) centre = ac;
    }
    // (ii) left limits at every breakpoint within the window. The window is the
    // base width plus the widest tooth of any observation at the centre
    // (log((K+1)/K) at ceiling K), so that units with small ceilings, whose
    // teeth are wide and tall, are searched across whole teeth.
    double wmax = 0.0;
    for (int t = U.lo; t < U.hi; t++) {
      int K = (int)std::floor(std::exp(centre + U.o[t]) / (1.0 - U.alpha) + 1e-9);
      if (K < 1) K = 1;
      double wt = std::log1p(1.0 / (double)K); if (wt > wmax) wmax = wt;
    }
    double W = std::min(1.0, BP_WINDOW + wmax);
    double WL = std::max(L0, centre - W), WH = centre + W;
    best_a = centre; best_v = U.ll(centre); best_bp = false; best_t = -1;
    scan_breakpoints(U, WL - 1e-12, WH, L0, best_a, best_v, best_bp, best_t);
    if (best_a >= WH - 0.02) { scan_breakpoints(U, WH, WH + W, L0, best_a, best_v, best_bp, best_t); WH += W; }
    else if (best_a <= WL + 0.02 && WL > L0 + 1e-9) {
      double WL2 = std::max(L0, WL - W);
      scan_breakpoints(U, WL2 - 1e-12, WL, L0, best_a, best_v, best_bp, best_t); WL = WL2;
    }
    bool at_edge = (best_a >= WH - 0.02) || (best_a <= WL + 0.02 && WL > L0 + 1e-9);
    if (warm && at_edge) { warm = false; continue; }         // peak left the warm window: cold pass
    break;
  }
  if (best_v <= -1e299) { a_star = L0; return U.ll(L0); }
  // (iii) interior critical points: golden section inside the tooth holding the
  // peak and, if the best is a breakpoint, inside the tooth ending there; never
  // below the incumbent
  double prev, next, ag, vg, prev2, next2;
  U.tooth(centre, prev, next);
  double L = std::max(L0, prev + BP_EPS), H = next - BP_EPS;
  if (H > L) { vg = golden_max(U, L, H, inner_it, ag); if (vg > best_v) { best_v = vg; best_a = ag; best_bp = false; } }
  if (best_bp) {
    U.tooth(best_a - BP_EPS, prev2, next2);
    if (std::fabs(prev2 - prev) > 1e-12) {                  // a different tooth from the one just searched
      L = std::max(L0, prev2 + BP_EPS); H = best_a;
      if (H > L) { vg = golden_max(U, L, H, inner_it, ag); if (vg > best_v) { best_v = vg; best_a = ag; best_bp = false; } }
    }
  }
  // A golden-section optimum at an end of its tooth is an edge too: at the
  // right end it is the left limit of the next breakpoint; at the left end
  // (including the feasibility floor, where an observation's ceiling equals
  // its count) it sits on the breakpoint below and tracks it as beta moves.
  // Both need the breakpoint correction in the gradient.
  if (!best_bp) {
    int pt, nt; U.tooth(best_a, prev, next, pt, nt);
    if (best_a >= next - 1e-5 && nt >= 0) {
      double aH = next - BP_EPS, vH = U.ll(aH);
      if (vH >= best_v - 1e-12) { best_v = vH; best_a = aH; }
      best_bp = true; best_t = nt;
    } else if (best_a <= prev + 1e-5 && pt >= 0) {
      double aL = std::max(L0, prev + BP_EPS), vL = U.ll(aL);
      if (vL >= best_v - 1e-12) { best_v = vL; best_a = aL; }
      best_bp = true; best_t = pt;
    }
  }
  a_star = best_a; edge_t = best_bp ? best_t : -1;
  return best_v;
}

static void compute_offsets(NumericVector params, NumericMatrix X, NumericVector uoff, std::vector<double>& o) {
  int n = X.nrow(), p = X.ncol();
  for (int i = 0; i < n; i++) { double e = uoff[i]; for (int j = 0; j < p; j++) e += X(i, j) * params[j]; o[i] = e; }
}

// Concentrated negative log-likelihood at (beta, logit alpha). `awarm` carries the
// unit intercepts from the previous evaluation (length n_units; length 0 = cold
// start for every unit). Returns list(nll, a, edge): the intercepts, and for
// each unit the row whose breakpoint the intercept sits at (-1 if interior),
// which the gradient needs; on an infeasible point nll = 1e10 and `a` is
// `awarm` unchanged.
// [[Rcpp::export]]
List cpb_fe_nll_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                    NumericVector w, IntegerVector ustart, int n_units, int max_support,
                    bool truncated, int inner_it, NumericVector awarm) {
  int n = X.nrow(), p = X.ncol();
  if (ustart.size() != n_units + 1 || ustart[n_units] != n || Y.size() != n || w.size() != n || offset.size() != n)
    Rcpp::stop("internal error: unit index does not match the design (rows dropped after indexing)");
  bool warm = awarm.size() == n_units;
  IntegerVector edge(n_units, -1);
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  if (alpha <= 1e-12 || alpha >= 1.0 - 1e-12) return List::create(1e10, awarm, edge);
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  std::vector<double> lf(max_support + 2); lf[0] = 0.0;
  for (int k = 1; k < (int)lf.size(); k++) lf[k] = lf[k-1] + std::log((double)k);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector a(n_units);
  double total = 0.0, astar; int et;
  for (int u = 0; u < n_units; u++) {
    Unit U = { Y, o, w, ustart[u], ustart[u+1], alpha, la, l1a, lf, max_support, truncated };
    double a0 = warm ? awarm[u] : NA_REAL;
    double ll = maximize_unit(U, a0, astar, et, inner_it);
    if (ll <= -1e299) return List::create(1e10, awarm, edge);
    a[u] = astar; edge[u] = et; total += ll;
  }
  return List::create(-total, a, edge);
}

// Gradient of the concentrated negative log-likelihood with respect to
// (beta, logit alpha) at the solved intercepts `a` (see the header comment).
// [[Rcpp::export]]
NumericVector cpb_fe_grad_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                              NumericVector w, IntegerVector ustart, int n_units, int max_support,
                              bool truncated, NumericVector a, IntegerVector edge) {
  int n = X.nrow(), p = X.ncol();
  if (ustart.size() != n_units + 1 || ustart[n_units] != n || a.size() != n_units || edge.size() != n_units)
    Rcpp::stop("internal error: gradient inputs do not match the design");
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  std::vector<double> lf(max_support + 2); lf[0] = 0.0;
  for (int k = 1; k < (int)lf.size(); k++) lf[k] = lf[k-1] + std::log((double)k);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector G(p + 1);
  for (int u = 0; u < n_units; u++) {
    double Ga = 0.0, Galpha = 0.0; std::vector<double> Gb(p, 0.0);
    for (int t = ustart[u]; t < ustart[u+1]; t++) {
      double de, da;
      if (!cpb_dlogpmf(Y[t], std::exp(a[u] + o[t]), alpha, la, l1a, lf, max_support, truncated, de, da))
        return NumericVector(p + 1, NA_REAL);
      Ga += w[t] * de; Galpha += w[t] * da;
      for (int j = 0; j < p; j++) Gb[j] += w[t] * de * X(t, j);
    }
    if (edge[u] >= 0) {                                   // the intercept tracks a breakpoint
      int ts = edge[u];
      for (int j = 0; j < p; j++) Gb[j] -= Ga * X(ts, j);
      Galpha -= Ga / (1.0 - alpha);
    }
    for (int j = 0; j < p; j++) G[j] -= Gb[j];
    G[p] -= Galpha * alpha * (1.0 - alpha);                // chain rule to logit alpha
  }
  return G;
}

// Per-observation derivatives of the CPB log pmf (columns: d/d eta, d/d alpha),
// exported for verification against numerical differences of dcpb().
// [[Rcpp::export]]
NumericMatrix cpb_dlogpmf_cpp(IntegerVector y, NumericVector lam, double alpha, int max_support, bool truncated) {
  int n = y.size();
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  std::vector<double> lf(max_support + 2); lf[0] = 0.0;
  for (int k = 1; k < (int)lf.size(); k++) lf[k] = lf[k-1] + std::log((double)k);
  NumericMatrix out(n, 2);
  for (int i = 0; i < n; i++) {
    double de, da;
    if (cpb_dlogpmf(y[i], lam[i], alpha, la, l1a, lf, max_support, truncated, de, da)) { out(i, 0) = de; out(i, 1) = da; }
    else { out(i, 0) = NA_REAL; out(i, 1) = NA_REAL; }
  }
  return out;
}

// Unit intercepts at (beta, logit alpha) from a cold search for every unit.
// [[Rcpp::export]]
NumericVector cpb_fe_intercepts_cpp(NumericVector params, NumericMatrix X, IntegerVector Y, NumericVector offset,
                                    NumericVector w, IntegerVector ustart, int n_units, int max_support,
                                    bool truncated, int inner_it) {
  int n = X.nrow(), p = X.ncol();
  if (ustart.size() != n_units + 1 || ustart[n_units] != n || Y.size() != n || w.size() != n || offset.size() != n)
    Rcpp::stop("internal error: unit index does not match the design (rows dropped after indexing)");
  double alpha = 1.0 / (1.0 + std::exp(-params[p]));
  double la = std::log(alpha), l1a = std::log1p(-alpha);
  std::vector<double> lf(max_support + 2); lf[0] = 0.0;
  for (int k = 1; k < (int)lf.size(); k++) lf[k] = lf[k-1] + std::log((double)k);
  std::vector<double> o(n); compute_offsets(params, X, offset, o);
  NumericVector fe(n_units); double astar; int et;
  for (int u = 0; u < n_units; u++) {
    Unit U = { Y, o, w, ustart[u], ustart[u+1], alpha, la, l1a, lf, max_support, truncated };
    maximize_unit(U, NA_REAL, astar, et, inner_it);
    fe[u] = astar;
  }
  return fe;
}
