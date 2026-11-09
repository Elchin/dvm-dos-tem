/*! \file EnthalpySolver.cpp
 * Stage B standalone enthalpy solver (not used by the model). See EnthalpySolver.h.
 */
#include "../include/EnthalpySolver.h"
#include "../include/physicalconst.h"
#include "../include/SJCurve.h"

#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

namespace enthalpy {

namespace {
const double CW = SHCLIQ;   // J kg-1 K-1 (physicalconst.h, same values as layerEnthalpy)
const double CI = SHCICE;
const double LF = LHFUS;    // J kg-1

// The SJ2016 curve is the shared sjcurve:: implementation (SJCurve.cpp), the same code
// Soil_Bgc::sjLayerPhi calls; tests/unit checks bitwise equality.
double sjDPhiPure(const double t, const double b) {   // d/dt of sjPhiPure for t <= 0
  if (t > 0.0) { return 0.0; }
  return -b / 0.01 * pow((0.1 - t) / 0.01, b - 1.0);
}
bool isFree(const Cell& c) { return c.kind != SOIL || c.sfcc == SFCC_FREE; }
double cth(const Cell& c) { return c.Cs + c.W * CW; }          // thawed slope
// enthalpy at a given liquid fraction
double Hs(const Cell& c, double T, double phi) {
  return c.Cs * T + c.W * phi * CW * T + c.W * (1.0 - phi) * (CI * T - LF);
}
double sjDPhi(const Cell& c, double t) {
  const double fo = fmin(1.0, fmax(0.0, c.forg));
  const double cc = fmax(0.0, c.clay), si = fmax(0.0, c.silt), sa = fmax(0.0, c.sand);
  const double sum = cc + si + sa;
  const double dpo = sjDPhiPure(t, -1.0);
  if (!(sum > 0.0)) { return dpo; }
  const double dpm = (cc * sjDPhiPure(t, -0.3) + si * sjDPhiPure(t, -0.5) + sa * sjDPhiPure(t, -0.9)) / sum;
  return (1.0 - fo) * dpm + fo * dpo;
}
// SFCC 1 liquid fraction and derivative
double phi1(const Cell& c, double T, double* dphi) {
  const double pm = phimin(c);
  if (T > 0.0) { if (dphi) *dphi = 0.0; return 1.0; }
  if (T >= -c.eps) {
    const double pe = fmax(sjPhi(-c.eps, c.forg, c.clay, c.silt, c.sand), pm);
    if (dphi) *dphi = (1.0 - pe) / c.eps;
    return pe + (1.0 - pe) * (T + c.eps) / c.eps;
  }
  const double s = sjPhi(T, c.forg, c.clay, c.silt, c.sand);
  if (s > pm) { if (dphi) *dphi = sjDPhi(c, T); return s; }
  if (dphi) *dphi = 0.0;
  return pm;
}
double Hm0(const Cell& c) { return Hs(c, -c.eps, phimin(c)); } // sfcc0: H at -eps
} // namespace

Cell::Cell() : kind(SOIL), dz(0.1), Cs(0.0), W(0.0), liqmin(0.0), sfcc(SFCC_FREE), eps(0.01),
  clay(0.0), silt(0.0), sand(0.0), forg(1.0), lam_u(1.0), lam_f(1.0) {}
BC::BC() : top_type(1), Ts(0.0), cn_top(100.0), q_bot(0.0) {}
Options::Options() : face_mode(0), tol(1.0e-3), max_outer(50), max_inner(50), max_katz(2000) {}
Stats::Stats() : outer(0), inner(0), solves(0), converged(false), res_max(0), dH(0), gtop(0),
  gbot(0), src(0), closure(0) {}

double sjPhi(double t, double forg, double fc, double fsi, double fsa) {
  return sjcurve::layerPhi(t, forg, fc, fsi, fsa);
}

double phimin(const Cell& c) {
  if (c.kind == ROCK || !(c.W > 0.0)) { return 1.0; }
  return fmin(1.0, fmax(0.0, c.liqmin) / c.W);
}

bool hasLatent(const Cell& c) { return c.kind != ROCK && c.W > 0.0 && phimin(c) < 1.0; }

double H(const Cell& c, double T) {
  if (!hasLatent(c)) { return cth(c) * T; }
  if (T > 0.0) { return cth(c) * T; }
  if (isFree(c)) {
    const double hm = Hm0(c);
    if (T >= -c.eps) { return hm * (-T / c.eps); }   // linear H on the ramp (piecewise affine)
    return Hs(c, T, phimin(c));
  }
  return Hs(c, T, phi1(c, T, 0));
}

double Ca(const Cell& c, double T) {
  if (!hasLatent(c) || T > 0.0) { return cth(c); }
  const double pm = phimin(c);
  if (isFree(c)) {
    if (T > -c.eps) { return -Hm0(c) / c.eps; }
    return c.Cs + c.W * pm * CW + c.W * (1.0 - pm) * CI;
  }
  double dphi = 0.0;
  const double p = phi1(c, T, &dphi);
  return c.Cs + c.W * p * CW + c.W * (1.0 - p) * CI + c.W * dphi * ((CW - CI) * T + LF);
}

double liq(const Cell& c, double T) {
  if (c.kind == ROCK) { return 0.0; }
  if (!hasLatent(c) || T > 0.0) { return c.W; }
  if (isFree(c)) {
    if (T >= -c.eps) {   // liquid consistent with the linear H: H = Hs(T, liq/W)
      const double h = H(c, T);
      const double l = (h - c.Cs * T - c.W * (CI * T - LF)) / ((CW - CI) * T + LF);
      return fmin(c.W, fmax(c.W * phimin(c), l));
    }
    return c.W * phimin(c);
  }
  return c.W * phi1(c, T, 0);
}

double conductivity(const Cell& c, double T) {
  if (c.kind == ROCK || !(c.W > 0.0)) { return c.lam_u; }
  const double f = fmin(1.0, fmax(0.0, 1.0 - liq(c, T) / c.W));   // ice fraction
  if (f <= 0.0) { return c.lam_u; }
  if (f >= 1.0) { return c.lam_f; }
  return pow(c.lam_u, 1.0 - f) * pow(c.lam_f, f);
}

double invertH(const Cell& c, double h) {
  if (!hasLatent(c) || h >= 0.0) {
    const double k = cth(c);
    return k > 0.0 ? h / k : 0.0;
  }
  // bracket [lo, hi] with H(lo) <= h <= H(hi), H monotone increasing
  double hi = 0.0, lo = -c.eps;
  if (H(c, lo) > h) { hi = lo; lo = -1.0; while (H(c, lo) > h) { hi = lo; lo *= 2.0; if (lo < -1.0e4) break; } }
  double T = 0.5 * (lo + hi);
  for (int it = 0; it < 200; ++it) {
    const double f = H(c, T) - h;
    if (f > 0.0) { hi = T; } else { lo = T; }
    const double d = Ca(c, T);
    double Tn = (d > 0.0) ? T - f / d : 0.5 * (lo + hi);
    if (!(Tn > lo && Tn < hi)) { Tn = 0.5 * (lo + hi); }
    if (fabs(Tn - T) <= 1.0e-15 * (1.0 + fabs(T)) || hi - lo <= 1.0e-15 * (1.0 + fabs(T))) { return Tn; }
    T = Tn;
  }
  return T;
}

// ---------------------------------------------------------------------------
ColumnSolver::ColumnSolver() : n(0), dt_(0.0), ktop_(0.0) {}

void ColumnSolver::assemble(const std::vector<Cell>& cells, const std::vector<double>& T,
                            double dt, const BC& bc, const std::vector<double>* S) {
  n = (int)cells.size();
  K.assign(n + 1, 0.0); lam.assign(n, 0.0); lo.assign(n, 0.0); up.assign(n, 0.0);
  di.assign(n, 0.0); b.assign(n, 0.0); Hn.assign(n, 0.0); rhs.assign(n, 0.0);
  w1.assign(n, 0.0); w2.assign(n, 0.0); x.assign(n, 0.0); v.assign(n, 0.0); r.assign(n, 0.0);
  dt_ = dt;
  for (int i = 0; i < n; ++i) { lam[i] = conductivity(cells[i], T[i]); }
  const double half0 = 0.5 * cells[0].dz / lam[0];
  K[0] = (bc.top_type == 0) ? 1.0 / half0 : 1.0 / (1.0 / bc.cn_top + half0);
  for (int i = 1; i < n; ++i) {
    if (opt.face_mode == 1) {   // Tubini et al. 2021 Eq. 13 (max rule) - optional fallback
      K[i] = std::max(lam[i - 1], lam[i]) / (0.5 * cells[i - 1].dz + 0.5 * cells[i].dz);
    } else {                    // harmonic / series: exact steady flux across layered media
      K[i] = 1.0 / (0.5 * cells[i - 1].dz / lam[i - 1] + 0.5 * cells[i].dz / lam[i]);
    }
  }
  K[n] = 0.0;   // bottom: prescribed flux
  ktop_ = K[0];
  for (int i = 0; i < n; ++i) {
    di[i] = dt * (K[i] + K[i + 1]);
    lo[i] = (i > 0) ? -dt * K[i] : 0.0;
    up[i] = (i < n - 1) ? -dt * K[i + 1] : 0.0;
    Hn[i] = H(cells[i], T[i]);
    b[i] = Hn[i] + (S ? dt * (*S)[i] : 0.0);
  }
  b[0] += dt * K[0] * bc.Ts;
  b[n - 1] += dt * bc.q_bot;
}

void ColumnSolver::applyA(const std::vector<double>& T, std::vector<double>& out) const {
  for (int i = 0; i < n; ++i) {
    double s = di[i] * T[i];
    if (i > 0) s += lo[i] * T[i - 1];
    if (i < n - 1) s += up[i] * T[i + 1];
    out[i] = s;
  }
}

void ColumnSolver::thomas(const std::vector<double>& d, std::vector<double>& sol) {
  // (A + diag(d)) sol = rhs ; w1 = modified upper, w2 = modified rhs
  double m = di[0] + d[0];
  w1[0] = up[0] / m; w2[0] = rhs[0] / m;
  for (int i = 1; i < n; ++i) {
    m = di[i] + d[i] - lo[i] * w1[i - 1];
    w1[i] = (i < n - 1) ? up[i] / m : 0.0;
    w2[i] = (rhs[i] - lo[i] * w2[i - 1]) / m;
  }
  sol[n - 1] = w2[n - 1];
  for (int i = n - 2; i >= 0; --i) { sol[i] = w2[i] - w1[i] * sol[i + 1]; }
}

void ColumnSolver::finish(const std::vector<Cell>& cells, const std::vector<double>& T,
                          double dt, const BC& bc, const std::vector<double>* S, Stats& st) {
  applyA(T, r);
  st.res_max = 0.0; st.dH = 0.0; st.src = 0.0;
  for (int i = 0; i < n; ++i) {
    const double h = H(cells[i], T[i]);
    st.res_max = std::max(st.res_max, fabs(h + r[i] - b[i]));
    st.dH += h - Hn[i];
    if (S) st.src += dt * (*S)[i];
  }
  st.gtop = dt * ktop_ * (bc.Ts - T[0]);
  st.gbot = dt * bc.q_bot;
  st.closure = st.dH - (st.gtop + st.gbot + st.src);
}

// round-off floor of a residual H(T) + A T - b: a few ulp of the largest term
double ColumnSolver::rfloor(const std::vector<double>& AT, const std::vector<double>& T) const {
  double m = 0.0;
  for (int i = 0; i < n; ++i) { m = std::max(m, std::max(fabs(b[i]), fabs(AT[i]) + fabs(di[i] * T[i]))); }
  return 8.0 * 2.2e-16 * m;
}

bool ColumnSolver::stepNCZ(const std::vector<Cell>& cells, std::vector<double>& T, double dt,
                           const BC& bc, const std::vector<double>* S, Stats& st) {
  st = Stats();
  assemble(cells, T, dt, bc, S);
  std::vector<double> ca0(n), q(n), p(n), dg(n), y(n), yn(n), Ay(n);
  for (int i = 0; i < n; ++i) { ca0[i] = Ca(cells[i], 0.0); x[i] = std::min(T[i], 0.0); }  // T0 <= T* = 0
  bool ok = false;
  for (int k = 0; k < opt.max_outer && !ok; ++k) {
    st.outer++;
    // outer linearisation of h2 at x (Tubini Eq. 21): h1(T) - [h2(x) + q (T - x)] + A T = b;
    // h2 = (Ca(0-) - Cth) T for T > 0 is linear, so h2(x) - q x = 0 and the right side stays b
    for (int i = 0; i < n; ++i) { q[i] = (x[i] > 0.0) ? ca0[i] - cth(cells[i]) : 0.0; }
    y = x;
    for (int l = 0; l < opt.max_inner; ++l) {
      applyA(y, Ay);
      double rin = 0.0;
      for (int i = 0; i < n; ++i) {
        const double h1 = (y[i] > 0.0) ? ca0[i] * y[i] : H(cells[i], y[i]);
        rin = std::max(rin, fabs(h1 - q[i] * y[i] + Ay[i] - b[i]));
        p[i] = (y[i] > 0.0) ? ca0[i] : Ca(cells[i], y[i]);
        rhs[i] = b[i] - h1 + p[i] * y[i];
        dg[i] = p[i] - q[i];
      }
      if (rin < opt.tol + rfloor(Ay, y)) break;
      thomas(dg, yn);
      st.inner++;
      double dmax = 0.0;
      for (int i = 0; i < n; ++i) { dmax = std::max(dmax, fabs(yn[i] - y[i]) / (1.0 + fabs(y[i]))); }
      y.swap(yn);
      if (dmax < 1.0e-15) break;   // round-off floor
    }
    x = y;
    applyA(x, Ay);
    double rout = 0.0;
    for (int i = 0; i < n; ++i) { rout = std::max(rout, fabs(H(cells[i], x[i]) + Ay[i] - b[i])); }
    if (rout < opt.tol + rfloor(Ay, x)) ok = true;
  }
  st.solves = st.inner;
  T = x;
  finish(cells, T, dt, bc, S, st);
  st.converged = ok;
  return ok;
}

bool ColumnSolver::stepKatz(const std::vector<Cell>& cells, std::vector<double>& T, double dt,
                            const BC& bc, const std::vector<double>* S, Stats& st) {
  st = Stats();
  assemble(cells, T, dt, bc, S);
  std::vector<int> seg(n), lat(n);
  std::vector<double> sl0(n), sl1(n), sl2(n), dg(n), Ax(n);
  for (int i = 0; i < n; ++i) {
    const Cell& c = cells[i];
    lat[i] = hasLatent(c) ? 1 : 0;
    if (lat[i] && !isFree(c)) { return false; }   // Katzenelson needs piecewise-affine H
    sl2[i] = cth(c);
    if (lat[i]) { sl1[i] = -Hm0(c) / c.eps; sl0[i] = Ca(c, -2.0 * c.eps); } else { sl0[i] = sl1[i] = sl2[i]; }
    x[i] = T[i];
    seg[i] = !lat[i] ? 2 : (x[i] > 0.0 ? 2 : (x[i] >= -c.eps ? 1 : 0));
  }
  bool ok = false;
  for (int it = 0; it <= opt.max_katz; ++it) {
    applyA(x, Ax);
    double rm = 0.0;
    for (int i = 0; i < n; ++i) { r[i] = H(cells[i], x[i]) + Ax[i] - b[i]; rm = std::max(rm, fabs(r[i])); }
    if (rm < opt.tol + rfloor(Ax, x)) { ok = true; break; }
    if (it == opt.max_katz) break;
    for (int i = 0; i < n; ++i) {
      dg[i] = (seg[i] == 0) ? sl0[i] : (seg[i] == 1 ? sl1[i] : sl2[i]);
      rhs[i] = -r[i];
    }
    thomas(dg, v);   // Newton direction in the current affine region
    st.solves++;
    double tmax = 1.0; int imin = -1, nseg = 0; double bnd = 0.0;
    for (int i = 0; i < n; ++i) {
      if (!lat[i] || v[i] == 0.0) continue;
      const double e = cells[i].eps;
      if (v[i] < 0.0 && seg[i] > 0) {
        const double lowb = (seg[i] == 1) ? -e : 0.0;
        const double t = (lowb - x[i]) / v[i];
        if (t < tmax) { tmax = std::max(t, 0.0); imin = i; nseg = seg[i] - 1; bnd = lowb; }
      } else if (v[i] > 0.0 && seg[i] < 2) {
        const double hib = (seg[i] == 0) ? -e : 0.0;
        const double t = (hib - x[i]) / v[i];
        if (t < tmax) { tmax = std::max(t, 0.0); imin = i; nseg = seg[i] + 1; bnd = hib; }
      }
    }
    for (int i = 0; i < n; ++i) { x[i] += tmax * v[i]; }
    if (imin >= 0) { x[imin] = bnd; seg[imin] = nseg; }
  }
  st.outer = st.solves; st.inner = st.solves;
  T = x;
  finish(cells, T, dt, bc, S, st);
  st.converged = ok;
  return ok;
}

// ---------------------------------------------------------------------------
// Stage C helpers
double Hstate(const Cell& c, double T, double liqv, double icev) {
  return c.Cs * T + liqv * CW * T + icev * (CI * T - LF);
}

bool onCurve(const Cell& c, double T, double liqv, double tol_kg) {
  return fabs(liqv - liq(c, T)) <= tol_kg;
}

double equilibrate(const Cell& c, double h, double& liqv, double& icev) {
  const double T = invertH(c, h);
  liqv = liq(c, T);
  icev = fmax(0.0, c.W - liqv);
  return T;
}

double mergeT(double Cu, double Tu, double Cl, double Tl) {
  const double C = Cu + Cl;
  return (C > 0.0) ? (Cu * Tu + Cl * Tl) / C : 0.5 * (Tu + Tl);
}

void reconstructFronts(const FrontIn& in, FrontOut& out) {
  const int n = (int)in.T.size();
  out.fz.clear(); out.ftype.clear();
  out.frozen.assign(n, 1); out.frozenfrac.assign(n, 1.0);
  std::vector<int> st(n); std::vector<double> f(n, 1.0);
  for (int i = 0; i < n; ++i) {
    if (!in.latent[i]) { st[i] = (in.T[i] > 0.0) ? -1 : 1; f[i] = (st[i] == 1) ? 1.0 : 0.0; continue; }
    const double fi = (in.ice_eps[i] > 0.0) ? fmin(1.0, fmax(0.0, in.ice[i] / in.ice_eps[i])) : 0.0;
    f[i] = fi;
    if (in.ice[i] <= 1.0e-12 * fmax(in.W[i], 1.0) || fi <= 1.0e-12) { st[i] = -1; f[i] = 0.0; }
    else if (in.T[i] < -in.eps || fi >= 1.0 - 1.0e-12) { st[i] = 1; f[i] = 1.0; }
    else { st[i] = 0; }
  }
  std::vector<double> fz; std::vector<int> ft;
  int prev_state = (in.Tabove > 0.0) ? -1 : 1, prev_idx = -1, ra = -1;
  for (int i = 0; i <= n; ++i) {
    const bool end = (i == n);
    if (!end && st[i] == 0) { if (ra < 0) ra = i; continue; }
    const int s = end ? ((in.Tbelow > 0.0) ? -1 : 1) : st[i];
    if (ra >= 0) {   // ramp run [ra, i-1] between prev_state (above) and s (below)
      const int rb = i - 1;
      const double ztop = in.z[ra], zbot = in.z[rb] + in.dz[rb];
      double F = 0.0, U = 0.0;
      for (int k = ra; k <= rb; ++k) { F += f[k] * in.dz[k]; U += (1.0 - f[k]) * in.dz[k]; }
      const double zc = 0.5 * (ztop + zbot), lo = ztop + 1.0e-9 * (zbot - ztop);
      if (prev_state == -1 && s == 1)      { fz.push_back(fmax(lo, zbot - F)); ft.push_back(-1); }
      else if (prev_state == 1 && s == -1) { fz.push_back(fmin(zbot, fmax(lo, ztop + F))); ft.push_back(1); }
      else if (prev_state == -1 && s == -1) {
        if (F > 0.0) { fz.push_back(fmax(lo, zc - 0.5 * F)); ft.push_back(-1); fz.push_back(fmin(zbot, zc + 0.5 * F)); ft.push_back(1); }
      } else {
        if (U > 0.0) { fz.push_back(fmax(lo, zc - 0.5 * U)); ft.push_back(1); fz.push_back(fmin(zbot, zc + 0.5 * U)); ft.push_back(-1); }
      }
      ra = -1;
    } else if (!end && prev_idx >= 0 && prev_state != s) {   // sharp crossing between two cells
      const int p = prev_idx;
      const double zp = in.z[p] + 0.5 * in.dz[p], zi = in.z[i] + 0.5 * in.dz[i];
      const double d = in.T[p] - in.T[i];
      double w = (d != 0.0) ? in.T[p] / d : 0.5;
      w = fmin(1.0, fmax(0.0, w));
      double z0 = zp + w * (zi - zp);
      const double zt = in.z[p], zb = in.z[i] + in.dz[i];
      z0 = fmin(zb, fmax(zt + 1.0e-9 * (zb - zt), z0));
      fz.push_back(z0); ft.push_back(prev_state == 1 ? 1 : -1);
    }
    if (!end) { prev_state = s; prev_idx = i; }
  }
  // cap the number of fronts: remove the closest adjacent pair
  while ((int)fz.size() > in.maxfronts && fz.size() >= 2) {
    size_t kmin = 0; double dmin = 1.0e300;
    for (size_t k = 0; k + 1 < fz.size(); ++k) { if (fz[k + 1] - fz[k] < dmin) { dmin = fz[k + 1] - fz[k]; kmin = k; } }
    fz.erase(fz.begin() + kmin, fz.begin() + kmin + 2); ft.erase(ft.begin() + kmin, ft.begin() + kmin + 2);
  }
  out.fz = fz; out.ftype = ft;
  // cell states from the fronts (legacy convention, cf. Ground::getLayerFrozenstatusByFronts)
  for (int i = 0; i < n; ++i) {
    const double zt = in.z[i], zb = in.z[i] + in.dz[i];
    double pprev = zt, frz = 0.0; int nin = 0, last = 0;
    for (size_t k = 0; k < fz.size(); ++k) {
      if (fz[k] > zt && fz[k] <= zb) {
        if (ft[k] == 1) frz += fz[k] - pprev;
        pprev = fz[k]; last = ft[k]; nin++;
      }
    }
    if (nin > 0) {
      if (last == -1) frz += zb - pprev;
      out.frozen[i] = 0; out.frozenfrac[i] = fmin(1.0, fmax(0.0, frz / in.dz[i]));
      continue;
    }
    int s = st[i];
    int above = 0, below = 0;
    for (size_t k = 0; k < fz.size(); ++k) { if (fz[k] <= zt) above = ft[k]; }
    for (size_t k = fz.size(); k-- > 0;) { if (fz[k] > zb) below = ft[k]; }
    if (above != 0) s = -above; else if (below != 0) s = below;
    else if (s == 0) s = (f[i] >= 0.5) ? 1 : -1;
    out.frozen[i] = s; out.frozenfrac[i] = (s == 1) ? 1.0 : 0.0;
  }
}

} // namespace enthalpy
