// Stage B: standalone enthalpy solver (include/EnthalpySolver.h), NOT used by the model.
// Design: /workspace/tem-dev/enthalpy_solver_design.md. Refs: Tubini et al. 2021,
// doi:10.5194/tc-15-2541-2021 (NCZ, Neumann solution App. E); Casulli & Zanolli 2010,
// doi:10.1137/100786320; Hoetten et al. 2026, doi:10.1553/etna_vol65s553 (Katzenelson);
// Schaefer & Jafarov 2016, doi:10.5194/bg-13-1991-2016 (unfrozen-water curve).
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/EnthalpySolver.h"
#include <chrono>
#include <cmath>
#include <vector>

using namespace enthalpy;

namespace {
const double DAY = 86400.0;

Cell soil(double dz, double poro, double vhc, double sat, int sfcc, double clay, double silt,
          double sand, double forg, double lu, double lf) {
  Cell c; c.kind = SOIL; c.dz = dz; c.Cs = vhc * (1.0 - poro) * dz; c.W = sat * poro * 1000.0 * dz;
  c.liqmin = 0.05 * poro * 1000.0 * dz; c.sfcc = sfcc; c.eps = 0.01;
  c.clay = clay; c.silt = silt; c.sand = sand; c.forg = forg; c.lam_u = lu; c.lam_f = lf; return c;
}
// dvmdostem-like column: 2 snow, moss, fibric, 2 humic, 15 mineral (MINETHICK), 5 rock (ROCKTHICK)
std::vector<Cell> column(int sfcc, bool snow = true) {
  std::vector<Cell> v;
  if (snow) {
    for (int i = 0; i < 2; ++i) {
      Cell s; s.kind = SNOW; s.dz = 0.15; s.W = 250.0 * 0.15; s.lam_u = s.lam_f = 0.25; s.eps = 0.01; v.push_back(s);
    }
  }
  v.push_back(soil(0.04, 0.98, 2.5e6, 0.3, sfcc, 0, 0, 0, 1.0, 0.10, 0.25));
  v.push_back(soil(0.08, 0.95, 2.5e6, 0.5, sfcc, 0, 0, 0, 1.0, 0.20, 0.60));
  v.push_back(soil(0.10, 0.85, 2.5e6, 0.8, sfcc, 0, 0, 0, 1.0, 0.35, 1.20));
  v.push_back(soil(0.10, 0.85, 2.5e6, 0.9, sfcc, 0, 0, 0, 1.0, 0.40, 1.40));
  const double mt[15] = {0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.5,0.5,1.0,1.0};
  for (int i = 0; i < 15; ++i) v.push_back(soil(mt[i], 0.45, 2.0e6, 0.9, sfcc, 20, 60, 20, 0.05, 1.3, 2.1));
  const double rt[5] = {2.0, 4.0, 8.0, 16.0, 20.0};
  for (int i = 0; i < 5; ++i) { Cell r; r.kind = ROCK; r.dz = rt[i]; r.Cs = 2.0e6 * rt[i]; r.lam_u = r.lam_f = 2.0; v.push_back(r); }
  return v;
}
double tsurf(int d) { return -10.0 + 20.0 * sin(2.0 * M_PI * (d - 110) / 365.0); }
// depth of the deepest 0 degC crossing below a thawed top (thaw depth), cell-centre interpolation
double zeroIsoDepth(const std::vector<Cell>& c, const std::vector<double>& T, int first) {
  double z = 0.0, zc_prev = 0.0, Tp = 0.0; bool have = false; double out = 0.0;
  for (size_t i = first; i < c.size(); ++i) {
    const double zc = z + 0.5 * c[i].dz;
    if (have && ((Tp > 0.0) != (T[i] > 0.0))) { out = zc_prev + (zc - zc_prev) * Tp / (Tp - T[i]); }
    have = true; Tp = T[i]; zc_prev = zc; z += c[i].dz;
  }
  return out;
}
double frozenDepth(const std::vector<Cell>& c, const std::vector<double>& T) {   // sum ice-fraction * dz
  double z = 0.0;
  for (size_t i = 0; i < c.size(); ++i) {
    if (!(c[i].W > 0.0)) continue;
    const double fr = (c[i].W - liq(c[i], T[i])) / (c[i].W - c[i].W * phimin(c[i]));
    z += c[i].dz * std::min(1.0, std::max(0.0, fr));
  }
  return z;
}
// Neumann (Stefan) solution for freezing of water, Tubini et al. 2021 App. E
double neumannGamma(double li, double lw, double ci, double cw, double rho, double L, double Ts, double Tm, double T0) {
  const double ai = li / (rho * ci), aw = lw / (rho * cw);
  auto f = [&](double g) {
    const double gp = g * sqrt(ai / aw);
    return li * (Tm - Ts) * exp(-g * g) / (erf(g) * sqrt(M_PI * ai))
         - lw * (T0 - Tm) * exp(-gp * gp) / (erfc(gp) * sqrt(M_PI * aw)) - rho * L * g * sqrt(ai);
  };
  double a = 1e-6, b = 5.0;
  for (int i = 0; i < 200; ++i) { const double m = 0.5 * (a + b); if (f(m) > 0) a = m; else b = m; }
  return 0.5 * (a + b);
}
} // namespace

TEST_CASE("enthalpy SFCC 1: sjPhi is bitwise identical to Soil_Bgc::sjLayerPhi; liq = W*max(phi, phimin) below -eps") {
  const double tex[4][4] = {{20, 60, 20, 0.0}, {0, 0, 0, 1.0}, {0.4, 0.3, 0.3, 0.3}, {1, 0, 0, 0.7}};
  for (int k = 0; k < 4; ++k) {
    for (double t = -40.0; t < 2.0; t += 0.0137) {
      CHECK(sjPhi(t, tex[k][3], tex[k][0], tex[k][1], tex[k][2]) == Soil_Bgc::sjLayerPhi(t, tex[k][3], tex[k][0], tex[k][1], tex[k][2]));
    }
    Cell c = soil(0.1, 0.5, 2e6, 0.9, SFCC_SJ2016, tex[k][0], tex[k][1], tex[k][2], tex[k][3], 1, 2);
    for (double t = -30.0; t < -0.0101; t += 0.031) {
      const double ref = std::max(Soil_Bgc::sjLayerPhi(t, tex[k][3], tex[k][0], tex[k][1], tex[k][2]), phimin(c));
      CHECK(liq(c, t) / c.W == doctest::Approx(ref).epsilon(1e-14));
    }
    CHECK(liq(c, 0.0) == doctest::Approx(c.W));
  }
}

TEST_CASE("enthalpy curve: continuous, monotone, C1/C2 (Ca >= 0, non-decreasing to T* = 0), invertH round trip") {
  for (int sf = 0; sf < 2; ++sf) {
    std::vector<Cell> col = column(sf);
    for (size_t i = 0; i < col.size(); ++i) {
      const Cell& c = col[i];
      double hp = H(c, -50.0), cap = Ca(c, -50.0);
      for (double t = -50.0 + 7e-4; t < 20.0; t += (t > -0.05 && t < 0.01) ? 1e-5 : 7e-3) {
        const double h = H(c, t), a = Ca(c, t);
        CHECK(h >= hp);
        CHECK(a >= 0.0);
        if (t <= 0.0 && hasLatent(c)) { CHECK(a >= cap * (1.0 - 1e-12)); }   // non-decreasing up to T*
        if (t > 0.0) { CHECK(a <= Ca(c, 0.0)); }                            // non-increasing above
        hp = h; cap = a;
      }
      for (double d : {1e-9, 1e-12}) {   // continuity at the kinks: jump bounded by slope * width
        const double cmax = std::max(Ca(c, 0.0), Ca(c, -c.eps - d));
        CHECK(fabs(H(c, -c.eps + d) - H(c, -c.eps - d)) <= 2.0 * d * cmax * 1.0001 + 1e-6);
        CHECK(fabs(H(c, d) - H(c, -d)) <= 2.0 * d * cmax * 1.0001 + 1e-6);
      }
      for (double t = -40.0; t < 20.0; t += 0.0071) { CHECK(invertH(c, H(c, t)) == doctest::Approx(t).epsilon(1e-10).scale(1.0)); }
      CHECK(fabs(H(c, 0.0)) < 1e-6);
    }
  }
}

TEST_CASE("enthalpy NCZ: Neumann analytical freezing front (water, dz 5 mm, dt 1 h)") {
  const int N = 400; const double dz = 0.005, Ts = -10.0, T0 = 3.0, rho = 1000.0;
  std::vector<Cell> col(N);
  for (int i = 0; i < N; ++i) { col[i].kind = SOIL; col[i].dz = dz; col[i].W = rho * dz; col[i].liqmin = 0;
    col[i].sfcc = SFCC_FREE; col[i].eps = 1e-4; col[i].lam_u = TCLIQ; col[i].lam_f = TCICE; }
  std::vector<double> T(N, T0);
  ColumnSolver s; BC bc; bc.top_type = 0; bc.Ts = Ts; Stats st;
  const double g = neumannGamma(TCICE, TCLIQ, SHCICE, SHCLIQ, rho, LHFUS, Ts, 0.0, T0);
  const double ai = TCICE / (rho * SHCICE);
  double maxerr = 0.0;
  for (int h = 1; h <= 240; ++h) {
    REQUIRE(s.stepNCZ(col, T, 3600.0, bc, 0, st));
    const double X = 2.0 * g * sqrt(ai * h * 3600.0);
    if (h >= 24) maxerr = std::max(maxerr, fabs(frozenDepth(col, T) - X));
  }
  MESSAGE("Neumann: gamma=" << g << " front after 10 d=" << frozenDepth(col, T) << " m, max |error| (days 1-10)=" << maxerr << " m");
  CHECK(maxerr < 2.0 * dz);
}

TEST_CASE("enthalpy NCZ: energy conservation <= 1e-6 W m-2 per daily step for a year (sfcc 0 and 1, with sources)") {
  for (int sf = 0; sf < 2; ++sf) {
    std::vector<Cell> col = column(sf);
    std::vector<double> T(col.size(), -2.0), S(col.size(), 0.0);
    ColumnSolver s; BC bc; Stats st; double worst = 0.0, cum = 0.0;
    for (int d = 0; d < 365; ++d) {
      bc.Ts = tsurf(d); bc.q_bot = 0.06; S[0] = (d % 7 == 0) ? 15.0 : 0.0;
      REQUIRE(s.stepNCZ(col, T, DAY, bc, &S, st));
      worst = std::max(worst, fabs(st.closure) / DAY); cum += st.closure;
    }
    MESSAGE("sfcc " << sf << ": max |daily closure| = " << worst << " W m-2, annual mean = " << fabs(cum) / (365 * DAY));
    CHECK(worst <= 1e-6);
  }
}

TEST_CASE("enthalpy NCZ: converges for large dt (1 h, 1 d, 10 d); 1 d thaw depth within one cell of 1 h") {
  std::vector<Cell> col = column(SFCC_SJ2016, false);
  const double dts[3] = {3600.0, DAY, 10 * DAY};
  double zt[3];
  for (int k = 0; k < 3; ++k) {
    std::vector<double> T(col.size(), -3.0);
    ColumnSolver s; BC bc; bc.Ts = 8.0; Stats st;
    const int nst = (int)(120 * DAY / dts[k]);
    for (int i = 0; i < nst; ++i) { REQUIRE(s.stepNCZ(col, T, dts[k], bc, 0, st)); CHECK(fabs(st.closure) / dts[k] <= 1e-6); }
    zt[k] = zeroIsoDepth(col, T, 0);
  }
  MESSAGE("thaw depth after 120 d: 1h " << zt[0] << " 1d " << zt[1] << " 10d " << zt[2]);
  CHECK(fabs(zt[1] - zt[0]) < 0.3);   // mineral cells there are 0.2-0.3 m
  CHECK(fabs(zt[2] - zt[0]) < 0.5);
}

TEST_CASE("enthalpy NCZ: freezing front is monotone under constant cold surface") {
  for (int sf = 0; sf < 2; ++sf) {
    std::vector<Cell> col = column(sf, false);
    std::vector<double> T(col.size(), 2.0);
    ColumnSolver s; BC bc; bc.Ts = -15.0; Stats st; double prev = 0.0;
    for (int d = 0; d < 200; ++d) {
      REQUIRE(s.stepNCZ(col, T, DAY, bc, 0, st));
      const double z = frozenDepth(col, T);
      CHECK(z >= prev - 1e-12);
      prev = z;
    }
    CHECK(prev > 0.5);
  }
}

TEST_CASE("enthalpy faces: harmonic mean gives the exact two-layer steady flux (Dirichlet and Robin top)") {
  std::vector<Cell> col(2);
  col[0].kind = ROCK; col[0].dz = 0.2; col[0].Cs = 2e6 * 0.2; col[0].lam_u = col[0].lam_f = 0.2;
  col[1].kind = ROCK; col[1].dz = 0.6; col[1].Cs = 2e6 * 0.6; col[1].lam_u = col[1].lam_f = 2.0;
  const double q = 0.3, Ts = 5.0;
  for (int top = 0; top < 2; ++top) {
    ColumnSolver s; BC bc; bc.top_type = top; bc.Ts = Ts; bc.q_bot = q; bc.cn_top = 100.0; Stats st;
    std::vector<double> T(2, 0.0);
    for (int i = 0; i < 5; ++i) REQUIRE(s.stepNCZ(col, T, 1e12, bc, 0, st));
    const double r0 = (top == 1) ? 1.0 / 100.0 : 0.0;
    CHECK(T[0] == doctest::Approx(Ts + q * (r0 + 0.1 / 0.2)).epsilon(1e-8));
    CHECK(T[1] == doctest::Approx(Ts + q * (r0 + 0.2 / 0.2 + 0.3 / 2.0)).epsilon(1e-8));
    s.opt.face_mode = 1;   // documented fallback: Tubini max rule overestimates the flux through the insulating layer
    for (int i = 0; i < 5; ++i) REQUIRE(s.stepNCZ(col, T, 1e12, bc, 0, st));
    CHECK(T[1] < Ts + q * (r0 + 0.2 / 0.2 + 0.3 / 2.0) - 0.05);
  }
}

TEST_CASE("enthalpy Katzenelson (free water) agrees with NCZ; refuses SFCC 1") {
  unsigned seed = 12345;
  auto rnd = [&]() { seed = seed * 1103515245u + 12345u; return ((seed >> 8) & 0xffff) / 65535.0; };
  int nk = 0, nn = 0;
  for (int trial = 0; trial < 200; ++trial) {
    std::vector<Cell> col = column(SFCC_FREE, trial % 2 == 0);
    std::vector<double> T0(col.size());
    for (size_t i = 0; i < col.size(); ++i) T0[i] = (rnd() - 0.6) * 6.0 * (i < 8 ? 1.0 : 0.3);
    if (trial % 5 == 0) for (size_t i = 3; i < 10; ++i) T0[i] = -0.005;   // cells on the ramp
    BC bc; bc.Ts = (rnd() - 0.5) * 40.0; bc.q_bot = 0.05;
    const double dt = (trial % 3 == 0) ? 10 * DAY : DAY;
    std::vector<double> Tn = T0, Tk = T0; ColumnSolver a, b; Stats sa, sb;
    REQUIRE(a.stepNCZ(col, Tn, dt, bc, 0, sa));
    REQUIRE(b.stepKatz(col, Tk, dt, bc, 0, sb));
    for (size_t i = 0; i < col.size(); ++i) CHECK(Tk[i] == doctest::Approx(Tn[i]).epsilon(1e-8).scale(1.0));
    CHECK(fabs(sb.closure) / dt <= 1e-6);
    nk += sb.solves; nn += sa.solves;
  }
  MESSAGE("random columns: mean solves/step NCZ " << nn / 200.0 << ", Katzenelson " << nk / 200.0);
  std::vector<Cell> c1 = column(SFCC_SJ2016); std::vector<double> T(c1.size(), -1.0);
  ColumnSolver k; BC bc; Stats st;
  CHECK_FALSE(k.stepKatz(c1, T, DAY, bc, 0, st));
}

TEST_CASE("enthalpy microbenchmark: solves per day and us per column-day on a synthetic 31-cell column (3 years)") {
  for (int mode = 0; mode < 3; ++mode) {   // 0 NCZ sfcc0, 1 NCZ sfcc1, 2 Katzenelson sfcc0
    std::vector<Cell> col = column(mode == 1 ? SFCC_SJ2016 : SFCC_FREE);
    std::vector<double> T(col.size(), -2.0);
    ColumnSolver s; BC bc; Stats st; long solves = 0; int mx = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (int d = 0; d < 3 * 365; ++d) {
      bc.Ts = tsurf(d);
      const bool ok = (mode == 2) ? s.stepKatz(col, T, DAY, bc, 0, st) : s.stepNCZ(col, T, DAY, bc, 0, st);
      REQUIRE(ok);
      solves += st.solves; mx = std::max(mx, st.solves);
    }
    const double us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / (3 * 365);
    MESSAGE(std::string(mode == 0 ? "NCZ sfcc0" : mode == 1 ? "NCZ sfcc1" : "Katz sfcc0") << ": mean solves/day " << solves / (3.0 * 365)
            << ", max " << mx << ", " << us << " us/column-day (this build)");
  }
}
