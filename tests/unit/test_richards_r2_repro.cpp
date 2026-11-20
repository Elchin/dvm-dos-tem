// Bug R2 regression (TEM-dev 2026-11-20): start-of-day soil column of CMT06 S0_base sc day 9267 (2041-05-23),
// dumped by TEM_RICH_DIAG from the DIAG_1 rerun (6 significant digits). Fibric (3 layers) over a deep layer split
// into 2 humic layers on 2041-02-01, saturated mineral below; no ice; qinfil = evap = 0.
// Legacy: the day does not close (DIAG_1: +1124 mm at 8640 substeps). richards_mass_guard 1 rejects it.
#include "doctest.h"
#include "test_helpers.h"
#include <lapacke.h>
#include <cmath>
#include <json/value.h>

using ut::make_soil_layer;

namespace {
// z, dz, liq0, maxliq, minliq, poro, hksat, bsw, psisat, type(1 fib, 2 hum, 3 min)
const double R2COL[20][10] = {
 {0.01,0.02,3.7864,17.7421,0.887103,0.887103,191.321,2.7,-10,1}, {0.03,0.04,13.6329,35.4841,1.77421,0.887103,191.321,2.7,-10,1},
 {0.07,0.0483044,42.851,42.851,2.14255,0.887103,191.321,2.7,-10,1}, {0.118304,0.0492935,43.7087,43.7087,2.18544,0.886703,0.0103873,8,-12,2},
 {0.167598,0.0985871,87.4174,87.4174,4.37087,0.886703,0.0103873,8,-12,2}, {0.266185,0.1,42.8482,42.8482,2.14241,0.428482,0.00500493,4.19032,-178.157,3},
 {0.366185,0.1,42.8482,42.8482,2.14241,0.428482,0.00500493,4.19032,-178.157,3}, {0.466185,0.1,42.8482,42.8482,2.14241,0.428482,0.00500493,4.19032,-178.157,3},
 {0.566185,0.1,42.8482,42.8482,2.14241,0.428482,0.00500493,4.19032,-178.157,3}, {0.666185,0.2,85.6964,85.6964,4.28482,0.428482,0.00500493,4.19032,-178.157,3},
 {0.866185,0.2,85.6964,85.6964,4.28482,0.428482,0.00500493,4.19032,-178.157,3}, {1.06618,0.2,85.6964,85.6964,4.28482,0.428482,0.00500493,4.19032,-178.157,3},
 {1.26618,0.3,128.545,128.545,6.42723,0.428482,0.00500493,4.19032,-178.157,3}, {1.56618,0.3,128.545,128.545,6.42723,0.428482,0.00500493,4.19032,-178.157,3},
 {1.86618,0.3,128.545,128.545,6.42723,0.428482,0.00500493,4.19032,-178.157,3}, {2.16618,0.3,128.545,128.545,6.42723,0.428482,0.00500493,4.19032,-178.157,3},
 {2.46618,0.5,214.241,214.241,10.712,0.428482,0.00500493,4.19032,-178.157,3}, {2.96618,0.5,214.241,214.241,10.712,0.428482,0.00500493,4.19032,-178.157,3},
 {3.46618,1,428.482,428.482,21.4241,0.428482,0.00500493,4.19032,-178.157,3}, {4.46618,1,428.482,428.482,21.4241,0.428482,0.00500493,4.19032,-178.157,3}};
std::vector<SoilLayer*> r2col() {
  std::vector<SoilLayer*> v;
  for (int i = 0; i < 20; ++i) {
    const double* r = R2COL[i];
    SoilLayer* l = make_soil_layer(r[0], r[1], r[5], r[2], 0.0, 2.0, -1, 0.0);
    l->maxliq = r[3]; l->minliq = r[4]; l->hksat = r[6]; l->bsw = r[7]; l->psisat = r[8];
    v.push_back(l);
  }
  return v;
}
double r2sum(const ut::Column& c) { double w = 0; for (auto l : c.L) w += l->liq; return w; }
void r2day(Richards& rr, ut::Column& c) {
  double tr[MAX_SOI_LAY] = {0}; Layer* dl = c.L.back(); const double dd = dl->z + dl->dz;
  rr.update(c.L[0], dl, dd, 0.0, dd, tr, 0.0, 0.0, 0.0);
}
} // namespace

TEST_SUITE("richards_r2") {

TEST_CASE("bug R2: one 10 s substep of the real column - linear system, conservation, coefficients") {
  ut::Column c(r2col()); Richards rr; rr.dt_fix = 1;
  Layer* top = c.L[0]; Layer* dl = c.L.back(); const int ti = top->indl, di = dl->indl;
  rr.clearRichardsArrays(); rr.qinfil = 0.0; rr.qevap = 0.0;
  rr.prepareSoilColumn(top, di);
  const double dts = 10.0;
  for (int i = ti; i <= di + 1; i++) { rr.dt_dz[i] = rr.dzmm[i] > 0 ? dts / rr.dzmm[i] : 0.0;
    rr.vol_liq[i] = fmax(rr.soi_liq[i], 1.e-6) / (DENLIQ * (rr.dzmm[i] / 1.e3)); rr.vol_h2o[i] = fmax(fmin(rr.vol_liq[i] + rr.vol_ice[i], 1.0), 0.0); }
  rr.computeHydraulicProperties(top, di); rr.computeMoistureFluxesAndDerivs(top, ti, di);
  rr.computeLHS(top, ti, di); rr.computeRHS(top, ti, di);
  const int n = rr.num_al; double sub[64], dg[64], sup[64], res[64], A[64], B[64], C[64], Rr[64];
  for (int k = 0; k < n; k++) { dg[k] = rr.bmx[k + ti]; res[k] = rr.rmx[k + ti]; B[k] = dg[k]; Rr[k] = res[k]; }
  for (int k = 0; k < n - 1; k++) { sub[k] = rr.amx[k + ti + 1]; sup[k] = rr.cmx[k + ti]; A[k] = sub[k]; C[k] = sup[k]; }
  int info = LAPACKE_dgtsv(LAPACK_ROW_MAJOR, n, 1, sub, dg, sup, res, 1);
  double sumc = 0, maxres = 0;
  for (int k = 0; k < n; k++) {
    sumc += rr.dzmm[k + ti] * res[k];
    double r = B[k] * res[k] - Rr[k] + (k > 0 ? A[k - 1] * res[k - 1] : 0) + (k < n - 1 ? C[k] * res[k + 1] : 0);
    maxres = fmax(maxres, fabs(r));
  }
  for (int k = 0; k < n; k++) {
    int i = k + ti;
    MESSAGE("row " << i << " dz " << rr.dzmm[i] << " vol " << rr.vol_liq[i] << " lporo " << rr.liq_poro[i] << " s2 " << rr.s2[i] << " hk " << rr.hk[i] << " smp " << rr.smp[i]
            << " qin " << rr.qin[i] << " qout " << rr.qout[i] << " dqidw0 " << rr.dqidw0[i] << " dqidw1 " << rr.dqidw1[i] << " dqodw1 " << rr.dqodw1[i] << " dqodw2 " << rr.dqodw2[i]
            << " | a " << rr.amx[i] << " b " << rr.bmx[i] << " c " << rr.cmx[i] << " r " << rr.rmx[i] << " dtheta " << res[k]);
  }
  // telescoping check on the coefficients themselves: dqidw0[i] == dqodw1[i-1], dqidw1[i] == dqodw2[i-1]
  double tele = 0; for (int i = ti + 1; i <= di; i++) tele = fmax(tele, fabs(rr.dqidw0[i] - rr.dqodw1[i - 1]) + fabs(rr.dqidw1[i] - rr.dqodw2[i - 1]));
  MESSAGE("info " << info << " sum(dz*dtheta) " << sumc << " mm (expected 0), max row residual " << maxres << ", coefficient telescoping mismatch " << tele);
  CHECK(info == 0);
}

TEST_CASE("bug R2: the real column does not close in legacy; richards_mass_guard 1 rejects the day") {
  ut::Column a(r2col()), b(r2col());
  Richards ra; ra.dt_fix = 1; const double w0 = r2sum(a); r2day(ra, a);
  MESSAGE("legacy (dt_fix 1, no guard): mass_err " << ra.mass_err_mm << " mm, column change " << r2sum(a) - w0 << " mm, substeps " << ra.diag_nsub);
  CHECK(fabs(ra.mass_err_mm) > 1.0);   // documents bug R2 on the real state
  Richards rb; rb.dt_fix = 1; rb.mass_guard = 1; const double wb = r2sum(b); r2day(rb, b);
  CHECK(rb.guard_tripped); CHECK(rb.n_guard == 1); CHECK(r2sum(b) == wb);
}

}

TEST_SUITE("richards_r2_fix") {

TEST_CASE("richards_sat_deriv_fix: default 0, read, invalid FATAL") {
  Richards rr; CHECK(rr.sat_deriv_fix == 0);
  ModelData m0; CHECK(m0.richards_sat_deriv_fix == 0);
  Json::Value c1; c1["model_settings"]["richards_sat_deriv_fix"] = 1; ModelData m1(c1); CHECK(m1.richards_sat_deriv_fix == 1);
  CHECK(ut::exit_status_of([] { Json::Value c; c["model_settings"]["richards_sat_deriv_fix"] = 2; ModelData m(c); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { Json::Value c; c["model_settings"]["richards_sat_deriv_fix"] = true; ModelData m(c); }) == EXIT_FAILURE);
}

TEST_CASE("bug R2 root cause: with richards_sat_deriv_fix 1 the real column closes, stays physical and needs few substeps") {
  ut::Column a(r2col()); Richards ra; ra.dt_fix = 1; ra.sat_deriv_fix = 1;
  const double w0 = r2sum(a); r2day(ra, a);
  double lmin = 1e30, over = 0; for (auto l : a.L) { lmin = fmin(lmin, l->liq - l->minliq); over = fmax(over, l->liq - l->maxliq); }
  MESSAGE("sat_deriv_fix 1: mass_err " << ra.mass_err_mm << " mm, column change " << r2sum(a) - w0 << " mm, substeps " << ra.diag_nsub
          << ", min(liq - minliq) " << lmin << " mm, max(liq - maxliq) " << over << " mm");
  CHECK(fabs(ra.mass_err_mm) < 1e-6);
  CHECK(ra.diag_nsub < 8640);
  CHECK(lmin > -1e-6);
  // 30 more days on the same column with the fix: every day closes
  double worst = 0; for (int d = 0; d < 30; ++d) { r2day(ra, a); worst = fmax(worst, fabs(ra.mass_err_mm)); }
  MESSAGE("30 more days: worst |mass_err| " << worst << " mm");
  CHECK(worst < 1e-6);
}

TEST_CASE("richards_sat_deriv_fix 1 is bit-identical on an unsaturated column (no clamp active)") {
  auto mk = [] { std::vector<SoilLayer*> v; for (int i = 0; i < 20; ++i) { SoilLayer* l = make_soil_layer(i * 0.05, 0.05, 0.4, 0.3 * 0.4 * DENLIQ * 0.05, 0.0, 5.0, -1, 0.0);
                 l->hksat = 0.005; l->psisat = -100.0; l->bsw = 5.0; v.push_back(l); } return v; };
  ut::Column a(mk()), b(mk()); Richards ra; ra.dt_fix = 1; Richards rb; rb.dt_fix = 1; rb.sat_deriv_fix = 1;
  double tr[MAX_SOI_LAY] = {0};
  for (int d = 0; d < 3; ++d) {
    const double dd = a.L.back()->z + a.L.back()->dz;
    ra.update(a.L[0], a.L.back(), dd, 0.0, dd, tr, 0.0, 2.0 / SEC_IN_DAY, 0.0);
    rb.update(b.L[0], b.L.back(), dd, 0.0, dd, tr, 0.0, 2.0 / SEC_IN_DAY, 0.0);
  }
  for (int i = 0; i < 20; ++i) CHECK(a.L[i]->liq == b.L[i]->liq);
}

}
