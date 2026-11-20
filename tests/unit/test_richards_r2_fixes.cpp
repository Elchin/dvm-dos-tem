// Bug R2 offline candidates (TEM-dev 2026-11-20): richards_ksat_harmonic (A) and richards_lowsat_deriv_fix (B), tested on the two
// dumped real columns: C1 = S0_base sc 2041-05-23 (day 9267, legacy run DIAG_1, qinfil 0) and C2 = S0_base + sat_deriv_fix sc 2041-05-09
// (DIAG_6, qinfil 4.09768 mm/d, fibric bottom 88% saturated over saturated humic). Plus unsaturated columns (no-op / closure checks).
#include "doctest.h"
#include "test_helpers.h"
#include <cmath>
#include <json/value.h>
#include <string>

using ut::make_soil_layer;

namespace {
// z, dz, liq0, maxliq, minliq, poro, type(1 fib, 2 hum, 3 min)
const double C1[20][7] = {
 {0.01,0.02,3.7864,17.7421,0.887103,0.887103,1}, {0.03,0.04,13.6329,35.4841,1.77421,0.887103,1}, {0.07,0.0483044,42.851,42.851,2.14255,0.887103,1},
 {0.118304,0.0492935,43.7087,43.7087,2.18544,0.886703,2}, {0.167598,0.0985871,87.4174,87.4174,4.37087,0.886703,2},
 {0.266185,0.1,42.8482,42.8482,2.14241,0.428482,3}, {0.366185,0.1,42.8482,42.8482,2.14241,0.428482,3}, {0.466185,0.1,42.8482,42.8482,2.14241,0.428482,3},
 {0.566185,0.1,42.8482,42.8482,2.14241,0.428482,3}, {0.666185,0.2,85.6964,85.6964,4.28482,0.428482,3}, {0.866185,0.2,85.6964,85.6964,4.28482,0.428482,3},
 {1.06618,0.2,85.6964,85.6964,4.28482,0.428482,3}, {1.26618,0.3,128.545,128.545,6.42723,0.428482,3}, {1.56618,0.3,128.545,128.545,6.42723,0.428482,3},
 {1.86618,0.3,128.545,128.545,6.42723,0.428482,3}, {2.16618,0.3,128.545,128.545,6.42723,0.428482,3}, {2.46618,0.5,214.241,214.241,10.712,0.428482,3},
 {2.96618,0.5,214.241,214.241,10.712,0.428482,3}, {3.46618,1,428.482,428.482,21.4241,0.428482,3}, {4.46618,1,428.482,428.482,21.4241,0.428482,3}};
const double C2[20][7] = {
 {0.01,0.02,5.77584,17.7421,0.887103,0.887103,1}, {0.03,0.04,12.2394,35.4841,1.77421,0.887103,1}, {0.07,0.0489279,38.3604,43.4041,2.17021,0.887103,1},
 {0.118928,0.0492928,43.7081,43.7081,2.1854,0.886703,2}, {0.168221,0.0985856,87.4161,87.4161,4.37081,0.886703,2},
 {0.266806,0.1,42.8482,42.8482,2.14241,0.428482,3}, {0.366806,0.1,42.8482,42.8482,2.14241,0.428482,3}, {0.466806,0.1,42.8482,42.8482,2.14241,0.428482,3},
 {0.566806,0.1,42.8482,42.8482,2.14241,0.428482,3}, {0.666806,0.2,85.6964,85.6964,4.28482,0.428482,3}, {0.866806,0.2,85.6964,85.6964,4.28482,0.428482,3},
 {1.06681,0.2,85.6964,85.6964,4.28482,0.428482,3}, {1.26681,0.3,128.545,128.545,6.42723,0.428482,3}, {1.56681,0.3,128.545,128.545,6.42723,0.428482,3},
 {1.86681,0.3,128.545,128.545,6.42723,0.428482,3}, {2.16681,0.3,128.545,128.545,6.42723,0.428482,3}, {2.46681,0.5,214.241,214.241,10.712,0.428482,3},
 {2.96681,0.5,214.241,214.241,10.712,0.428482,3}, {3.46681,1,428.482,428.482,21.4241,0.428482,3}, {4.46681,1,428.482,428.482,21.4241,0.428482,3}};
const double HK[4] = {0, 191.321, 0.0103873, 0.00500493}, BS[4] = {0, 2.7, 8.0, 4.19032}, PS[4] = {0, -10.0, -12.0, -178.157};
std::vector<SoilLayer*> col(const double (*c)[7], double satscale = 1.0) {
  std::vector<SoilLayer*> v;
  for (int i = 0; i < 20; ++i) {
    const double* r = c[i]; const int t = (int)r[6];
    SoilLayer* l = make_soil_layer(r[0], r[1], r[5], r[2] * satscale, 0.0, 2.0, -1, 0.0);
    l->maxliq = r[3]; l->minliq = r[4]; l->hksat = HK[t]; l->bsw = BS[t]; l->psisat = PS[t];
    v.push_back(l);
  }
  return v;
}
double wsum(const ut::Column& c) { double w = 0; for (auto l : c.L) w += l->liq; return w; }
void day(Richards& rr, ut::Column& c, double infil_mm) {
  double tr[MAX_SOI_LAY] = {0}; Layer* dl = c.L.back(); const double dd = dl->z + dl->dz;
  rr.update(c.L[0], dl, dd, 0.0, dd, tr, 0.0, infil_mm / SEC_IN_DAY, 0.0);
}
struct Res { double err, err30, minfree, over; int nsub, nsub30, nfloor, nfloor30, overlay; };
Res run(const double (*c)[7], double infil, int sd, int kh, int ld) {
  ut::Column a(col(c)); Richards r; r.dt_fix = 1; r.sat_deriv_fix = sd; r.ksat_harmonic = kh; r.lowsat_deriv_fix = ld;
  Res o{}; day(r, a, infil); o.err = r.mass_err_mm; o.nsub = r.diag_nsub; o.nfloor = r.diag_nfloor;
  o.minfree = 1e30; o.over = -1e30; o.overlay = -1; for (int i = 0; i < 20; ++i) { Layer* l = a.L[i]; o.minfree = fmin(o.minfree, l->liq - l->minliq); if (l->liq - l->maxliq > o.over) { o.over = l->liq - l->maxliq; o.overlay = i; } }
  o.err30 = 0; o.nsub30 = 0; o.nfloor30 = 0;
  for (int d = 0; d < 30; ++d) { day(r, a, infil); o.err30 = fmax(o.err30, fabs(r.mass_err_mm)); o.nsub30 = std::max(o.nsub30, r.diag_nsub); o.nfloor30 += r.diag_nfloor;
    for (auto l : a.L) { o.minfree = fmin(o.minfree, l->liq - l->minliq); } }
  return o;
}
const char* NAMES[5] = {"legacy", "sat_deriv_fix", "A ksat_harmonic", "B lowsat_deriv_fix", "A+B"};
const int SW[5][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 1}};
} // namespace

TEST_SUITE("richards_r2_candidates") {

TEST_CASE("R2 candidates: settings default 0, read, invalid FATAL") {
  Richards rr; CHECK(rr.ksat_harmonic == 0); CHECK(rr.lowsat_deriv_fix == 0);
  ModelData m0; CHECK(m0.richards_ksat_harmonic == 0); CHECK(m0.richards_lowsat_deriv_fix == 0);
  Json::Value c1; c1["model_settings"]["richards_ksat_harmonic"] = 1; c1["model_settings"]["richards_lowsat_deriv_fix"] = 1;
  ModelData m1(c1); CHECK(m1.richards_ksat_harmonic == 1); CHECK(m1.richards_lowsat_deriv_fix == 1);
  CHECK(ut::exit_status_of([] { Json::Value c; c["model_settings"]["richards_ksat_harmonic"] = 2; ModelData m(c); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { Json::Value c; c["model_settings"]["richards_lowsat_deriv_fix"] = "1"; ModelData m(c); }) == EXIT_FAILURE);
}

TEST_CASE("R2 candidates on the two real columns: day closure, substeps, floor hits, 30 more days") {
  const double (*cols[2])[7] = {C1, C2}; const double inf[2] = {0.0, 4.09768}; const char* cn[2] = {"C1 2041-05-23 (legacy dump)", "C2 2041-05-09 (sat_deriv_fix dump)"};
  Res R[2][5];
  for (int c = 0; c < 2; ++c) for (int v = 0; v < 5; ++v) {
    R[c][v] = run(cols[c], inf[c], SW[v][0], SW[v][1], SW[v][2]); const Res& o = R[c][v];
    MESSAGE("R2TABLE | " << std::string(cn[c]) << " | " << std::string(NAMES[v]) << " | day err " << o.err << " mm | nsub " << o.nsub << " | floor hits " << o.nfloor
            << " | min(liq-minliq) " << o.minfree << " mm | max(liq-maxliq) " << o.over << " mm | 30d worst |err| " << o.err30 << " | 30d max nsub " << o.nsub30 << " | 30d floor hits " << o.nfloor30);
  }
  CHECK(fabs(R[0][0].err) > 100.0);   // C1: legacy creates water (DIAG_1 +1124 mm in the model, +870 mm here)
  // C2 from a 6-digit dump: legacy and sat_deriv_fix do not reproduce the model's +14798 mm, but the same instability shows as
  // thousands of in-solver floor hits (negative humic liquid clamped to 1e-4 mm) and a layer driven below minliq
  CHECK(R[1][0].nfloor > 1000); CHECK(R[1][1].nfloor > 1000); CHECK(R[1][0].minfree < 0.0);
  for (int c = 0; c < 2; ++c) for (int v : {2, 3, 4}) {
    CAPTURE(c); CAPTURE(v);
    CHECK(fabs(R[c][v].err) < 1e-6); CHECK(R[c][v].err30 < 1e-6); CHECK(R[c][v].nfloor == 0); CHECK(R[c][v].nfloor30 == 0);
    CHECK(R[c][v].nsub < 8640); CHECK(R[c][v].minfree > -1e-9);
  }
}

TEST_CASE("R2 candidates on unsaturated columns: B is a no-op (bit-identical), A closes and the size of its change is reported") {
  for (int c = 0; c < 2; ++c) {
    const double (*cc)[7] = c == 0 ? C1 : C2;
    ut::Column a(col(cc, 0.6)), b(col(cc, 0.6)), h(col(cc, 0.6));
    Richards ra; ra.dt_fix = 1; Richards rb; rb.dt_fix = 1; rb.lowsat_deriv_fix = 1; Richards rh; rh.dt_fix = 1; rh.ksat_harmonic = 1;
    double worst_h = 0; int nsa = 0, nsh = 0;
    for (int d = 0; d < 10; ++d) { day(ra, a, 2.0); day(rb, b, 2.0); day(rh, h, 2.0); worst_h = fmax(worst_h, fabs(rh.mass_err_mm)); nsa = std::max(nsa, ra.diag_nsub); nsh = std::max(nsh, rh.diag_nsub); }
    bool same = true; double dmax = 0; for (int i = 0; i < 20; ++i) { same = same && (a.L[i]->liq == b.L[i]->liq); dmax = fmax(dmax, fabs(a.L[i]->liq - h.L[i]->liq)); }
    MESSAGE("R2UNSAT | column " << (c + 1) << " at 60% of its liquid, 2 mm/d infil, 10 days | B bit-identical " << same << " | A worst |err| " << worst_h
            << " mm, max |layer liq A - legacy| " << dmax << " mm, max nsub legacy " << nsa << " A " << nsh);
    CHECK(same); CHECK(worst_h < 1e-6);
  }
  // uniform column: harmonic of equal ksat = ksat (A changes nothing beyond round-off)
  auto mk = [] { std::vector<SoilLayer*> v; for (int i = 0; i < 20; ++i) { SoilLayer* l = make_soil_layer(i * 0.05, 0.05, 0.4, 0.3 * 0.4 * DENLIQ * 0.05, 0.0, 5.0, -1, 0.0);
                 l->hksat = 0.005; l->psisat = -100.0; l->bsw = 5.0; v.push_back(l); } return v; };
  ut::Column a(mk()), h(mk()); Richards ra; ra.dt_fix = 1; Richards rh; rh.dt_fix = 1; rh.ksat_harmonic = 1;
  for (int d = 0; d < 3; ++d) { day(ra, a, 2.0); day(rh, h, 2.0); }
  double dm = 0; for (int i = 0; i < 20; ++i) dm = fmax(dm, fabs(a.L[i]->liq - h.L[i]->liq));
  MESSAGE("R2UNSAT | uniform column: max |A - legacy| " << dm << " mm");
  CHECK(dm < 1e-9);
}

}
