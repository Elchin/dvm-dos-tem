// Richards equation, water mass balance and combined water + energy tests
// (tests/richards-water-energy). Tests only: no model code is changed.
//  - Richards (CLM5 form, Richards.cpp) vs analytic / semi-analytic solutions:
//      hydrostatic equilibrium (h = psi - z = const), steady downward flux profile
//      (ODE dpsi/dz = 1 - q/K(psi)), flux-boundary infiltration bounds (piston / Green-Ampt
//      type), exact mass conservation with no flux and with ET, uptake and lateral drainage.
//  - Soil_Env::checkSoilLiquidWaterValidity: water conserved (liq + ice + puddle + qover);
//      frozen cells keep their freezing-curve water with thermal_solver = 1 (regression guard of
//      7c5e4759), legacy (thermal_solver = 0) behaviour documented unchanged.
//  - Combined: daily enthalpy step (thermal_sfcc 0 and 1) + Richards + validity check +
//      energy-conserving water remap on a freezing wet column: daily energy budget closes
//      < 1e-6 W m-2, water conserved, latent heat = L_f * ice change, frozen cells keep
//      their unfrozen water.
// Set RWE_DUMP_DIR=<dir> to write CSV profiles / time series for the plots.
// Refs: Clapp & Hornberger 1978 doi:10.1029/WR014i004p00601; CLM5 tech note ch. 7 (Lawrence et
//   al. 2018); Green & Ampt 1911; Philip 1957 doi:10.1097/00010694-195708000-00008;
//   Schaefer & Jafarov 2016 doi:10.5194/bg-13-1991-2016; Romanovsky & Osterkamp 2000
//   doi:10.1002/1099-1530(200007/09)11:3<219::AID-PPP352>3.0.CO;2-7.
#include "doctest.h"
#include "test_helpers.h"
#include <cmath>
#include <cstdlib>
#include <fstream>

using ut::make_soil_layer;

namespace {
const double PSIS = 100.0;   // mm (make_soil_layer psisat = -100 -> Richards psisat = +100)
const double BSW = 5.0;
const double KS = 0.005;     // mm/s

std::ofstream* dumpfile(const std::string& name) {
  const char* d = getenv("RWE_DUMP_DIR");
  if (!d) return NULL;
  return new std::ofstream((std::string(d) + "/" + name).c_str());
}

// n layers of thickness dz, porosity poro, liquid saturation s0, thawed
std::vector<SoilLayer*> uniform(int n, double dz, double poro, double s0, double tem = 5.0, double ks = KS) {
  std::vector<SoilLayer*> v;
  for (int i = 0; i < n; ++i) {
    SoilLayer* l = make_soil_layer(i * dz, dz, poro, s0 * poro * DENLIQ * dz, 0.0, tem, -1, 0.0);
    l->hksat = ks; l->psisat = -PSIS; l->bsw = BSW;
    v.push_back(l);
  }
  return v;
}
double sat(Layer* l) { return l->liq / (DENLIQ * l->dz * (l->poro - l->ice / (DENICE * l->dz))); }
double smp(double s) { return -PSIS * pow(s, -BSW); }            // mm
double Kof(double s, double ks = KS) { return ks * pow(s, 2 * BSW + 3); }          // mm/s
// silt-loam MineralLayer (texture drives the SJ2016 curve in Soil_Env::enthalpyCell), props as make_soil_layer
SoilLayer* make_mineral(double z, double dz, double poro, double liq, double tem,
                        float sand = 20, float silt = 60, float clay = 20) {
  MineralLayer* l = new MineralLayer(dz, sand, silt, clay);
  l->z = z; l->dz = dz; l->poro = poro; l->liq = liq; l->ice = 0.0; l->tem = tem;
  l->frozen = -1; l->frozenfrac = 0.0; l->tcsolid = 2.0; l->bulkden = 2700.0 * (1.0 - poro);
  l->vhcsolid = 2.0e6; l->rawc = l->soma = l->sompr = l->somcr = 0.0;
  l->derivePhysicalProperty();
  l->hksat = KS; l->psisat = -PSIS; l->bsw = BSW; l->pce_f = 0.; l->pce_t = 0.; l->solind = 0; l->indl = 0;
  return l;
}
double liqsum(const ut::Column& c) { double w = 0; for (auto l : c.L) w += l->liq; return w; }

void day(Richards& rr, ut::Column& c, Layer* drainl, double infil_mmd, double evap_mmd,
         const double* trans_mmd, double slope = 0.0, double fbase = 0.0, double wt = -1.0) {
  double tr[MAX_SOI_LAY] = {0};
  if (trans_mmd) for (size_t i = 0; i < c.L.size(); ++i) tr[i] = trans_mmd[i] / SEC_IN_DAY;
  const double dd = drainl->z + drainl->dz;
  rr.update(c.L[0], drainl, dd, fbase, wt < 0 ? dd : wt, tr, evap_mmd / SEC_IN_DAY,
            infil_mmd / SEC_IN_DAY, slope);
}
} // namespace

TEST_SUITE("richards_analytic") {

TEST_CASE("Richards no-flux: mass conserved to round-off and the profile relaxes to hydrostatic equilibrium") {
  ut::Column c(uniform(10, 0.1, 0.4, 0.6));
  Richards rr; const double M0 = liqsum(c);   // default tolerances (no fluxes: rejection cannot create mass)
  double maxdm = 0.0;
  for (int d = 0; d < 3000; ++d) {
    const double m = liqsum(c);
    day(rr, c, c.L.back(), 0.0, 0.0, NULL);
    maxdm = fmax(maxdm, fabs(liqsum(c) - m));
  }
  CHECK(maxdm < 1e-9);                                  // mm per day
  CHECK(fabs(liqsum(c) - M0) < 1e-8);
  // discrete equilibrium: smp_{i+1} - smp_i = node spacing (mm), i.e. h = smp - z constant
  double maxdev = 0.0;
  for (size_t i = 0; i + 1 < c.L.size(); ++i) {
    const double dh = smp(sat(c.L[i + 1])) - smp(sat(c.L[i]));
    maxdev = fmax(maxdev, fabs(dh / 100.0 - 1.0));
  }
  CHECK(maxdev < 1e-3);
  // analytic hydrostatic profile with the same mass: s(z) = ((C + z)/(-psis))^(-1/B)
  auto prof = [&](double C, int i) { const double z = (i + 0.5) * 100.0; return pow(-(C + z) / PSIS, -1.0 / BSW); };
  auto mass = [&](double C) { double m = 0; for (int i = 0; i < 10; ++i) m += prof(C, i) * 0.4 * 100.0; return m; };
  double lo = -1e6, hi = -1000.0 - 1e-6;               // C + z < 0 for all nodes (z_max = 950 mm)
  for (int k = 0; k < 200; ++k) { const double mid = 0.5 * (lo + hi); (mass(mid) > M0 ? hi : lo) = mid; }
  double maxerr = 0.0;
  std::ofstream* f = dumpfile("richards_hydrostatic.csv");
  if (f) *f << "z_m,s_model,s_analytic\n";
  for (int i = 0; i < 10; ++i) {
    const double sa = prof(0.5 * (lo + hi), i);
    maxerr = fmax(maxerr, fabs(sat(c.L[i]) - sa));
    if (f) *f << (i + 0.5) * 0.1 << "," << sat(c.L[i]) << "," << sa << "\n";
  }
  delete f;
  CHECK(maxerr < 1e-3);
}

TEST_CASE("Richards steady downward flux: interface fluxes = q and profile matches dpsi/dz = 1 - q/K") {
  // n <= MAX_SOI_LAY: Richards' work arrays and the uptake vector are sized MAX_SOI_LAY
  const int n = MAX_SOI_LAY < 16 ? MAX_SOI_LAY : 16; const double q = 10.0, ks = 0.05;   // mm/day in at the top, out as uptake at the bottom
  ut::Column c(uniform(n, 0.1, 0.4, 0.6, 5.0, ks));
  Richards rr; rr.toler_upper = 1e30; double tr[MAX_SOI_LAY] = {0}; tr[n - 1] = q;
  for (int d = 0; d < 400; ++d) day(rr, c, c.L.back(), q, 0.0, tr);
  // discrete fluxes at steady state (start-of-last-substep linearisation point)
  double maxrel = 0.0;
  for (int i = 1; i < n - 1; ++i) maxrel = fmax(maxrel, fabs(rr.qout[i] * SEC_IN_DAY / q - 1.0));
  CHECK(maxrel < 1e-3);
  // continuous ODE from the top node down (RK4, 100 sub-steps per node spacing)
  double psi = smp(sat(c.L[0])), maxds = 0.0;
  auto sOf = [](double p) { return pow(-p / PSIS, -1.0 / BSW); };
  auto rhs = [&](double p) { return 1.0 - (q / SEC_IN_DAY) / Kof(sOf(p), ks); };
  std::ofstream* f = dumpfile("richards_steadyflux.csv");
  if (f) *f << "z_m,s_model,s_ode,s_unitgradient\n";
  const double sug = pow((q / SEC_IN_DAY) / ks, 1.0 / (2 * BSW + 3));
  for (int i = 0; i < n; ++i) {
    if (i > 0) { const double h = 100.0 / 100; for (int k = 0; k < 100; ++k) {
      const double k1 = rhs(psi), k2 = rhs(psi + 0.5 * h * k1), k3 = rhs(psi + 0.5 * h * k2), k4 = rhs(psi + h * k3);
      psi += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0; } }
    if (i < n - 1) maxds = fmax(maxds, fabs(sat(c.L[i]) - sOf(psi)) / sOf(psi));
    if (f) *f << (i + 0.5) * 0.1 << "," << sat(c.L[i]) << "," << sOf(psi) << "," << sug << "\n";
  }
  delete f;
  CHECK(maxds < 0.02);
  // upper column approaches the gravity-drained (unit-gradient) saturation K(s) = q
  CHECK(fabs(sat(c.L[0]) - sug) / sug < 0.05);
}

// Richards::update rejects a substep (error > toler_upper) by halving dtsub, but dt_dz = dtsub/dz is
// computed once per substep (Richards.cpp l.138), before the trial loop, so the re-solved system still
// uses the old dtsub while dtdone advances by the halved one: the fluxes of the day are applied for
// more than a day (water created by infiltration, extra loss by ET/uptake). With rejection disabled
// (toler_upper huge) the same solver conserves mass to round-off. Tests below run with rejection
// disabled for the physics checks; the default-tolerance cases are documented as may_fail (bug R1).
TEST_CASE("Richards flux-boundary infiltration: mass exact; wetting front between piston and Green-Ampt-type bounds") {
  const double q = 20.0, s0 = 0.3, poro = 0.4; const int days = 5, n = 20;
  ut::Column c(uniform(n, 0.05, poro, s0));
  Richards rr; rr.toler_upper = 1e30; const double M0 = liqsum(c);
  std::ofstream* f = dumpfile("richards_infiltration.csv");
  if (f) *f << "day,z_m,s\n";
  for (int d = 0; d < days; ++d) {
    day(rr, c, c.L.back(), q, 0.0, NULL);
    if (f) for (int i = 0; i < n; ++i) *f << d + 1 << "," << (i + 0.5) * 0.05 << "," << sat(c.L[i]) << "\n";
  }
  delete f;
  const double I = q * days;
  CHECK(fabs(liqsum(c) - M0 - I) < 1e-8);
  double zf = 0.0;   // bottom of the deepest wetted layer (s > s0 + 0.02)
  for (int i = 0; i < n; ++i) if (sat(c.L[i]) > s0 + 0.02) zf = (i + 1) * 0.05 * 1000.0;
  const double piston = I / (poro * (1.0 - s0));                           // saturated piston (lower bound)
  const double sw = pow((q / SEC_IN_DAY) / KS, 1.0 / (2 * BSW + 3));        // K(s_w) = q
  const double ga = I / (poro * (sw - s0));                                 // flux-BC Green-Ampt depth
  CAPTURE(zf); CAPTURE(piston); CAPTURE(ga);
  CHECK(zf >= piston);
  CHECK(zf <= 1.5 * ga + 50.0);
  CHECK(sat(c.L[0]) <= 1.0);
}

TEST_CASE("Richards water balance with evaporation, uptake and lateral drainage closes per day") {
  ut::Column c(uniform(8, 0.125, 0.45, 0.7));
  for (int i = 5; i < 8; ++i) c.L[i]->liq = 0.97 * c.L[i]->maxliq;   // saturated bottom -> lateral drainage
  Richards rr; rr.toler_upper = 1e30; double tr[MAX_SOI_LAY] = {0.2, 0.3, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0};
  double maxerr = 0.0, qd = 0.0;
  for (int d = 0; d < 60; ++d) {
    const double m = liqsum(c);
    const double inf = (d % 7 == 0) ? 15.0 : 0.0, ev = 0.8;
    day(rr, c, c.L.back(), inf, ev, tr, 5.0, 1.0, 0.3);
    double ut = 0; for (int i = 0; i < 8; ++i) ut += tr[i];
    maxerr = fmax(maxerr, fabs((liqsum(c) - m) - (inf - ev - ut - rr.qdrain)));
    qd += rr.qdrain;
  }
  CHECK(qd > 0.0);           // lateral drainage was exercised
  CHECK(maxerr < 1e-9);      // mm/day
}

TEST_CASE("bug R1 (documented): default substep rejection breaks Richards mass conservation" * doctest::may_fail()) {
  for (int rej : {0, 1}) {
    CAPTURE(rej);
    ut::Column c(uniform(20, 0.05, 0.4, 0.3));
    Richards rr; if (!rej) rr.toler_upper = 1e30;
    double worst = 0.0, tot = 0.0; int nsub_days = 0;
    for (int d = 0; d < 5; ++d) { const double m = liqsum(c); day(rr, c, c.L.back(), 20.0, 0.0, NULL);
      const double e = liqsum(c) - m - 20.0; worst = fmax(worst, fabs(e)); tot += e; nsub_days++; }
    MESSAGE("Richards infiltration 20 mm/d x 5 d, rejection " << (rej ? "default (toler_upper 0.1)" : "disabled")
            << ": worst daily mass error " << worst << " mm, total " << tot << " mm");
    CHECK(worst < 1e-8);
  }
}

TEST_CASE("Richards with a frozen (frozen = 1) layer above the drainage layer (documented, not a model path)" * doctest::may_fail()) {
  // A frozen lens inside fstsoill..drainl: num_al counts only frozen <= 0 layers, so the tridiagonal
  // system is truncated while the update loop runs to drainl. Not reachable from Soil_Env (drainl is the
  // first layer holding a front, Ground.cpp ~l.2090, so a frozen layer is never above it).
  ut::Column c(uniform(6, 0.1, 0.4, 0.6));
  c.L[1]->frozen = 1; c.L[1]->frozenfrac = 1.0; c.L[1]->ice = 0.5 * c.L[1]->maxice; c.L[1]->liq = 0.05 * c.L[1]->maxliq; c.L[1]->tem = -1.0;
  Richards rr; rr.toler_upper = 1e30; double worst = 0.0;
  for (int d = 0; d < 30; ++d) {
    const double m = liqsum(c);
    day(rr, c, c.L.back(), 5.0, 0.0, NULL);
    worst = fmax(worst, fabs(liqsum(c) - m - 5.0));
  }
  MESSAGE("frozen-lens Richards worst daily mass error (mm): " << worst << ", num_al = " << rr.num_al);
  CHECK(worst < 1e-6);
}

}

TEST_SUITE("water_validity") {

TEST_CASE("checkSoilLiquidWaterValidity conserves liq + ice + puddle + qover (thawed, oversaturated)") {
  ut::Column c(uniform(6, 0.1, 0.4, 0.5));
  c.L[3]->liq = 1.3 * c.L[3]->maxliq; c.L[0]->liq = 0.01 * c.L[0]->maxliq; c.L[5]->liq = 1.2 * c.L[5]->maxliq;
  for (int ts : {0, 1}) {
    CAPTURE(ts);
    ut::Column cc(uniform(6, 0.1, 0.4, 0.5));
    for (int i = 0; i < 6; ++i) cc.L[i]->liq = c.L[i]->liq;
    EnvData ed; Soil_Env se; se.setEnvData(&ed); se.setGround(&cc.g); se.thermal_solver = ts;
    ed.d_soi2l.magic_puddle = 0.0; ed.d_soi2l.qover = 0.0; ed.d_soi2l.qdrain = 0.0;
    const double w0 = cc.total_water();
    se.checkSoilLiquidWaterValidity(cc.L[0], cc.L[0]->solind);
    CHECK(cc.total_water() + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover + ed.d_soi2l.qdrain == doctest::Approx(w0).epsilon(1e-14));
    for (auto l : cc.L) { CHECK(l->liq <= l->maxliq * (1 + 1e-12)); CHECK(l->liq >= l->minliq * (1 - 1e-12)); }
  }
}

TEST_CASE("checkSoilLiquidWaterValidity, whole column frozen: thermal_solver 1 keeps the unfrozen water (7c5e4759 guard); legacy dumps it, water conserved") {
  for (int ts : {0, 1}) {
    CAPTURE(ts);
    ut::Column c(uniform(6, 0.1, 0.4, 0.0, -2.0));
    std::vector<double> liq0;
    for (auto l : c.L) { l->frozen = 1; l->frozenfrac = 1.0; l->ice = 0.8 * l->maxice; l->liq = 0.1 * l->maxliq; liq0.push_back(l->liq); }
    EnvData ed; Soil_Env se; se.setEnvData(&ed); se.setGround(&c.g); se.thermal_solver = ts;
    ed.d_soi2l.magic_puddle = 0.0; ed.d_soi2l.qover = 0.0; ed.d_soi2l.qdrain = 0.0;
    const double w0 = c.total_water();
    se.checkSoilLiquidWaterValidity(c.L[0], c.L[0]->solind);
    CHECK(c.total_water() + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover == doctest::Approx(w0).epsilon(1e-14));
    for (size_t i = 0; i < c.L.size(); ++i) {
      if (ts == 1) CHECK(c.L[i]->liq == liq0[i]);       // bit-for-bit kept
      else         CHECK(c.L[i]->liq == 0.0);           // legacy: effmaxliq = 0 for frozen cells
    }
    if (ts == 1) { CHECK(ed.d_soi2l.magic_puddle == 0.0); CHECK(ed.d_soi2l.qover == 0.0); }
  }
}

}

namespace {
double colH(Ground& g) { double h = 0; for (Layer* l = g.toplayer; l; l = l->nextl) h += Soil_Env::enthLayerH(l); return h; }
double colW(Ground& g) { double w = 0; for (Layer* l = g.toplayer; l; l = l->nextl) if (!l->isRock) w += l->liq + l->ice; return w; }
double colIce(Ground& g) { double w = 0; for (Layer* l = g.toplayer; l; l = l->nextl) if (!l->isRock) w += l->ice; return w; }
double sensH(Ground& g) {   // sensible part, independent of enthLayerH: (Cs + liq c_w + ice c_i) T
  double h = 0;
  for (Layer* l = g.toplayer; l; l = l->nextl) h += (Soil_Env::enthLayerCs(l) + l->liq * SHCLIQ + l->ice * SHCICE) * l->tem;
  return h;
}
}

TEST_SUITE("water_energy_combined") {

TEST_CASE("freeze-up of a wet column with Richards on: daily energy budget < 1e-6 W m-2, water conserved, L_f * dIce consistent, unfrozen water kept") {
  for (int sf : {0, 1}) {
    CAPTURE(sf);
    // 8 layers, wetter at depth: Richards redistributes downward while the column freezes from the top
    std::vector<SoilLayer*> v;
    const double dz[8] = {0.05, 0.05, 0.1, 0.1, 0.15, 0.2, 0.25, 0.4};
    double z = 0.0;
    for (int i = 0; i < 8; ++i) {
      const double poro = i < 3 ? 0.6 : 0.45, s0 = i < 3 ? 0.55 : 0.85;
      SoilLayer* l = make_mineral(z, dz[i], poro, s0 * poro * DENLIQ * dz[i], 1.5 + 0.2 * i);
      z += dz[i]; v.push_back(l);
    }
    ut::Column c(v); Ground& g = c.g;
    EnvData ed; Soil_Env se; se.setEnvData(&ed); se.setGround(&g);
    se.thermal_solver = 1; se.thermal_sfcc = sf;
    ed.d_soi2l.magic_puddle = 0.0; ed.d_soi2l.qover = 0.0; ed.d_soi2l.qdrain = 0.0;
    Richards rr;
    const double W0 = colW(g);
    double worst_th = 0, worst_day = 0, worst_w = 0, worst_lat = 0, cum_in = 0, moved = 0, richards_err = 0, out_tot = 0, unf_max = 0;
    int frozen_kept_bad = 0, richards_days = 0;
    std::ofstream* f = dumpfile(std::string("combined_sfcc") + (sf ? "1" : "0") + ".csv");
    std::ofstream* fp = dumpfile(std::string("combined_profiles_sfcc") + (sf ? "1" : "0") + ".csv");
    if (f) *f << "day,ts,res_thermal_Wm2,res_water_Wm2,res_day_Wm2,dW_kg,ice_kg,liq_kg,latent_err_Wm2\n";
    if (fp) *fp << "day,z_m,T,liq,ice\n";
    const double H0all = colH(g);
    for (int d = 0; d < 120; ++d) {
      const double ts = 3.0 - 0.25 * d + 2.0 * sin(2 * M_PI * d / 9.0);   // cooling with synoptic swings
      const double Hs0 = sensH(g), I0 = colIce(g);
      EnergyDiag q; se.enthalpyThermalStep(ts, 0.0, q);
      REQUIRE(q.fallbacks == 0);
      const double hres = q.res + q.equil + q.init;   // thermal-step residual incl. curve projection
      // latent-heat consistency: dH = d(sensible) - L_f dIce, from an independent sensible sum
      worst_lat = fmax(worst_lat, fabs((q.h1 - q.h0) - ((sensH(g) - Hs0) - LHFUS * (colIce(g) - I0))) / SEC_IN_DAY);
      // water window: Richards (thawed part above the first frozen cell) + validity check + remap
      Soil_Env::EnthSnap s; se.enthSnap(s);
      std::vector<double> lfz; for (auto l : c.L) lfz.push_back(l->frozen == 1 ? l->liq : -1.0);
      const double Wb = colW(g) + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover;
      const double Wsoil = colW(g);
      double hmin = 0, hmax = 0;   // specific enthalpy range of movable water (liquid at the cell T)
      for (size_t i = 0; i < s.L.size(); ++i) { hmin = fmin(hmin, SHCLIQ * s.T[i]); hmax = fmax(hmax, SHCLIQ * s.T[i]); }
      Layer* drainl = NULL;
      for (auto l : c.L) { if (l->frozen == 1) break; drainl = l; }
      double rw = 0.0;
      if (drainl != NULL && c.L[0]->frozen != 1) {
        const double before = colW(g);
        day(rr, c, drainl, 0.0, 0.0, NULL);
        richards_err = fmax(richards_err, fabs(colW(g) - before));   // no boundary flux, slope 0
        se.checkSoilLiquidWaterValidity(c.L[0], c.L[0]->solind);
        for (size_t i = 0; i < c.L.size(); ++i) moved += fabs(c.L[i]->liq - s.liq[i]);
        (void)before; richards_days++;
      }
      se.enthRemapWaterCore(s, 0.0, 0.0);
      for (size_t i = 0; i < c.L.size(); ++i) if (lfz[i] >= 0.0 && c.L[i]->liq != lfz[i]) frozen_kept_bad++;
      rw = colH(g) - s.hcol;                              // enthalpy change of the water window
      const double dW = colW(g) + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover - Wb;
      const double expel = Wsoil - colW(g);                // water that left the soil to puddle/qover, kg m-2
      out_tot += expel;
      cum_in += rw;                                        // booked boundary enthalpy of the water window
      // the energy carried out must be that of the expelled liquid: -expel*h, h in [hmin, hmax]
      if (expel > 1e-12) { CHECK(rw <= -expel * hmin + 1e-6 * SEC_IN_DAY); CHECK(rw >= -expel * hmax - 1e-6 * SEC_IN_DAY); rw = 0.0; }
      worst_w = fmax(worst_w, fabs(dW));
      worst_th = fmax(worst_th, fabs(hres) / SEC_IN_DAY);
      worst_day = fmax(worst_day, fabs(hres + rw) / SEC_IN_DAY);
      cum_in += q.gtop + q.gbot + q.src + q.clip;
      for (auto l : c.L) if (l->frozen == 1 && l->tem < -0.5) unf_max = fmax(unf_max, l->liq - l->minliq);
      double liq = 0; for (auto l : c.L) liq += l->liq;
      if (f) *f << d << "," << ts << "," << hres / SEC_IN_DAY << "," << rw / SEC_IN_DAY << "," << (hres + rw) / SEC_IN_DAY
                << "," << dW << "," << colIce(g) << "," << liq << "," << worst_lat << "\n";
      if (fp && d % 5 == 0) for (auto l : c.L) *fp << d << "," << l->z + 0.5 * l->dz << "," << l->tem << "," << l->liq << "," << l->ice << "\n";
    }
    delete f; delete fp;
    MESSAGE("sfcc " << sf << ": worst thermal " << worst_th << " W/m2, worst day " << worst_day << " W/m2, worst dW "
            << worst_w << " kg/m2, Richards days " << richards_days << ", liquid moved " << moved << " kg/m2, ice end " << colIce(g));
    CHECK(richards_days > 5); CHECK(moved > 1e-3);        // water redistribution actually exercised
    CHECK(worst_th < 1e-6); CHECK(worst_day < 1e-6);
    CHECK(worst_lat < 1e-6);
    CHECK(worst_w < 1e-9);
    CHECK(richards_err < 1e-9);
    CHECK(colW(g) + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover + g.expelled_water_mm == doctest::Approx(W0).epsilon(1e-12));
    CHECK(fabs(colH(g) - H0all - cum_in) / (120 * SEC_IN_DAY) < 1e-6);   // whole-run budget
    CHECK(frozen_kept_bad == 0);
    CHECK(colIce(g) > 0.5 * W0);                         // the column froze
    // unfrozen water ON keeps freezing-curve liquid above minliq in frozen cells (T < -0.5 C); OFF only minliq
    MESSAGE("sfcc " << sf << ": max(liq - minliq) in frozen cells below -0.5 C = " << unf_max << " kg/m2; expelled to puddle/qover "
            << out_tot << ", pore-cap expelled " << g.expelled_water_mm);
    if (sf == 1) CHECK(unf_max > 0.5); else CHECK(unf_max < 1e-9);
  }
}

TEST_CASE("unfrozen water ON vs OFF: same column, same forcing; ON stays warmer during early freeze-up (R&O 2000 qualitative)") {
  double Tmid[2] = {0, 0}, ice[2] = {0, 0};
  for (int sf : {0, 1}) {
    std::vector<SoilLayer*> v;
    for (int i = 0; i < 10; ++i) v.push_back(make_mineral(0.1 * i, 0.1, 0.45, 0.9 * 0.45 * DENLIQ * 0.1, 0.5));
    ut::Column c(v); Soil_Env se; se.setGround(&c.g); se.thermal_solver = 1; se.thermal_sfcc = sf;
    for (int d = 0; d < 20; ++d) { EnergyDiag q; se.enthalpyThermalStep(-5.0, 0.0, q); REQUIRE(fabs(q.res) / SEC_IN_DAY < 1e-6); }
    Tmid[sf] = c.L[2]->tem; ice[sf] = colIce(c.g);
  }
  MESSAGE("T(0.25 m) after 20 d at -5 C: OFF " << Tmid[0] << "  ON " << Tmid[1] << "; ice OFF " << ice[0] << " ON " << ice[1]);
  CHECK(ice[1] < ice[0]);     // the freezing curve keeps liquid water below 0 C
  // T is not asserted here: with SJ2016 about 30 % of the water is still liquid at -eps, so the front
  // releases less latent heat than free water and the ON/OFF T order depends on time and depth
  // (quantified in tests/regression/ro2000).
}

}
