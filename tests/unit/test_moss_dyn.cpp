// model_settings.moss_dyn: dynamic live-moss thickness (default 0 = legacy static 0.01 m).
// See include/MossDyn.h and docs/moss_dyn_worknotes.md.
// Refs: Yi et al. 2009 doi:10.1029/2008JG000841 (regrowth curve); Jafarov et al. 2013 doi:10.1088/1748-9326/8/3/035030.
#include "doctest.h"
#include "../../include/MossDyn.h"
#include "../../include/ModelData.h"
#include <json/json.h>

namespace {
Json::Value ms(std::initializer_list<std::pair<const char*, Json::Value>> kv) {
  Json::Value o(Json::objectValue); for (auto& p : kv) o[p.first] = p.second; return o;
}
}

TEST_CASE("moss_dyn: absent keys give mode 0 and documented defaults") {
  mossdyn::Params p;
  CHECK(mossdyn::parse(Json::Value(Json::objectValue), p) == "");
  CHECK(p.mode == 0); CHECK(p.dmax == -1.0); CHECK(p.dmin == 0.001);
  CHECK(p.t_half == 5.0); CHECK(p.bulkden == 25.0); CHECK(p.cfrac == 0.45);
  Json::Value c; c["model_settings"]["deep_carbon"] = 1;
  ModelData md(c);
  CHECK(md.moss.mode == 0);
}

TEST_CASE("moss_dyn: valid modes and parameters parse") {
  mossdyn::Params p;
  for (int m = 0; m <= 2; m++) { CHECK(mossdyn::parse(ms({{"moss_dyn", m}}), p) == ""); CHECK(p.mode == m); }
  CHECK(mossdyn::parse(ms({{"moss_dyn", 1}, {"moss_dmax", 0.05}, {"moss_t_half", 0.0}, {"moss_min_thick", 0.002}}), p) == "");
  CHECK(p.dmax == 0.05); CHECK(p.t_half == 0.0); CHECK(p.dmin == 0.002);
}

TEST_CASE("moss_dyn: invalid settings are rejected") {
  mossdyn::Params p;
  CHECK(mossdyn::parse(ms({{"moss_dyn", 3}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_dyn", -1}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_dyn", 1.5}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_dyn", true}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_dyn", "1"}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_min_thick", 0.0}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_min_thick", 0.1}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_dmax", 0.0005}}), p) != "");   // below moss_min_thick
  CHECK(mossdyn::parse(ms({{"moss_dmax", 2.0}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_t_half", -1.0}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_bulkden", 0.0}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_cfrac", 1.5}}), p) != "");
}

TEST_CASE("moss_dyn mode 1: regrowth curve d = dmax t/(t+t_half), floored and capped") {
  const double dmax = 0.04, dmin = 0.001, th = 5.0;
  CHECK(mossdyn::regrow_thick(dmax, dmin, 0.0, th) == dmin);                       // just removed
  CHECK(mossdyn::regrow_thick(dmax, dmin, th, th) == doctest::Approx(dmax / 2));   // half at t_half
  CHECK(mossdyn::regrow_thick(dmax, dmin, -9999.0, th) == dmax);                   // never disturbed
  CHECK(mossdyn::regrow_thick(dmax, dmin, 3.0, 0.0) == dmax);                      // t_half 0 = static dmax
  double prev = 0.0;
  for (int t = 0; t <= 500; t++) {
    double d = mossdyn::regrow_thick(dmax, dmin, t, th);
    CHECK(d >= prev); CHECK(d <= dmax); CHECK(d >= dmin); prev = d;
  }
  // equivalence case: dmax 0.01, t_half 0 gives exactly the legacy 0.01 m
  CHECK(mossdyn::regrow_thick(0.01, dmin, 0.0, 0.0) == 0.01);
}

TEST_CASE("moss_dyn mode 2: thickness from moss C, floored and capped") {
  // 25 kg/m3 x 0.45 = 11.25 kgC/m3 -> 112.5 gC/m2 per cm
  CHECK(mossdyn::carbon_thick(112.5, 25.0, 0.45, 0.001, 0.05) == doctest::Approx(0.01));
  CHECK(mossdyn::carbon_thick(0.0, 25.0, 0.45, 0.001, 0.05) == 0.001);
  CHECK(mossdyn::carbon_thick(-5.0, 25.0, 0.45, 0.001, 0.05) == 0.001);
  CHECK(mossdyn::carbon_thick(1.0e4, 25.0, 0.45, 0.001, 0.05) == 0.05);
}

// ---- moss/dyn-moss task 2: DOS-TEM dmax default and Chadburn et al. 2015 moss conductivity ----
#include "../../include/SoilLayer.h"
#include "../../include/physicalconst.h"
namespace {
struct MossTcGuard {
  int s = SoilLayer::moss_tc_scheme, o = SoilLayer::org_tc_scheme;
  double kd = SoilLayer::moss_tc_kdry, ks = SoilLayer::moss_tc_ksat;
  ~MossTcGuard() { SoilLayer::moss_tc_scheme = s; SoilLayer::org_tc_scheme = o;
                   SoilLayer::moss_tc_kdry = kd; SoilLayer::moss_tc_ksat = ks; }
};
// Moss-flagged layer, poro 0.95, water theta_liq + theta_ice = sat*poro split by frac_ice.
SoilLayer moss_at(double sat, double frac_ice = 0.0, bool moss = true) {
  SoilLayer l; l.isMoss = moss; l.isOrganic = !moss; l.isFibric = !moss;
  l.dz = 0.04; l.z = 0.0; l.poro = 0.95;
  const double w = sat * l.poro * DENLIQ * l.dz;
  l.liq = (1.0 - frac_ice) * w; l.ice = frac_ice * w * DENICE / DENLIQ;
  l.tcsatunf = 1.0; l.tcsatfrz = 2.0; l.tcdry = 0.05; l.tcmin = 0.05;
  return l;
}
}

TEST_CASE("moss_dyn dmax default: DOS-TEM 4 cm upland, 6 cm lowland (Jafarov et al. 2013)") {
  CHECK(mossdyn::default_dmax(0) == 0.04);
  CHECK(mossdyn::default_dmax(1) == 0.06);
  CHECK(mossdyn::default_dmax(-1) == 0.04);   // unknown -> upland
  mossdyn::Params p;
  CHECK(mossdyn::parse(Json::Value(Json::objectValue), p) == "");
  CHECK(p.dmax < 0.0);                        // absent key = use the default
  // regrowth with the default: half of dmax at t = t_half
  CHECK(mossdyn::regrow_thick(mossdyn::default_dmax(1), 0.001, 5.0, 5.0) == doctest::Approx(0.03));
  CHECK(mossdyn::regrow_thick(mossdyn::default_dmax(0), 0.001, 1e6, 5.0) == doctest::Approx(0.04).epsilon(1e-4));
}

TEST_CASE("moss_tc_scheme: parse, defaults and rejection") {
  mossdyn::Params p;
  CHECK(mossdyn::parse(Json::Value(Json::objectValue), p) == "");
  CHECK(p.tc_scheme == 0); CHECK(p.tc_kdry == 0.06); CHECK(p.tc_ksat == 0.5);
  CHECK(mossdyn::parse(ms({{"moss_tc_scheme", 1}}), p) == ""); CHECK(p.tc_scheme == 1);
  CHECK(mossdyn::parse(ms({{"moss_tc_scheme", 3}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_tc_scheme", 0.5}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_tc_scheme", true}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_tc_kdry", 0.6}}), p) != "");   // kdry > ksat
  CHECK(mossdyn::parse(ms({{"moss_tc_kdry", 0.0}}), p) != "");
  Json::Value c; c["model_settings"]["moss_tc_scheme"] = 1;
  ModelData md(c); CHECK(md.moss.tc_scheme == 1);
}

TEST_CASE("chadburn_k: linear 0.06 -> 0.5 W/m/K in theta/theta_sat, clamped") {
  CHECK(mossdyn::chadburn_k(0.0, 0.06, 0.5) == doctest::Approx(0.06));
  CHECK(mossdyn::chadburn_k(1.0, 0.06, 0.5) == doctest::Approx(0.5));
  CHECK(mossdyn::chadburn_k(0.5, 0.06, 0.5) == doctest::Approx(0.28));
  CHECK(mossdyn::chadburn_k(-0.2, 0.06, 0.5) == doctest::Approx(0.06));
  CHECK(mossdyn::chadburn_k(1.3, 0.06, 0.5) == doctest::Approx(0.5));
}

TEST_CASE("moss_tc_scheme 1: moss layer only, frozen and unfrozen, over org_tc_scheme 0 and 1") {
  MossTcGuard g;
  for (int ots : {0, 1}) {
    SoilLayer::org_tc_scheme = ots;
    SoilLayer::moss_tc_scheme = 0;
    const double leg_u = moss_at(0.5).getUnfThermCond(), leg_f = moss_at(0.5, 1.0).getFrzThermCond();
    const double fib_u = moss_at(0.5, 0.0, false).getUnfThermCond(), fib_f = moss_at(0.5, 1.0, false).getFrzThermCond();
    SoilLayer::moss_tc_scheme = 1; SoilLayer::moss_tc_kdry = 0.06; SoilLayer::moss_tc_ksat = 0.5;
    CHECK(moss_at(0.0).getUnfThermCond() == doctest::Approx(0.06));
    CHECK(moss_at(1.0).getUnfThermCond() == doctest::Approx(0.5));
    CHECK(moss_at(0.5).getUnfThermCond() == doctest::Approx(0.28));
    CHECK(moss_at(0.5, 1.0).getFrzThermCond() == doctest::Approx(0.28));   // ice counts as water content
    CHECK(moss_at(0.5, 0.5).getFrzThermCond() == doctest::Approx(0.28));
    CHECK(moss_at(0.5).getUnfThermCond() != doctest::Approx(leg_u));
    (void)leg_f;
    // non-moss organic layer unchanged
    CHECK(moss_at(0.5, 0.0, false).getUnfThermCond() == fib_u);
    CHECK(moss_at(0.5, 1.0, false).getFrzThermCond() == fib_f);
  }
}

// ---- moss/dyn-moss task 4: moss_water_scheme 1 (Chadburn 2015 Eq. 1) and moss_tc_scheme 2 (Johansen/Farouki) ----
namespace {
struct MossTc2Guard {
  int s = SoilLayer::moss_tc_scheme, o = SoilLayer::org_tc_scheme;
  double kd = SoilLayer::moss_tc_kdry, ks = SoilLayer::moss_tc_ksolid;
  double tf = SoilLayer::org_tc_thre_f, ku = SoilLayer::org_tc_kunf, kf = SoilLayer::org_tc_kfrz, od = SoilLayer::org_tc_dry,
         dep = SoilLayer::org_tc_depth;
  ~MossTc2Guard() { SoilLayer::moss_tc_scheme = s; SoilLayer::org_tc_scheme = o; SoilLayer::moss_tc_kdry = kd;
                    SoilLayer::moss_tc_ksolid = ks; SoilLayer::org_tc_thre_f = tf; SoilLayer::org_tc_kunf = ku;
                    SoilLayer::org_tc_kfrz = kf; SoilLayer::org_tc_dry = od; SoilLayer::org_tc_depth = dep; }
};
double mass(const mossdyn::WaterIn& w, const mossdyn::WaterOut& o) {   // residual of the moss balance
  return (o.liq - w.liq) - (o.intercept - o.drain + o.rise - o.evap);
}
}

TEST_CASE("moss_water_scheme / moss_tc_scheme 2: parse, defaults, rejection") {
  mossdyn::Params p;
  CHECK(mossdyn::parse(Json::Value(Json::objectValue), p) == "");
  CHECK(p.water_scheme == 0); CHECK(p.w_bsw == 1.0); CHECK(p.w_psisat == 0.12); CHECK(p.w_thetasat == -1.0);
  CHECK(p.tc_ksolid == 0.25);
  CHECK(mossdyn::parse(ms({{"moss_water_scheme", 1}, {"moss_tc_scheme", 2}, {"moss_tc_ksolid", 0.3}}), p) == "");
  CHECK(p.water_scheme == 1); CHECK(p.tc_scheme == 2); CHECK(p.tc_ksolid == 0.3);
  CHECK(mossdyn::parse(ms({{"moss_water_scheme", 2}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_tc_scheme", 3}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_tc_ksolid", 0.0}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_psisat", -1.0}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_thetasat", 1.5}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_thetasat", -0.5}}), p) != "");
  CHECK(mossdyn::parse(ms({{"moss_thetasat", 0.9}}), p) == ""); CHECK(p.w_thetasat == 0.9);
  CHECK(mossdyn::parse(ms({{"moss_bsw", 0.0}}), p) != "");
}

TEST_CASE("theta_eq: Chadburn 2015 Eq. 1 Brooks-Corey, b = 1, psi_sat 0.12 m, theta_sat 0.9") {
  CHECK(mossdyn::theta_eq(0.0, 0.9, 0.12, 1.0) == 0.9);
  CHECK(mossdyn::theta_eq(0.12, 0.9, 0.12, 1.0) == 0.9);
  CHECK(mossdyn::theta_eq(0.24, 0.9, 0.12, 1.0) == doctest::Approx(0.45));
  CHECK(mossdyn::theta_eq(1.2, 0.9, 0.12, 1.0) == doctest::Approx(0.09));
  CHECK(mossdyn::theta_eq(0.48, 0.9, 0.12, 2.0) == doctest::Approx(0.45));   // (4)^(-1/2)
}

TEST_CASE("water_step: mass balance closes and limits hold") {
  // 3 cm moss, poro 0.95, theta_sat 0.9 -> 27 mm capacity; minliq 0.05*0.95*30 = 1.425 mm
  for (double liq : {0.0, 3.0, 10.0, 26.0})
    for (double ice : {0.0, 5.0})
      for (double infil : {0.0, 2.0, 40.0})
        for (double theq : {0.05, 0.3, 0.9})
          for (double fav : {0.0, 3.0, 50.0})
            for (double ev : {0.0, 1.0, 30.0}) {
              mossdyn::WaterIn w{liq, ice, 0.03, 0.95, 0.9, theq, 1.425, infil, ev, fav, true};
              mossdyn::WaterOut o = mossdyn::water_step(w);
              CHECK(std::fabs(mass(w, o)) < 1e-12);
              CHECK(o.intercept >= 0.0); CHECK(o.intercept <= infil + 1e-12);
              CHECK(o.drain >= 0.0); CHECK(o.rise >= 0.0); CHECK(o.rise <= fav + 1e-12);
              CHECK(o.evap >= 0.0); CHECK(o.evap <= ev + 1e-12);
              CHECK(o.liq + ice * 1000.0 / 917.0 <= std::max(27.0, liq + ice * 1000.0 / 917.0) + 1e-9);
              CHECK(o.liq >= std::min(liq, 1.425) - 1e-12);
              CHECK(!(o.drain > 0.0 && o.rise > 0.0));
            }
  mossdyn::WaterIn w{5.0, 0.0, 0.03, 0.95, 0.9, 0.3, 1.425, 10.0, 1.0, 10.0, false};   // frozen: nothing
  mossdyn::WaterOut o = mossdyn::water_step(w);
  CHECK(o.liq == 5.0); CHECK(o.intercept == 0.0); CHECK(o.drain == 0.0); CHECK(o.rise == 0.0); CHECK(o.evap == 0.0);
  // rain fills to theta_sat then drains to theta_eq: 0 -> intercept 27, drain 18, liq 9 (theta 0.3)
  mossdyn::WaterIn r{0.0, 0.0, 0.03, 0.95, 0.9, 0.3, 1.425, 40.0, 0.0, 0.0, true};
  o = mossdyn::water_step(r);
  CHECK(o.intercept == doctest::Approx(27.0)); CHECK(o.drain == doctest::Approx(18.0)); CHECK(o.liq == doctest::Approx(9.0));
}

TEST_CASE("johansen_k (moss_tc_scheme 2): saturated end members, dry, Kersten") {
  const double eta = 0.95, ks = 0.25, kl = 0.6, ki = 2.29, kd = 0.06;
  const double ku = std::pow(ks, 1 - eta) * std::pow(kl, eta), kf = std::pow(ks, 1 - eta) * std::pow(ki, eta);
  CHECK(mossdyn::johansen_k(1.0, eta, ks, kl, ki, kd, false) == doctest::Approx(ku));
  CHECK(mossdyn::johansen_k(1.0, eta, ks, kl, ki, kd, true) == doctest::Approx(kf));
  CHECK(kf / ku == doctest::Approx(std::pow(ki / kl, eta)));           // GIPL (Nicolsky 2007) ratio, ~3.6
  CHECK(kf / ku > 3.0);
  CHECK(mossdyn::johansen_k(0.0, eta, ks, kl, ki, kd, false) == kd);
  CHECK(mossdyn::johansen_k(0.1, eta, ks, kl, ki, kd, false) == doctest::Approx(kd));   // Ke = 0 at Sr = 0.1
  CHECK(mossdyn::johansen_k(0.5, eta, ks, kl, ki, kd, false) == doctest::Approx((std::log10(0.5) + 1) * ku + (-std::log10(0.5)) * kd));
  CHECK(mossdyn::johansen_k(0.5, eta, ks, kl, ki, kd, true) == doctest::Approx(0.5 * kf + 0.5 * kd));
}

TEST_CASE("moss_tc_scheme 2: moss layer only, S21H (org_tc_scheme 1) not applied to moss") {
  MossTc2Guard g;
  SoilLayer::org_tc_scheme = 1; SoilLayer::org_tc_kunf = 0.27; SoilLayer::org_tc_kfrz = 0.5134;
  SoilLayer::org_tc_dry = 0.08571; SoilLayer::org_tc_thre_f = 0.1795; SoilLayer::org_tc_depth = 10.0;
  SoilLayer::moss_tc_scheme = 0;
  const double s21h_moss = moss_at(0.25).getUnfThermCond();
  const double fib_u = moss_at(0.25, 0.0, false).getUnfThermCond(), fib_f = moss_at(0.25, 1.0, false).getFrzThermCond();
  CHECK(s21h_moss == doctest::Approx(moss_at(0.25).getOrgThermCond(false)));   // S21H on moss with scheme 0
  SoilLayer::moss_tc_scheme = 2; SoilLayer::moss_tc_ksolid = 0.25; SoilLayer::moss_tc_kdry = 0.06;
  SoilLayer l = moss_at(0.25);
  const double sr = (l.getVolLiq() + l.getVolIce()) / l.poro;
  CHECK(l.getUnfThermCond() == doctest::Approx(mossdyn::johansen_k(sr, 0.95, 0.25, TCLIQ, TCICE, 0.06, false)));
  CHECK(l.getUnfThermCond() != doctest::Approx(s21h_moss));
  SoilLayer fz = moss_at(1.0, 1.0);
  CHECK(fz.getFrzThermCond() == doctest::Approx(std::pow(0.25, 0.05) * std::pow((double)TCICE, 0.95)).epsilon(1e-3));
  CHECK(moss_at(0.25, 0.0, false).getUnfThermCond() == fib_u);     // fibric keeps S21H
  CHECK(moss_at(0.25, 1.0, false).getFrzThermCond() == fib_f);
}

TEST_CASE("refreeze_take (moss_water_scheme 1, frozen moss): supply, pore space and cold content limits") {
  const double lf = 3.337e5;
  // 3 cm moss, poro 0.95: 28.5 mm pores; hcap 1e5 J/m2/K at -5 C -> cold 1.498 kg
  CHECK(mossdyn::refreeze_take(10.0, 0.03, 0.95, 1.0, 5.0, 1.0e5, -5.0, 0.0, lf) == doctest::Approx(5.0e5 / lf));
  CHECK(mossdyn::refreeze_take(1.0, 0.03, 0.95, 1.0, 5.0, 1.0e5, -5.0, 0.0, lf) == doctest::Approx(1.0));   // supply
  // nearly full of ice: room 28.5 - 1 - 25*1000/917 = 0.2370 mm -> 0.2173 kg as ice
  CHECK(mossdyn::refreeze_take(10.0, 0.03, 0.95, 1.0, 25.0, 1.0e6, -20.0, 0.0, lf)
        == doctest::Approx((28.5 - 1.0 - 25.0 * 1000.0 / 917.0) * 0.917));
  CHECK(mossdyn::refreeze_take(10.0, 0.03, 0.95, 1.0, 30.0, 1.0e6, -20.0, 0.0, lf) == 0.0);   // no room
  CHECK(mossdyn::refreeze_take(10.0, 0.03, 0.95, 1.0, 5.0, 1.0e5, 0.0, 0.0, lf) == 0.0);      // no cold content
  CHECK(mossdyn::refreeze_take(0.0, 0.03, 0.95, 1.0, 5.0, 1.0e5, -5.0, 0.0, lf) == 0.0);      // no supply
  // warm rain (h_in = 4188 * 2 C) needs more cold per kg
  CHECK(mossdyn::refreeze_take(10.0, 0.03, 0.95, 1.0, 5.0, 1.0e5, -5.0, 8376.0, lf) == doctest::Approx(5.0e5 / (lf + 8376.0)));
  // refreezing the take at Lf + h_in uses at most the cold content: the layer stays at <= 0 C
  for (double T : {-0.01, -1.0, -10.0})
    for (double sup : {0.5, 5.0, 50.0}) {
      const double hcap = 6.0e4, t = mossdyn::refreeze_take(sup, 0.03, 0.95, 0.5, 3.0, hcap, T, 0.0, lf);
      CHECK(t * lf <= hcap * (-T) * (1.0 + 1e-12));
      CHECK(t <= sup);
    }
}
