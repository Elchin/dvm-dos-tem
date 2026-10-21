// Area 2: Kersten number in SoilLayer::getUnfThermCond (kersten_form),
// Area 10: mineral thermal options and org_tc_scheme (Cote & Konrad 2005).
#include "doctest.h"
#include "test_helpers.h"

namespace {
// Restores the SoilLayer static switches after each test.
struct StaticGuard {
  int kf = SoilLayer::kersten_form, ots = SoilLayer::org_tc_scheme, mtf = SoilLayer::mineral_tcdry_form;
  double kunf = SoilLayer::org_tc_kunf, kfrz = SoilLayer::org_tc_kfrz, kdry = SoilLayer::org_tc_dry,
         tf = SoilLayer::org_tc_thre_f, th = SoilLayer::org_tc_thre_h, dep = SoilLayer::org_tc_depth,
         mts = SoilLayer::mineral_tcsolid, mvh = SoilLayer::mineral_vhcsolid;
  ~StaticGuard() {
    SoilLayer::kersten_form = kf; SoilLayer::org_tc_scheme = ots; SoilLayer::mineral_tcdry_form = mtf;
    SoilLayer::org_tc_kunf = kunf; SoilLayer::org_tc_kfrz = kfrz; SoilLayer::org_tc_dry = kdry;
    SoilLayer::org_tc_thre_f = tf; SoilLayer::org_tc_thre_h = th; SoilLayer::org_tc_depth = dep;
    SoilLayer::mineral_tcsolid = mts; SoilLayer::mineral_vhcsolid = mvh;
  }
};

const double KSAT = 2.0, KSATF = 3.0, KDRY = 0.25;

// Layer with saturation Sr = (theta_liq + theta_ice)/poro, all liquid.
SoilLayer layer_at_sr(double sr, double poro = 0.5, int solind = 5) {
  SoilLayer l;
  l.isMineral = true;
  l.dz = 0.1; l.z = 0.5; l.poro = poro;
  l.liq = sr * poro * DENLIQ * l.dz; l.ice = 0.0;
  l.tcsatunf = KSAT; l.tcsatfrz = KSATF; l.tcdry = KDRY; l.tcmin = KDRY;
  l.solind = solind;
  return l;
}
double expected_k(double ke) { return ke * KSAT + (1.0 - ke) * KDRY; }
}

TEST_SUITE("kersten") {

TEST_CASE("legacy ln form (kersten_form=0): Ke = 1 + ln(Sr), floored at 0") {
  StaticGuard g; SoilLayer::kersten_form = 0; SoilLayer::org_tc_scheme = 0;
  CHECK(layer_at_sr(1.0).getUnfThermCond() == doctest::Approx(KSAT));
  CHECK(layer_at_sr(std::exp(-1.0)).getUnfThermCond() == doctest::Approx(KDRY));   // Ke = 0 at Sr = 1/e
  CHECK(layer_at_sr(0.1).getUnfThermCond() == doctest::Approx(KDRY));              // Ke < 0 -> 0
  CHECK(layer_at_sr(0.0).getUnfThermCond() == doctest::Approx(KDRY));              // dry
  CHECK(layer_at_sr(0.5).getUnfThermCond() == doctest::Approx(expected_k(1.0 + std::log(0.5))));  // Ke = 0.3069
  CHECK(layer_at_sr(0.8).getUnfThermCond() == doctest::Approx(expected_k(1.0 + std::log(0.8))));
}

TEST_CASE("Johansen fine-soil log10 form (kersten_form=1): Ke = 1 + log10(Sr), floored at 0") {
  StaticGuard g; SoilLayer::kersten_form = 1; SoilLayer::org_tc_scheme = 0;
  CHECK(layer_at_sr(1.0).getUnfThermCond() == doctest::Approx(KSAT));
  CHECK(layer_at_sr(0.1).getUnfThermCond() == doctest::Approx(KDRY));              // Ke = 0 at Sr = 0.1
  CHECK(layer_at_sr(0.0).getUnfThermCond() == doctest::Approx(KDRY));
  // Sr = 1/e: Ke = 1 - log10(e) = 0.5657 (ln form gives 0 here)
  CHECK(layer_at_sr(std::exp(-1.0)).getUnfThermCond() == doctest::Approx(expected_k(1.0 - std::log10(std::exp(1.0)))));
  CHECK(layer_at_sr(std::exp(-1.0)).getUnfThermCond() == doctest::Approx(0.5657 * KSAT + 0.4343 * KDRY).epsilon(1e-4));
  CHECK(layer_at_sr(0.5).getUnfThermCond() == doctest::Approx(expected_k(1.0 + std::log10(0.5))));  // Ke = 0.699
}

TEST_CASE("log10 form is never below the ln form; both equal at Sr = 1") {
  StaticGuard g; SoilLayer::org_tc_scheme = 0;
  for (int i = 0; i <= 100; ++i) {
    double sr = i / 100.0;
    SoilLayer::kersten_form = 0; double kln = layer_at_sr(sr).getUnfThermCond();
    SoilLayer::kersten_form = 1; double klg = layer_at_sr(sr).getUnfThermCond();
    CHECK(klg >= kln - 1e-12);
  }
}

TEST_CASE("frozen conductivity uses Ke = Sr (linear in saturation)") {
  StaticGuard g; SoilLayer::org_tc_scheme = 0;
  for (double sr : {0.05, 0.1, std::exp(-1.0), 0.5, 0.9, 1.0}) {
    SoilLayer l = layer_at_sr(sr);
    l.ice = l.liq; l.liq = 0.0;   // all water frozen (same mass)
    double sr_eff = std::min(1.0, (l.ice / DENICE / l.dz) / l.poro);
    CHECK(l.getFrzThermCond() == doctest::Approx(sr_eff * KSATF + (1.0 - sr_eff) * KDRY));
  }
  SoilLayer dry = layer_at_sr(0.0);
  CHECK(dry.getFrzThermCond() == doctest::Approx(KDRY));
}

TEST_CASE("k stays within [kdry, ksat] and is non-decreasing in Sr (both forms, unfrozen and frozen)") {
  StaticGuard g; SoilLayer::org_tc_scheme = 0;
  for (int form = 0; form <= 1; ++form) {
    SoilLayer::kersten_form = form;
    double prevu = -1, prevf = -1;
    for (int i = 0; i <= 130; ++i) {          // includes Sr > 1 (capped at 1)
      double sr = i / 100.0;
      SoilLayer l = layer_at_sr(sr);
      double ku = l.getUnfThermCond();
      CHECK(ku >= KDRY - 1e-12); CHECK(ku <= KSAT + 1e-12);
      CHECK(ku >= prevu - 1e-12); prevu = ku;
      double kf = l.getFrzThermCond();
      CHECK(kf >= KDRY - 1e-12); CHECK(kf <= KSATF + 1e-12);
      CHECK(kf >= prevf - 1e-12); prevf = kf;
    }
  }
}

TEST_CASE("high-porosity layers are floored at tcmin (organic convection floor)") {
  StaticGuard g; SoilLayer::org_tc_scheme = 0; SoilLayer::kersten_form = 0;
  SoilLayer l = layer_at_sr(0.2, 0.92, 5);
  l.tcmin = 0.3;    // above kdry
  CHECK(l.getUnfThermCond() == doctest::Approx(0.3));
  SoilLayer m = layer_at_sr(0.2, 0.5, 5);
  m.tcmin = 0.3;    // poro < 0.8: no floor
  CHECK(m.getUnfThermCond() == doctest::Approx(KDRY));
}

TEST_CASE("derivePhysicalProperty: geometric-mean saturated/dry conductivities (Johansen 1975)") {
  StaticGuard g; SoilLayer::mineral_tcdry_form = 0;
  SoilLayer l; l.isMineral = true; l.dz = 0.2; l.poro = 0.4; l.tcsolid = 3.0; l.bulkden = 1620.;
  l.derivePhysicalProperty();
  CHECK(l.tcsatunf == doctest::Approx(std::pow(3.0, 0.6) * std::pow(0.6, 0.4)));
  CHECK(l.tcsatfrz == doctest::Approx(std::pow(3.0, 0.6) * std::pow(2.29, 0.4)));
  CHECK(l.tcdry == doctest::Approx(std::pow(3.0, 0.6) * std::pow(0.023, 0.4)).epsilon(1e-6));
  CHECK(l.maxliq == doctest::Approx(0.4 * 1000 * 0.2));
  CHECK(l.maxice == doctest::Approx(0.4 * 917 * 0.2));
  CHECK(l.minliq == doctest::Approx(0.05 * 0.4 * 1000 * 0.2));
}

} // TEST_SUITE kersten

TEST_SUITE("thermal_options") {

TEST_CASE("MineralLayer texture formulas (CLM4.5): porosity, tcsolid, vhcsolid") {
  StaticGuard g; SoilLayer::mineral_tcsolid = -1; SoilLayer::mineral_vhcsolid = -1;
  MineralLayer m(0.1, 40.0f, 40.0f, 20.0f);
  CHECK(m.poro == doctest::Approx(0.489 - 0.00126 * 40));
  CHECK(m.tcsolid == doctest::Approx((8.8 * 40 + 2.92 * 20) / 60.0));      // 6.84 W/m/K
  CHECK(m.vhcsolid == doctest::Approx((2.128 * 40 + 2.385 * 20) / 60.0 * 1e6));
  CHECK(m.bulkden == doctest::Approx(2700 * (1 - m.poro)));
}

TEST_CASE("mineral_tcsolid / mineral_vhcsolid replace the texture values") {
  StaticGuard g; SoilLayer::mineral_tcsolid = 4.044; SoilLayer::mineral_vhcsolid = 1.834e6;
  MineralLayer m(0.1, 40.0f, 40.0f, 20.0f);
  CHECK(m.tcsolid == doctest::Approx(4.044));
  CHECK(m.vhcsolid == doctest::Approx(1.834e6));
  CHECK(m.tcsatunf == doctest::Approx(std::pow(4.044, 1 - m.poro) * std::pow(0.6, m.poro)));
  CHECK(m.tcsatfrz == doctest::Approx(std::pow(4.044, 1 - m.poro) * std::pow(2.29, m.poro)));
}

TEST_CASE("mineral_tcdry_form=1: Johansen dry conductivity, ~0.21 W/m/K at rho_b = 1500 kg/m3") {
  StaticGuard g; SoilLayer::mineral_tcsolid = -1;
  MineralLayer m(0.1, 40.0f, 40.0f, 20.0f);
  m.bulkden = 1500.;
  SoilLayer::mineral_tcdry_form = 0; m.derivePhysicalProperty();
  double legacy = std::pow(m.tcsolid, 1 - m.poro) * std::pow(0.023, m.poro);
  CHECK(m.tcdry == doctest::Approx(legacy).epsilon(1e-6));
  SoilLayer::mineral_tcdry_form = 1; m.derivePhysicalProperty();
  CHECK(m.tcdry == doctest::Approx((0.135 * 1500 + 64.7) / (2700 - 0.947 * 1500)));
  CHECK(m.tcdry == doctest::Approx(0.2088).epsilon(1e-3));
  CHECK(m.tcmin == doctest::Approx(m.tcdry));
  // organic layers are not affected
  SoilLayer o; o.isOrganic = true; o.dz = 0.1; o.poro = 0.9; o.tcsolid = 0.25; o.bulkden = 80.;
  o.derivePhysicalProperty();
  CHECK(o.tcdry == doctest::Approx(std::pow(0.25, 0.1) * std::pow(0.023, 0.9)).epsilon(1e-6));
}

// Cote & Konrad (2005), CLASS organic soils: Ke = kappa*Sr/(1+(kappa-1)*Sr),
// kappa_unf = 0.6, kappa_frz = 0.25 for organic fibrous soils; kdry = 0.30 exp(-2 poro).
SoilLayer org_layer(double theta, bool humic = false, double z = 0.05) {
  SoilLayer l; l.isOrganic = true; l.isFibric = !humic; l.isHumic = humic;
  l.z = z; l.dz = 0.05; l.poro = 0.9;
  l.liq = theta * DENLIQ * l.dz; l.ice = 0.0;
  l.tcsatunf = 0.55; l.tcsatfrz = 2.0; l.tcdry = 0.05; l.tcmin = 0.05; l.solind = 3;
  return l;
}

TEST_CASE("org_tc_scheme=1 (Cote & Konrad): known values") {
  StaticGuard g;
  SoilLayer::org_tc_scheme = 1; SoilLayer::org_tc_kunf = 0.6; SoilLayer::org_tc_kfrz = 0.25;
  SoilLayer::org_tc_dry = -1.0; SoilLayer::org_tc_thre_f = 0.275; SoilLayer::org_tc_thre_h = 0.62;
  SoilLayer::org_tc_depth = 10.0;
  const double kdry = 0.30 * std::exp(-2.0 * 0.9);     // 0.0496 W/m/K
  CHECK(kdry == doctest::Approx(0.04959).epsilon(1e-3));
  // Sr = 0 -> kdry ; Sr >= 1 (theta >= theta_re) -> ksat
  CHECK(org_layer(0.0).getUnfThermCond() == doctest::Approx(kdry));
  CHECK(org_layer(0.275).getUnfThermCond() == doctest::Approx(0.55));
  CHECK(org_layer(0.6).getUnfThermCond() == doctest::Approx(0.55));
  // Sr = 0.5: Ke = 0.6*0.5/(1-0.4*0.5) = 0.375
  CHECK(org_layer(0.1375).getUnfThermCond() == doctest::Approx(0.375 * 0.55 + 0.625 * kdry));
  // frozen, Sr = 0.5: Ke = 0.25*0.5/(1-0.75*0.5) = 0.2, uses tcsatfrz
  SoilLayer f = org_layer(0.0); f.ice = 0.1375 * DENICE * f.dz;
  CHECK(f.getFrzThermCond() == doctest::Approx(0.2 * 2.0 + 0.8 * kdry));
  // humic layers use theta_re = 0.62
  CHECK(org_layer(0.31, true).getUnfThermCond() == doctest::Approx(0.375 * 0.55 + 0.625 * kdry));
  // kappa = 1 reduces to linear interpolation (JULES moss, Chadburn et al. 2015)
  SoilLayer::org_tc_kunf = 1.0;
  CHECK(org_layer(0.1375).getUnfThermCond() == doctest::Approx(0.5 * 0.55 + 0.5 * kdry));
  // user dry conductivity
  SoilLayer::org_tc_dry = 0.08571;
  CHECK(org_layer(0.0).getUnfThermCond() == doctest::Approx(0.08571));
}

TEST_CASE("org_tc_scheme: not applied to mineral layers, below org_tc_depth, or when off") {
  StaticGuard g; SoilLayer::kersten_form = 0;
  SoilLayer::org_tc_scheme = 1; SoilLayer::org_tc_kunf = 0.6; SoilLayer::org_tc_dry = -1.0;
  SoilLayer::org_tc_thre_f = 0.275; SoilLayer::org_tc_depth = 0.1;
  SoilLayer deep = org_layer(0.1375, false, 0.2);      // z >= org_tc_depth
  CHECK(deep.getOrgThermCond(false) < 0.0);
  SoilLayer::org_tc_scheme = 0; double legacy = deep.getUnfThermCond();
  SoilLayer::org_tc_scheme = 1; CHECK(deep.getUnfThermCond() == doctest::Approx(legacy));
  SoilLayer mine = layer_at_sr(0.5);
  CHECK(mine.getOrgThermCond(false) < 0.0);
  SoilLayer::org_tc_scheme = 0;
  CHECK(org_layer(0.1375).getOrgThermCond(false) < 0.0);
}

TEST_CASE("org_tc_scheme: theta_re larger than porosity falls back to porosity") {
  StaticGuard g; SoilLayer::org_tc_scheme = 1; SoilLayer::org_tc_kunf = 1.0; SoilLayer::org_tc_dry = 0.05;
  SoilLayer::org_tc_thre_f = 0.95; SoilLayer::org_tc_depth = 10.;
  SoilLayer l = org_layer(0.45);           // Sr = 0.45/0.9 = 0.5
  CHECK(l.getUnfThermCond() == doctest::Approx(0.5 * 0.55 + 0.5 * 0.05));
}

TEST_CASE("org_tc_scheme: k within [kdry, ksat] for all moisture (unfrozen and frozen)") {
  StaticGuard g; SoilLayer::org_tc_scheme = 1; SoilLayer::org_tc_kunf = 0.8659; SoilLayer::org_tc_kfrz = 0.5134;
  SoilLayer::org_tc_dry = 0.08571; SoilLayer::org_tc_thre_f = 0.1795; SoilLayer::org_tc_depth = 10.;
  double prev = 0;
  for (int i = 0; i <= 90; ++i) {
    SoilLayer l = org_layer(i / 100.0);
    double k = l.getUnfThermCond();
    CHECK(k >= 0.08571 - 1e-12); CHECK(k <= 0.55 + 1e-12); CHECK(k >= prev - 1e-12); prev = k;
  }
}

} // TEST_SUITE thermal_options
