// v5k: switchable fixes for the known bugs B1-B4 (all default OFF = legacy).
// Each fix has tests that pass with its switch ON; the original expected-fail
// tests (test_climate.cpp B1, test_stefan_water.cpp B2/B3, test_seb.cpp B4) are
// kept unchanged and still document the legacy (switch OFF) behaviour.
//   B1 model_settings.svp_ice_fix        -> Climate::svp_ice_fix (svp + SEB esat)
//   B2 model_settings.richards_infil_fix -> Soil_Env::richards_infil_fix
//   B3 model_settings.stefan_fronts_fix  -> Stefan::fronts_fix
//   B4 model_settings.snow_albedo_fix    -> Snow_Env::albedo_fix
#include "doctest.h"
#include "test_helpers.h"
#include <unistd.h>

using ut::make_soil_layer;

namespace {
struct FixRestore {
  int svp, fr, alb, vp;
  FixRestore() : svp(Climate::svp_ice_fix), fr(Stefan::fronts_fix), alb(Snow_Env::albedo_fix), vp(Climate::vp_hpa_fix) {}
  ~FixRestore() { Climate::svp_ice_fix = svp; Stefan::fronts_fix = fr; Snow_Env::albedo_fix = alb; Climate::vp_hpa_fix = vp; }
};
double wmo_ice_pa(double t) { return 611.2 * std::exp(22.46 * t / (272.62 + t)); }

struct SnowSeb {   // same set-up as test_seb.cpp SebFixture(snow = true)
  ut::Column col; EnvData ed; CohortData cd; Snow_Env snowenv; Soil_Env se;
  SnowSeb()
    : col(std::vector<SoilLayer*>{make_soil_layer(0.0, 0.05, 0.6, 15, 0, -2.0, 1, 1),
                                  make_soil_layer(0.05, 0.1, 0.5, 20, 0, -3.0, 1, 1),
                                  make_soil_layer(0.15, 0.2, 0.4, 30, 0, -4.0, 1, 1)}) {
    SnowLayer* sn = new SnowLayer(); sn->dz = 0.3; sn->rho = 250.0; sn->ice = 75.0; sn->liq = 0.0; sn->tem = -10.0;
    sn->z = 0.0; sn->frozen = 1; col.g.insertFront(sn);
    se.setGround(&col.g); se.setEnvData(&ed); se.setCohortData(&cd);
    snowenv.snowenvpar.albmax = 0.8; snowenv.snowenvpar.albmin = 0.4; se.snowenvp = &snowenv;
    se.seb_scheme = 2; se.seb_z0_snow = 0.001; se.seb_z0_soil = 0.01; se.seb_zwind = 10.0; se.seb_ztemp = 2.0;
    se.seb_emis_snow = 0.99; se.seb_emis_soil = 0.97; se.seb_wind = -1.0; se.seb_minwind = 0.5;
    se.seb_stab_b = 10.0; se.seb_tran_frac = 1.0; se.seb_pres = 100000.0; se.seb_fc_last = 0.0;
    cd.m_vegd.fpcsum = 0.0;
    ed.d_v2g.swthfl = 0.0; ed.d_a2l.nirr = 0.0; ed.seb_lwdown = 230.0; ed.seb_wind = 3.0;
    ed.d_atmd.vp = 0.8; ed.d_soi2a.evap = 0.0; ed.d_v2a.tran = 0.0;   // vp in hPa (vp_hpa_fix 0)
  }
  double residual() const {
    const SebDiag& d = se.sebdiag;
    return (1.0 - d.alb) * d.swin + d.lwin - d.lwout - d.h - d.le - d.g - d.qmelt;
  }
};

std::vector<SoilLayer*> alternating(int n) {   // n x 0.1 m layers, frozen / unfrozen alternating
  std::vector<SoilLayer*> v;
  for (int i = 0; i < n; ++i) {
    bool fz = (i % 2 == 1);
    v.push_back(make_soil_layer(0.1 * i, 0.1, 0.4, fz ? 0.0 : 0.5 * 0.4 * DENLIQ * 0.1,
                                fz ? 0.5 * 0.4 * DENICE * 0.1 : 0.0, fz ? -1.0 : 1.0, fz ? 1 : -1, fz ? 1.0 : 0.0));
  }
  return v;
}
std::vector<SoilLayer*> wet3() {
  std::vector<SoilLayer*> v;
  for (int i = 0; i < 3; ++i) v.push_back(make_soil_layer(0.1 * i, 0.1, 0.4, 0.5 * 0.4 * DENLIQ * 0.1, 0.0, 2.0, -1, 0.0));
  return v;
}
}

TEST_SUITE("fixes_v5k") {

TEST_CASE("all B1-B4 switches default to 0 (legacy)") {
  CHECK(Climate::svp_ice_fix == 0);
  CHECK(Stefan::fronts_fix == 0);
  CHECK(Snow_Env::albedo_fix == 0);
  Soil_Env se; CHECK(se.richards_infil_fix == 0);
}

// ------------------------------------------------------------------ B1
TEST_CASE("B1 svp_ice_fix=1: svp over ice matches WMO / Murphy-Koop (-20 C 103.2 Pa, -40 C 12.84 Pa)") {
  FixRestore r; Climate::svp_ice_fix = 1;
  CHECK(calculate_saturated_vapor_pressure(-20.0f) == doctest::Approx(103.24).epsilon(0.005));
  CHECK(calculate_saturated_vapor_pressure(-40.0f) == doctest::Approx(12.84).epsilon(0.01));
  CHECK(calculate_saturated_vapor_pressure(-5.0f) == doctest::Approx(wmo_ice_pa(-5.0)).epsilon(1e-5));
}

TEST_CASE("B1 svp_ice_fix=1: continuous at 0 C, unchanged above 0 C, ice < water below 0 C") {
  FixRestore r;
  Climate::svp_ice_fix = 0; float w20 = calculate_saturated_vapor_pressure(20.0f);
  float leg = calculate_saturated_vapor_pressure(-10.0f);
  Climate::svp_ice_fix = 1;
  CHECK(calculate_saturated_vapor_pressure(0.0f) == doctest::Approx(611.2).epsilon(1e-6));
  CHECK(calculate_saturated_vapor_pressure(-1e-3f) == doctest::Approx(611.2).epsilon(1e-3));
  CHECK(calculate_saturated_vapor_pressure(20.0f) == w20);
  CHECK(calculate_saturated_vapor_pressure(-10.0f) < leg);
  CHECK(calculate_saturated_vapor_pressure(-10.0f) == doctest::Approx(259.9).epsilon(0.005));
}

TEST_CASE("B1 svp_ice_fix=1: daily svp/vpd series use the ice formula") {
  FixRestore r; Climate::svp_ice_fix = 1; Climate::vp_hpa_fix = 1;
  Climate a;
  a.tair_d = std::vector<float>(1, -20.0f); a.vapo_d = std::vector<float>(1, 0.5f);   // hPa
  a.prec_d = std::vector<float>(1, 0.0f); a.girr_d = std::vector<float>(1, 100.0f);
  a.nirr_d = std::vector<float>(1, 50.0f);
  a.derive_daily_secondary_vars();
  CHECK(a.svp_d[0] == doctest::Approx(wmo_ice_pa(-20.0)).epsilon(1e-5));
  CHECK(a.vpd_d[0] == doctest::Approx(wmo_ice_pa(-20.0) - 50.0).epsilon(1e-4));
}

TEST_CASE("B1 svp_ice_fix=1: SEB snow sublimation uses esat over ice (less LE, balance still closes)") {
  FixRestore r; Climate::vp_hpa_fix = 0;
  SnowSeb a; Climate::svp_ice_fix = 0; double ta = a.se.solveSEB(-15.0);
  SnowSeb b; Climate::svp_ice_fix = 1; double tb = b.se.solveSEB(-15.0);
  CHECK(std::fabs(a.residual()) < 1e-4);
  CHECK(std::fabs(b.residual()) < 1e-4);
  CHECK(b.se.sebdiag.le < a.se.sebdiag.le);          // smaller esat(Ts) -> less sublimation
  CHECK(b.se.sebdiag.le > 0.0);                      // still sublimating: ea 80 Pa < esat_ice(Ts)
  CHECK(tb >= ta);                                   // less latent cooling -> warmer (or equal) skin
}

// ------------------------------------------------------------------ B2
TEST_CASE("B2 richards_infil_fix: Richards flags the early return") {
  ut::Column c(wet3()); Richards rr; double trans[MAX_SOI_LAY] = {0};
  rr.update(c.L[0], c.L[2], 0.0, 0.0, 1.0, trans, 0.0, 5.0 / SEC_IN_DAY, 0.0);
  CHECK(rr.early_return);
  rr.update(c.L[0], c.L[2], 0.3, 0.0, 1.0, trans, 0.0, 0.0, 0.0);
  CHECK_FALSE(rr.early_return);
}

TEST_CASE("B2 richards_infil_fix=1: early-return infiltration is not lost (puddle to ponding_max, rest runoff)") {
  ut::Column c(wet3()); Richards rr; double trans[MAX_SOI_LAY] = {0};
  EnvData ed; Soil_Env se; se.setEnvData(&ed); se.richards_infil_fix = 1;
  ed.d_soi2l.magic_puddle = 1.0; ed.d_soi2l.qover = 0.0; ed.d_soi2l.qinfl = 5.0;
  double w0 = c.total_water() + ed.d_soi2l.magic_puddle;
  rr.update(c.L[0], c.L[2], 0.0, 0.0, 1.0, trans, 0.2 / SEC_IN_DAY, 5.0 / SEC_IN_DAY, 0.0);
  REQUIRE(rr.early_return);
  se.rebookUnappliedRichards(5.0 / SEC_IN_DAY, 0.2 / SEC_IN_DAY);
  CHECK(ed.d_soi2l.magic_puddle == doctest::Approx(se.ponding_max_mm));
  CHECK(ed.d_soi2l.qover == doctest::Approx(5.0 - (se.ponding_max_mm - 1.0)));
  CHECK(ed.d_soi2l.qinfl == 0.0);
  CHECK(se.evap_unapplied_mm == doctest::Approx(0.2));   // evaporation not applied -> not reported
  // closure: soil + puddle + runoff = initial + input
  CHECK(c.total_water() + ed.d_soi2l.magic_puddle + ed.d_soi2l.qover == doctest::Approx(w0 + 5.0));
}

TEST_CASE("B2 richards_infil_fix=1: over-full puddle sends all infiltration to runoff; no input -> no change") {
  EnvData ed; Soil_Env se; se.setEnvData(&ed);
  ed.d_soi2l.magic_puddle = se.ponding_max_mm + 3.0; ed.d_soi2l.qover = 0.5;
  se.rebookUnappliedRichards(2.0 / SEC_IN_DAY, 0.0);
  CHECK(ed.d_soi2l.magic_puddle == doctest::Approx(se.ponding_max_mm + 3.0));
  CHECK(ed.d_soi2l.qover == doctest::Approx(2.5));
  double p = ed.d_soi2l.magic_puddle, q = ed.d_soi2l.qover;
  se.rebookUnappliedRichards(0.0, 0.0);
  CHECK(ed.d_soi2l.magic_puddle == p); CHECK(ed.d_soi2l.qover == q); CHECK(se.evap_unapplied_mm == 0.0);
}

// ------------------------------------------------------------------ B3
TEST_CASE("B3 stefan_fronts_fix=1: combineExtraFronts terminates, <= MAX_NUM_FNT, sorted, alternating") {
  int st_code = ut::exit_status_of([]() {
    alarm(3);
    Stefan::fronts_fix = 1;
    ut::Column c(alternating(5));
    for (int i = 0; i < MAX_NUM_FNT + 4; ++i) { c.g.frontsz.push_back(0.01 + 0.03 * i + (i % 3) * 0.001); c.g.frontstype.push_back(i % 2 ? 1 : -1); }
    Stefan st; st.setGround(&c.g); st.combineExtraFronts();
    bool ok = (int)c.g.frontsz.size() <= MAX_NUM_FNT && c.g.frontsz.size() == c.g.frontstype.size();
    for (size_t i = 1; i < c.g.frontsz.size(); ++i) {
      ok = ok && c.g.frontsz[i] > c.g.frontsz[i - 1] && c.g.frontstype[i] == -c.g.frontstype[i - 1];
    }
    _exit(ok ? 0 : 1);
  });
  CHECK(st_code == 0);
}

TEST_CASE("B3 stefan_fronts_fix=1: removes the closest pair first (single surplus pair: same as legacy)") {
  FixRestore r;
  for (int fx = 0; fx <= 1; ++fx) {
    Stefan::fronts_fix = fx; CAPTURE(fx);
    ut::Column c(alternating(5));
    for (int i = 0; i < MAX_NUM_FNT + 2; ++i) { c.g.frontsz.push_back(0.1 * (i + 1)); c.g.frontstype.push_back(i % 2 ? 1 : -1); }
    c.g.frontsz[5] = c.g.frontsz[4] + 0.001;            // closest pair: indices 4,5
    double keep = c.g.frontsz[6];
    Stefan st; st.setGround(&c.g); st.combineExtraFronts();
    CHECK((int)c.g.frontsz.size() == MAX_NUM_FNT);
    CHECK(c.g.frontsz[4] == keep);
  }
}

TEST_CASE("B3 stefan_fronts_fix=1: rebuildFrontsFromLayerStates with > MAX_NUM_FNT+2 boundaries terminates") {
  int st_code = ut::exit_status_of([]() {
    alarm(3);
    Stefan::fronts_fix = 1;
    ut::Column c(alternating(MAX_NUM_FNT + 6));
    Stefan st; st.setGround(&c.g); st.rebuildFrontsFromLayerStates();
    _exit((int)c.g.frontsz.size() <= MAX_NUM_FNT ? 0 : 1);
  });
  CHECK(st_code == 0);
}

// ------------------------------------------------------------------ B4
TEST_CASE("B4 snow_albedo_fix=1: albmax below -10 C, continuous, unchanged above -10 C") {
  FixRestore r; Snow_Env se; se.snowenvpar.albmax = 0.8; se.snowenvpar.albmin = 0.4;
  Snow_Env::albedo_fix = 0; double leg_m5 = se.getAlbedoVis(-5.0); CHECK(se.getAlbedoVis(-30.0) == doctest::Approx(0.9));
  Snow_Env::albedo_fix = 1;
  CHECK(se.getAlbedoVis(-30.0) == doctest::Approx(0.8));
  CHECK(se.getAlbedoVis(-10.001) == doctest::Approx(se.getAlbedoVis(-10.0)).epsilon(1e-3));
  CHECK(se.getAlbedoVis(-5.0) == leg_m5);
  CHECK(se.getAlbedoNir(-30.0) == se.getAlbedoVis(-30.0));
  se.snowenvpar.albmax = 0.85;
  CHECK(se.getAlbedoVis(-25.0) == doctest::Approx(0.85));
}

} // TEST_SUITE fixes_v5k

// ------------------------------------------------------------------ N1-N3 (unreported N sinks/sources)
// Soil N bookkeeping of one Soil_Bgc::deltastate() call: with no N input, no leaching
// and no fire, the soil N change must equal litterfall N - root N extraction.
namespace {
struct SoilNFix {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  SoilNFix(const std::vector<int>& types, const std::vector<double>& dz) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.bgcpar.som2co2 = 0.2; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    cd.m_soil.numsl = (int)types.size(); cd.m_soil.deepthick = 0.1; sb.totdzliq = 0.0;   // no leaching
    double z = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      bool on = i < (int)types.size();
      cd.m_soil.type[i] = on ? types[i] : 3; cd.m_soil.dz[i] = on ? dz[i] : 0.1; cd.m_soil.z[i] = z; z += cd.m_soil.dz[i];
      sb.tmp_sois.rawc[i] = 100; sb.tmp_sois.soma[i] = 50; sb.tmp_sois.sompr[i] = 200; sb.tmp_sois.somcr[i] = 300;
      sb.tmp_sois.orgn[i] = 20; sb.tmp_sois.avln[i] = 1;
      sb.del_soi2a.rhrawc[i] = 3.0; sb.del_soi2a.rhsoma[i] = 0.5; sb.del_soi2a.rhsompr[i] = 0.1; sb.del_soi2a.rhsomcr[i] = 0.01;
      sb.ltrflc[i] = on ? 2.0 : 0.0; sb.ltrfln[i] = on ? 0.1 : 0.0; sb.rtnextract[i] = on ? 0.02 : 0.0;
      sb.del_soi2soi.netnmin[i] = on ? 0.05 : 0.0; bd.m_v2soi.rtlfalfrac[i] = on ? 0.2 : 0.0;
      ed.m_sois.liq[i] = 10.0;
    }
    sb.tmp_sois.wdebrisc = 0; sb.tmp_sois.wdebrisn = 0; sb.d2wdebrisc = 0; sb.d2wdebrisn = 0; sb.del_soi2a.rhwdeb = 0;
    bd.m_a2soi.orgninput = 0; bd.m_a2soi.avlninput = 0; sb.del_soi2l.avlnlost = 0; ed.m_sois.draindepth = 0;
  }
  double dN() { double s = 0; for (int i = 0; i < cd.m_soil.numsl; ++i) s += sb.del_sois.orgn[i] + sb.del_sois.avln[i]; return s; }
  double expected() { double s = 0; for (int i = 0; i < cd.m_soil.numsl; ++i) s += sb.ltrfln[i] - sb.rtnextract[i]; return s; }
};
const std::vector<int> PROFILE = {1, 1, 2, 2, 3, 3};                 // fibric x2, humic x2 (top 0.1 m), mineral x2
const std::vector<double> PROFILE_DZ = {0.05, 0.05, 0.05, 0.05, 0.1, 0.1};
}

TEST_CASE("N2 n_som_move_fix=1: SOM N movement conserves N (soil dN = litter N - extraction)") {
  SoilNFix f(PROFILE, PROFILE_DZ); f.sb.n_som_move_fix = 1;
  f.sb.deltastate();
  CHECK(f.dN() == doctest::Approx(f.expected()).epsilon(1e-9));
  // (c) equal fibric layers give equal N (legacy: the second gives N for both layers' C)
  CHECK(f.sb.del_sois.orgn[0] == doctest::Approx(f.sb.del_sois.orgn[1]).epsilon(1e-12));
}

// KNOWN BUG N2 (pre-existing, found by the v5k N trace): legacy deltastate() loses the
// shallow->deep organic N transfer (see Soil_Bgc.h n_som_move_fix).
TEST_CASE("N2 legacy (n_som_move_fix=0): SOM N movement conserves N [KNOWN BUG: xfail]" * doctest::should_fail(true)) {
  SoilNFix f(PROFILE, PROFILE_DZ); f.sb.n_som_move_fix = 0;
  f.sb.deltastate();
  CHECK(f.dN() == doctest::Approx(f.expected()).epsilon(1e-9));
}

TEST_CASE("N2 legacy: documents the size of the loss = shallow->deep N transfer (C:N of the donors)") {
  SoilNFix a(PROFILE, PROFILE_DZ); a.sb.n_som_move_fix = 0; a.sb.deltastate();
  CHECK(a.dN() < a.expected() - 1e-6);   // N disappears
}

TEST_CASE("N1 n_rootless_fix: rootless layers conserve N (1 = stays organic, 2 = to AVLN); legacy loses netnmin") {
  const std::vector<int> MIN = {3, 3, 3};
  const std::vector<double> MDZ = {0.1, 0.1, 0.1};
  for (int fx = 0; fx <= 2; ++fx) {
    CAPTURE(fx);
    SoilNFix f(MIN, MDZ); f.sb.n_rootless_fix = fx; f.sb.n_som_move_fix = 1;
    f.bd.m_v2soi.rtlfalfrac[2] = 0.0; f.sb.rtnextract[2] = 0.0;       // layer 2 has no roots
    f.sb.deltastate();
    if (fx == 0) {
      CHECK(f.dN() == doctest::Approx(f.expected() - 0.05).epsilon(1e-9));   // legacy: layer-2 netnmin lost
    } else {
      CHECK(f.dN() == doctest::Approx(f.expected()).epsilon(1e-9));
    }
    if (fx == 1) { CHECK(f.sb.del_sois.avln[2] == 0.0); CHECK(f.sb.del_soi2soi.netnmin[2] == 0.0); }
    if (fx == 2) { CHECK(f.sb.del_sois.avln[2] == doctest::Approx(0.05)); }
  }
}

TEST_CASE("N3 n_nonvasc_extract=1: non-vascular N uptake is extracted from moss/fibric layers (avln*liq weights)") {
  for (int fx = 0; fx <= 1; ++fx) {
    CAPTURE(fx);
    CohortData cd; BgcData bd; EnvData ed; Vegetation_Bgc vb;
    vb.setCohortData(&cd); vb.setBgcData(&bd); vb.setEnvData(&ed); vb.set_nfeed(true);
    vb.ipft = 0; cd.m_veg.nonvascular[0] = 1; vb.n_nonvasc_extract = fx;
    cd.m_soil.numsl = 4; int ty[4] = {0, 1, 2, 3}; double av[4] = {0.1, 0.2, 0.3, 0.4};
    for (int i = 0; i < MAX_SOI_LAY; ++i) { cd.m_soil.type[i] = i < 4 ? ty[i] : 3; bd.m_sois.avln[i] = i < 4 ? av[i] : 0; ed.m_sois.liq[i] = 10; vb.fracnuptake[i] = 0; }
    for (int i = 0; i < NUM_PFT_PART; ++i) { bd.m_vegs.c[i] = 1; bd.m_vegs.strn[i] = 0.1; bd.m_soi2v.snuptake[i] = 0; bd.m_v2soi.ltrfalc[i] = 0; bd.m_v2soi.ltrfaln[i] = 0; bd.m_a2v.ingpp[i] = 0; bd.m_a2v.innpp[i] = 0; bd.m_a2v.gpp[i] = 0; bd.m_a2v.npp[i] = 0; bd.m_v2a.rg[i] = 0; bd.m_v2a.rm[i] = 0; bd.m_v2v.nmobil[i] = 0; bd.m_v2v.nresorb[i] = 0; }
    bd.m_vegs.labn = 0.1; bd.m_soi2v.lnuptake = 0.3;
    vb.afterIntegration();
    double tot = 0; for (int i = 0; i < 4; ++i) tot += bd.m_soi2v.nextract[i];
    if (fx == 1) {
      CHECK(tot == doctest::Approx(0.3));
      CHECK(bd.m_soi2v.nextract[0] == doctest::Approx(0.1)); CHECK(bd.m_soi2v.nextract[1] == doctest::Approx(0.2));
      CHECK(bd.m_soi2v.nextract[2] == 0.0);
    } else {
      CHECK(tot == 0.0);   // legacy: uptake 0.3 enters the plant, nothing leaves the soil
    }
  }
}
