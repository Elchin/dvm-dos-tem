// v5q-gpp-fleaf: model_settings.gpp_fleaf_mode (0 = legacy, default, bit-identical: GPP has no
// leaf-phenology factor; 1 = gpp_all *= fleaf for every PFT; 2 = only for non-woody PFTs).
// Refs: Raich et al. 1991, doi:10.2307/1941899 (TEM GPP with f(PHENOLOGY)); McGuire et al. 1992,
//   doi:10.1029/92GB00219. Using the dvmdostem fleaf (normalised transpiration recursion,
//   Vegetation::phenology) as that factor: no published source.
//   Caveat: fleaf as a direct GPP multiplier is the original TEM formulation (EET-based relative leaf
//   area, Raich 1991 / McGuire 1992, for vascular vegetation types); applying it to the dvmdostem PFT
//   set, in particular non-vascular PFTs (Sphagnum, lichen) and CMT51 sedges, has NO published source
//   of its own (calibration choice, step22 / GPP_RECO_OPTIONS.md). DOIs checked on Crossref 2026-11-02.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>
#include <limits>

namespace {
Json::Value gfmsg(Json::Value mode) { Json::Value c; c["model_settings"]["gpp_fleaf_mode"] = mode; return c; }

// One PFT, one Vegetation_Bgc::prepareIntegration + delta() call (GPP is computed in delta()), non-limiting state.
struct GfVeg {
  CohortData cd; BgcData bd; EnvData ed; Vegetation_Bgc vb;
  GfVeg(int mode, int woody, double fleaf) {
    vb.ipft = 0; vb.setCohortData(&cd); vb.setBgcData(&bd); vb.setEnvData(&ed);
    vb.gpp_fleaf_mode = mode;
    ed.m_l2a.pet = 60.0; ed.m_l2a.eet = 30.0; ed.m_atms.co2 = 400.0; ed.m_atms.ta = 12.0;
    ed.m_a2l.par = 120.0; ed.m_soid.rtdpthawpct = 1.0;
    cd.m_veg.ifwoody[0] = woody; cd.m_veg.ifdeciwoody[0] = 0; cd.m_veg.ifperenial[0] = 1;
    cd.m_vegd.fleaf[0] = fleaf; cd.m_vegd.ffoliage[0] = 1.0; cd.m_vegd.maxleafc[0] = 50.0;
    cd.mthsdist = 1000; bd.m_vegs.deadc0 = 0; bd.m_vegs.deadc = 0;
    vb.bgcpar.tmin = -1; vb.bgcpar.tmax = 30; vb.bgcpar.toptmin = 10; vb.bgcpar.toptmax = 20;
    vb.bgcpar.ki = 75; vb.bgcpar.kc = 400; vb.bgcpar.raq10a0 = 2.35665; vb.bgcpar.raq10a1 = -0.053077;
    vb.bgcpar.raq10a2 = 0.00229; vb.bgcpar.raq10a3 = -0.0000412;
    vb.calpar.cmax = 300; vb.calpar.frg = 0.2;
    for (int i = 0; i < NUM_PFT_PART; ++i) {
      vb.bgcpar.cpart[i] = (i == I_leaf) ? 0.5 : 0.25; vb.calpar.kra = -8.06e-05;
      vb.calpar.krb[i] = -6.0; bd.m_vegs.c[i] = 10.0; vb.tmp_vegs.c[i] = 10.0;
    }
  }
  double ingpp() { vb.prepareIntegration(false); vb.delta(); double s = 0; for (int i = 0; i < NUM_PFT_PART; ++i) s += vb.del_a2v.ingpp[i]; return s; }
};
}

TEST_SUITE("gppfleaf_v5q") {

TEST_CASE("gpp_fleaf_mode default 0: ModelData(), ModelData(Json) without the key, Vegetation_Bgc()") {
  ModelData md; CHECK(md.gpp_fleaf_mode == 0);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c); CHECK(mj.gpp_fleaf_mode == 0);
  Vegetation_Bgc vb; CHECK(vb.gpp_fleaf_mode == 0);
}

TEST_CASE("gpp_fleaf_mode validation: integers 0..2 accepted, anything else FATAL") {
  CHECK(ModelData(gfmsg(0)).gpp_fleaf_mode == 0);
  CHECK(ModelData(gfmsg(1)).gpp_fleaf_mode == 1);
  CHECK(ModelData(gfmsg(2)).gpp_fleaf_mode == 2);
  CHECK(ut::exit_status_of([] { ModelData m(gfmsg(3)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(gfmsg(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(gfmsg(1.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(gfmsg("1")); }) == EXIT_FAILURE);
}

TEST_CASE("gppFleafScale: mode 0 -> 1; mode 1 -> clamp(fleaf,0,1); mode 2 -> 1 for woody, fleaf otherwise") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (double f : {0.0, 0.3, 1.0, 1.7, -0.2}) { CHECK(Vegetation_Bgc::gppFleafScale(0, 0, f) == 1.0); CHECK(Vegetation_Bgc::gppFleafScale(0, 1, f) == 1.0); }
  CHECK(Vegetation_Bgc::gppFleafScale(0, 0, nan) == 1.0);
  CHECK(Vegetation_Bgc::gppFleafScale(1, 0, 0.3) == 0.3); CHECK(Vegetation_Bgc::gppFleafScale(1, 1, 0.3) == 0.3);
  CHECK(Vegetation_Bgc::gppFleafScale(1, 0, 1.7) == 1.0); CHECK(Vegetation_Bgc::gppFleafScale(1, 0, -0.2) == 0.0);
  CHECK(Vegetation_Bgc::gppFleafScale(1, 0, nan) == 0.0);
  CHECK(Vegetation_Bgc::gppFleafScale(2, 1, 0.3) == 1.0); CHECK(Vegetation_Bgc::gppFleafScale(2, 0, 0.3) == 0.3);
}

TEST_CASE("delta(): mode 0 ignores fleaf (bit-identical); mode 1 scales INGPP by fleaf; mode 2 skips woody") {
  const double g0 = GfVeg(0, 0, 1.0).ingpp();
  REQUIRE(g0 > 0.0);
  // sum of INGPP over compartments = gpp_all up to rounding (fleaf moves C between leaf and stem/root via dleafc)
  CHECK(GfVeg(0, 0, 0.4).ingpp() == doctest::Approx(g0).epsilon(1e-13));   // legacy: fleaf has no effect on GPP
  CHECK(GfVeg(0, 1, 0.4).ingpp() == doctest::Approx(g0).epsilon(1e-13));
  CHECK(GfVeg(1, 0, 1.0).ingpp() == g0);                 // fleaf 1: no change (x 1.0 is exact)
  CHECK(GfVeg(1, 0, 0.4).ingpp() == doctest::Approx(0.4 * g0).epsilon(1e-12));
  CHECK(GfVeg(1, 1, 0.4).ingpp() == doctest::Approx(0.4 * g0).epsilon(1e-12));
  CHECK(GfVeg(2, 1, 0.4).ingpp() == doctest::Approx(g0).epsilon(1e-13));                 // woody PFT untouched in mode 2
  CHECK(GfVeg(2, 0, 0.4).ingpp() == doctest::Approx(0.4 * g0).epsilon(1e-12));
  CHECK(GfVeg(1, 0, 0.0).ingpp() == doctest::Approx(0.0).epsilon(1e-12));
}

}
