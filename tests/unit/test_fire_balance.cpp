// Fire-balance fixes (Elchin 2026-11-12), model_settings.fire_n_redep_fix and
// model_settings.fire_c_bookkeep_fix, both 0 (legacy) by default. See docs/WORK_NOTES_fire_balance.md.
//   N: the fire-emitted N re-deposition (FirData.fire_a2soi.orgn, annual rate) enters top-layer ORGN
//      at 1/12 per month only with fire_n_redep_fix 1 (legacy gate on the never-set CohortData.fri).
//   C: the standing-dead (DEADC/DEADN) fall to woody debris conserves C and N with fire_c_bookkeep_fix 1.
// Pure mass-bookkeeping fixes: no published source.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value fbm(const char* k, Json::Value v) { Json::Value c; c["model_settings"][k] = v; return c; }

// Soil column for one Soil_Bgc::deltastate() call (as test_bnf_v5q BnfSoil).
struct FbSoil {
  CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
  explicit FbSoil(int n, int fix) {
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.fire_n_redep_fix = fix;
    sb.bgcpar.som2co2 = 0.2; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    cd.m_soil.numsl = n; cd.m_soil.deepthick = 0.1; sb.totdzliq = 0.0;
    cd.fri = 0; cd.yrsdist = 3;   // legacy gate value in the runs: fri is never set (here 0)
    double z = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      bool on = i < n;
      cd.m_soil.type[i] = 3; cd.m_soil.dz[i] = 0.1; cd.m_soil.z[i] = z; z += 0.1;
      sb.tmp_sois.rawc[i] = 100; sb.tmp_sois.soma[i] = 50; sb.tmp_sois.sompr[i] = 200; sb.tmp_sois.somcr[i] = 300;
      sb.tmp_sois.orgn[i] = 20; sb.tmp_sois.avln[i] = 1; bd.m_sois.avln[i] = on ? 0.5 : 0.0;
      sb.del_soi2a.rhrawc[i] = 3.0; sb.del_soi2a.rhsoma[i] = 0.5; sb.del_soi2a.rhsompr[i] = 0.1; sb.del_soi2a.rhsomcr[i] = 0.01;
      sb.ltrflc[i] = on ? 2.0 : 0.0; sb.ltrfln[i] = on ? 0.1 : 0.0; sb.rtnextract[i] = on ? 0.02 : 0.0;
      sb.del_soi2soi.netnmin[i] = on ? 0.05 : 0.0; bd.m_v2soi.rtlfalfrac[i] = on ? 0.5 : 0.0;
      ed.m_sois.liq[i] = 10.0;
    }
    sb.tmp_sois.wdebrisc = 0; sb.tmp_sois.wdebrisn = 0; sb.d2wdebrisc = 0; sb.d2wdebrisn = 0; sb.del_soi2a.rhwdeb = 0;
    sb.del_soi2a.rhwdeb = 0;
    bd.m_a2soi.orgninput = 0; bd.m_a2soi.avlninput = 0; sb.del_soi2l.avlnlost = 0; ed.m_sois.draindepth = 0;
  }
  double dN() { double s = 0; for (int i = 0; i < cd.m_soil.numsl; ++i) s += sb.del_sois.orgn[i] + sb.del_sois.avln[i]; return s; }
  double internal() { double s = 0; for (int i = 0; i < cd.m_soil.numsl; ++i) s += sb.ltrfln[i] - sb.rtnextract[i]; return s; }
};

// One PFT for Vegetation_Bgc::prepareIntegration() (as test_gppfleaf_v5q GfVeg), standing dead set.
struct FbVeg {
  CohortData cd; BgcData bd; EnvData ed; Vegetation_Bgc vb;
  FbVeg(int fix, double deadc, double deadc0, int mthsdist) {
    vb.ipft = 0; vb.setCohortData(&cd); vb.setBgcData(&bd); vb.setEnvData(&ed);
    vb.fire_c_bookkeep_fix = fix;
    ed.m_l2a.pet = 60.0; ed.m_l2a.eet = 30.0; ed.m_atms.co2 = 400.0; ed.m_atms.ta = 12.0;
    ed.m_a2l.par = 120.0; ed.m_soid.rtdpthawpct = 1.0;
    cd.m_veg.ifwoody[0] = 1; cd.m_veg.ifdeciwoody[0] = 0; cd.m_veg.ifperenial[0] = 1;
    cd.m_vegd.fleaf[0] = 1.0; cd.m_vegd.ffoliage[0] = 1.0; cd.m_vegd.maxleafc[0] = 50.0;
    cd.mthsdist = mthsdist; bd.m_vegs.deadc = deadc; bd.m_vegs.deadc0 = deadc0;
    bd.m_vegs.deadn = deadc / 50.0; bd.m_vegs.deadn0 = deadc0 / 50.0;
    bd.m_v2soi.d2wdebrisc = 0.0; bd.m_v2soi.d2wdebrisn = 0.0;
    vb.bgcpar.tmin = -1; vb.bgcpar.tmax = 30; vb.bgcpar.toptmin = 10; vb.bgcpar.toptmax = 20;
    vb.bgcpar.ki = 75; vb.bgcpar.kc = 400; vb.bgcpar.raq10a0 = 2.35665; vb.bgcpar.raq10a1 = -0.053077;
    vb.bgcpar.raq10a2 = 0.00229; vb.bgcpar.raq10a3 = -0.0000412;
    vb.calpar.cmax = 300; vb.calpar.frg = 0.2;
    for (int i = 0; i < NUM_PFT_PART; ++i) {
      vb.bgcpar.cpart[i] = (i == I_leaf) ? 0.5 : 0.25; vb.calpar.kra = -8.06e-05;
      vb.calpar.krb[i] = -6.0; bd.m_vegs.c[i] = 10.0; vb.tmp_vegs.c[i] = 10.0;
    }
  }
  // one month: standing dead -> DWD transfer (prepareIntegration, nfeed off), then the month counter
  double month() { vb.prepareIntegration(false); double f = bd.m_v2soi.d2wdebrisc; cd.mthsdist++; return f; }
};
}

TEST_SUITE("fire_balance") {

TEST_CASE("fire-balance switches: default 0 everywhere; integers 0/1 only, anything else FATAL") {
  ModelData md; CHECK(md.fire_n_redep_fix == 0); CHECK(md.fire_c_bookkeep_fix == 0);
  ModelData mj(fbm("fire", 1)); CHECK(mj.fire_n_redep_fix == 0); CHECK(mj.fire_c_bookkeep_fix == 0);
  Soil_Bgc sb; CHECK(sb.fire_n_redep_fix == 0);
  Vegetation_Bgc vb; CHECK(vb.fire_c_bookkeep_fix == 0);
  for (const char* k : {"fire_n_redep_fix", "fire_c_bookkeep_fix"}) {
    CAPTURE(k);
    CHECK(ut::exit_status_of([k] { ModelData m(fbm(k, 0)); }) == 0);
    CHECK(ut::exit_status_of([k] { ModelData m(fbm(k, 1)); }) == 0);
    CHECK(ut::exit_status_of([k] { ModelData m(fbm(k, 2)); }) == EXIT_FAILURE);
    CHECK(ut::exit_status_of([k] { ModelData m(fbm(k, -1)); }) == EXIT_FAILURE);
    CHECK(ut::exit_status_of([k] { ModelData m(fbm(k, 0.5)); }) == EXIT_FAILURE);
    CHECK(ut::exit_status_of([k] { ModelData m(fbm(k, "1")); }) == EXIT_FAILURE);
  }
  CHECK(ModelData(fbm("fire_n_redep_fix", 1)).fire_n_redep_fix == 1);
  CHECK(ModelData(fbm("fire_c_bookkeep_fix", 1)).fire_c_bookkeep_fix == 1);
}

TEST_CASE("fire-month N balance: re-deposition reaches top-layer ORGN only with fire_n_redep_fix 1") {
  const double a2s = 11.5;                          // BURNAIR2SOILN, annual rate (g N/m2/yr)
  FbSoil leg(3, 0), fix(3, 1), none(3, 1);
  leg.fd.fire_a2soi.orgn = a2s; fix.fd.fire_a2soi.orgn = a2s; none.fd.fire_a2soi.orgn = 0.0;
  leg.sb.setFireNRedepInput(); fix.sb.setFireNRedepInput(); none.sb.setFireNRedepInput();
  CHECK(leg.bd.m_a2soi.orgninput == 0.0);           // legacy: reported in BURNAIR2SOILN, never applied
  CHECK(fix.bd.m_a2soi.orgninput == a2s / 12.0);
  CHECK(none.bd.m_a2soi.orgninput == 0.0);          // after one FRI Cohort zeroes fire_a2soi.orgn
  leg.sb.deltastate(); fix.sb.deltastate();
  // the whole input lands in ORGN of layer 0; the monthly N balance closes with it as an input
  CHECK(fix.sb.del_sois.orgn[0] - leg.sb.del_sois.orgn[0] == doctest::Approx(a2s / 12.0).epsilon(1e-14));
  for (int i = 1; i < 3; ++i) CHECK(fix.sb.del_sois.orgn[i] == leg.sb.del_sois.orgn[i]);
  CHECK(fix.dN() - fix.internal() == doctest::Approx(a2s / 12.0).epsilon(1e-12));
  CHECK(leg.dN() - leg.internal() == doctest::Approx(0.0).epsilon(1e-12));    // legacy gap = a2s/12 per month
  // 12 months of the fixed path return exactly the annual rate
  double yr = 0.0; for (int m = 0; m < 12; ++m) { fix.sb.setFireNRedepInput(); yr += fix.bd.m_a2soi.orgninput; }
  CHECK(yr == doctest::Approx(a2s).epsilon(1e-14));
}

TEST_CASE("standing-dead C: fix 1 conserves DEADC + DWD input in fire years and after; legacy does not") {
  const double d0 = 72.0;
  // fall window after a fire: d0/108 per month, total over 108 months = d0, nothing after
  for (int fix : {0, 1}) {
    CAPTURE(fix);
    FbVeg v(fix, d0, d0, 0);
    double fell = 0.0, after = 0.0;
    for (int m = 0; m < 108; ++m) fell += v.month();
    for (int m = 108; m < 168; ++m) after += v.month();
    CHECK(fell == doctest::Approx(d0 * 107.0 / 108.0 + (fix ? 0.0 : 0.0) + d0 / 108.0).epsilon(1e-9));
    if (fix == 1) {
      CHECK(v.bd.m_vegs.deadc == doctest::Approx(0.0).epsilon(1e-12));
      CHECK(fell + after + v.bd.m_vegs.deadc == doctest::Approx(d0).epsilon(1e-12));   // C conserved
      CHECK(after == doctest::Approx(0.0).epsilon(1e-9));
    } else {
      // legacy: the last month's fall (d0/108) stays in m_v2soi.d2wdebrisc after the window, and
      // Soil_Bgc adds it to DWD every month -> C created (the +8.5 g/yr of the classic-fire smoke case)
      CHECK(after == doctest::Approx(60.0 * d0 / 108.0).epsilon(1e-9));
    }
  }
}

TEST_CASE("standing-dead C after a restart (deadc0 unknown): legacy deletes it at 108 months, fix 1 lets it fall") {
  const double d = 72.0;
  FbVeg leg(0, d, 0.0, 0), fix(1, d, 0.0, 0);
  double fl = 0.0, ff = 0.0;
  for (int m = 0; m < 120; ++m) { fl += leg.month(); ff += fix.month(); }
  CHECK(fl == 0.0); CHECK(leg.bd.m_vegs.deadc == 0.0);            // legacy: 72 g vanish at month 108
  CHECK(ff + fix.bd.m_vegs.deadc == doctest::Approx(d).epsilon(1e-12));
  CHECK(fix.bd.m_vegs.deadc == doctest::Approx(0.0).epsilon(1e-12));
  // remainder falls linearly over the rest of the window: equal monthly amounts d/108
  FbVeg one(1, d, 0.0, 100); CHECK(one.month() == doctest::Approx(d / 8.0).epsilon(1e-14));
  // past the window, all that is left falls at once (no deletion)
  FbVeg late(1, 5.0, 3.0, 130); CHECK(late.month() == 5.0); CHECK(late.bd.m_vegs.deadc == 0.0);
  FbVeg lateleg(0, 5.0, 3.0, 130); CHECK(lateleg.month() == 0.0); CHECK(lateleg.bd.m_vegs.deadc == 0.0);
}

TEST_CASE("DWD side: Soil_Bgc::deltastate adds the monthly fall to woody debris (C and N)") {
  FbSoil s(3, 0);
  s.sb.d2wdebrisc = 0.6; s.sb.d2wdebrisn = 0.012; s.sb.del_soi2a.rhwdeb = 0.1;
  s.sb.deltastate();
  CHECK(s.sb.del_sois.wdebrisc == doctest::Approx(0.5).epsilon(1e-14));
}

}
