// v5p-daily-water-fix (output only): daily TRANSPIRATION / EET / PET NetCDF values.
// Expected values come from the model's own bookkeeping: the daily output of an accumulated
// flux must be the term EnvData adds into the monthly output, so the days sum to the month.
// Refs: no published source; rationale: output bookkeeping (see include/DailyWaterOut.h).
#include "doctest.h"
#include "test_helpers.h"
// Cohort is not in test_access.h; open its private helpers (getEd4land_daily etc.) the same way.
#define private public
#define protected public
#include "../../include/Cohort.h"
#undef private
#undef protected
#include "DailyWaterOut.h"

TEST_SUITE("daily_water_out") {

TEST_CASE("record_tran + veg_endOfDay: daily TRANSPIRATION days sum exactly to m_v2a.tran") {
  EnvData ed;
  ed.veg_beginOfYear();
  for (int m = 0; m < 12; ++m) {
    const int dinm = DINM[m];
    ed.veg_beginOfMonth();
    for (int d = 0; d < dinm; ++d) {
      ed.d_v2a.tran = (d % 5 == 0) ? 0.0 : 0.013 * (d + 1) + 0.37 * m;  // mm/day, incl. zero days
      ed.veg_endOfDay(dinm);
      dwout::record_tran(ed, d);
    }
    std::vector<double> out;
    dwout::append_days(ed.daily_tran, dinm, out);
    REQUIRE(out.size() == (size_t)dinm);
    double s = 0.0;
    for (double v : out) s += v;                  // same order as veg_endOfDay -> bit-equal
    CHECK(s == ed.m_v2a.tran);
    CHECK(out[dinm - 1] == ed.d_v2a.tran);
  }
}

TEST_CASE("record_tran ignores out-of-range day indices; daily_tran starts at 0") {
  EnvData ed;
  for (int d = 0; d < 31; ++d) CHECK(ed.daily_tran[d] == 0.0);
  ed.d_v2a.tran = 5.0;
  dwout::record_tran(ed, -1);
  dwout::record_tran(ed, 31);
  for (int d = 0; d < 31; ++d) CHECK(ed.daily_tran[d] == 0.0);
  dwout::record_tran(ed, 30);
  CHECK(ed.daily_tran[30] == 5.0);
}

TEST_CASE("append_days writes dinm values in day order and clamps dinm to [0,31]") {
  double a[31];
  for (int i = 0; i < 31; ++i) a[i] = i + 0.5;
  std::vector<double> out{-1.0};
  dwout::append_days(a, 28, out);
  REQUIRE(out.size() == 29u);
  CHECK(out[0] == -1.0);
  CHECK(out[1] == 0.5);
  CHECK(out[28] == 27.5);
  std::vector<double> o2; dwout::append_days(a, 40, o2); CHECK(o2.size() == 31u);
  std::vector<double> o3; dwout::append_days(a, -3, o3); CHECK(o3.empty());
}

TEST_CASE("cohort daily EET/PET (atm_endOfDay daily_eet/daily_pet) sum exactly to m_l2a.eet/pet") {
  EnvData ed;
  ed.atm_beginOfYear();
  ed.atm_beginOfMonth();
  const int dinm = 31;
  for (int d = 0; d < dinm; ++d) {
    ed.d_l2a.eet = 0.1 + 0.07 * d;
    ed.d_l2a.pet = 1.3 + 0.11 * d;
    ed.atm_endOfDay(dinm, d);
  }
  std::vector<double> e, p;
  dwout::append_days(ed.daily_eet, dinm, e);
  dwout::append_days(ed.daily_pet, dinm, p);
  double se = 0.0, sp = 0.0;
  for (int d = 0; d < dinm; ++d) { se += e[d]; sp += p[d]; }
  CHECK(se == ed.m_l2a.eet);
  CHECK(sp == ed.m_l2a.pet);
}

TEST_CASE("Cohort::getEd4land_daily: per-PFT EET/PET each hold the full ground terms, so their "
          "sum over PFTs (the old daily output) over-counts; the cohort value counts them once") {
  ModelData md;
  md.nv_rain_surface = 0;
  Cohort* c = new Cohort();                     // default constructor: no input files
  c->md = &md;
  EnvData all;
  c->edall = &all;
  for (int ip = 0; ip < NUM_PFT; ++ip) { c->cd.d_veg.vegcov[ip] = 0.0; c->ed[ip] = EnvData(); }
  const double veg_evap[3] = {0.20, 0.10, 0.05}, veg_tran[3] = {0.50, 0.30, 0.10};
  const double veg_evp_pet[3] = {0.40, 0.20, 0.10}, veg_tran_pet[3] = {1.0, 0.6, 0.2};
  // struct constructors fill UIN_D / NaN sentinels, so zero every term used explicitly
  auto zero_v2a = [](veg2atm_env& v) { v.tran = v.evap = v.tran_pet = v.evap_pet = v.sublim = v.evap_moss = 0.0; };
  zero_v2a(all.d_v2a);
  all.d_snw2a.sublim = 0.0; all.d_soi2a.evap = 0.0; all.d_soi2a.evap_pet = 0.0;
  for (int ip = 0; ip < NUM_PFT; ++ip) zero_v2a(c->ed[ip].d_v2a);
  for (int ip = 0; ip < 3; ++ip) {
    c->cd.d_veg.vegcov[ip] = 0.3;
    c->ed[ip].d_v2a.evap = veg_evap[ip];     c->ed[ip].d_v2a.tran = veg_tran[ip];
    c->ed[ip].d_v2a.evap_pet = veg_evp_pet[ip]; c->ed[ip].d_v2a.tran_pet = veg_tran_pet[ip];
    c->ed[ip].d_v2a.sublim = 0.0;
    all.d_v2a.evap += veg_evap[ip];  all.d_v2a.tran += veg_tran[ip];
    all.d_v2a.evap_pet += veg_evp_pet[ip]; all.d_v2a.tran_pet += veg_tran_pet[ip];
  }
  all.d_soi2a.evap = 0.8; all.d_soi2a.evap_pet = 2.0; all.d_snw2a.sublim = 0.05;   // ground terms
  c->assignGroundEd2pfts_daily();
  c->getEd4land_daily();
  const double ground_e = 0.8 + 0.05, ground_p = 2.0 + 0.05;
  const double veg_e = 0.35 + 0.90, veg_p = 0.70 + 1.80;
  CHECK(all.d_l2a.eet == doctest::Approx(veg_e + ground_e).epsilon(1e-14));
  CHECK(all.d_l2a.pet == doctest::Approx(veg_p + ground_p).epsilon(1e-14));
  double old_e = 0.0, old_p = 0.0;
  for (int ip = 0; ip < 3; ++ip) { old_e += c->ed[ip].d_l2a.eet; old_p += c->ed[ip].d_l2a.pet; }
  CHECK(old_e == doctest::Approx(veg_e + 3.0 * ground_e).epsilon(1e-14));  // ground x 3 PFTs
  CHECK(old_p == doctest::Approx(veg_p + 3.0 * ground_p).epsilon(1e-14));
  c->edall = nullptr; // 'all' is on the stack
  // Cohort is intentionally not deleted (its members were never initialised from input).
}

} // TEST_SUITE
