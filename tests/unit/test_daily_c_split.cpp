// v5p: statistical daily split of the monthly C fluxes (output only).
//   model_settings.daily_c_split (0 = off, default; 1 = on) -> ModelData::daily_c_split
//   outputs GPP_DS NPP_DS RM_DS RG_DS RHSOM_DS RHDWD_DS RECO_DS NEE_DS INGPP_DS INNPP_DS
//   (daily only, ecosystem totals, g/m2/day) and DS_FALLBACK (monthly flags).
// X_d = M * w_d / sum(w): GPP/INGPP weight = day-varying part of getGPP, RM = raq10(Ta_d),
// RG = max(GPP_d - RM_d, 0) (then GPP weight, then uniform), RH = rhmoist*rhq10 per layer,
// NPP/INNPP/RECO/NEE derived; sum over the calendar month == monthly value.
// Refs: no published source for this exact method. Olsen & Randerson 2004,
//   doi:10.1029/2003JD003968; Fisher et al. 2016, doi:10.5194/bg-13-4271-2016;
//   exact preservation of the total (principle only): Denton 1971,
//   doi:10.1080/01621459.1971.10482227; Chow & Lin 1971, doi:10.2307/1928739.
#include <random>   // system headers before test_access.h (via test_helpers.h)
#include <algorithm>
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/DailyCSplit.h"
#include <json/json.h>
#include <cmath>

namespace {
Json::Value ms5p(const char* k, Json::Value v) { Json::Value c; c["model_settings"][k] = v; return c; }
double sum(const double* x, int n) { double s = 0; for (int d = 0; d < n; ++d) s += x[d]; return s; }
const double TOL = 1e-13;   // relative: float round-off of n <= 31 terms
}

TEST_SUITE("daily_c_split") {

TEST_CASE("switch: default 0 in ModelData() and config; 0/1 parse; anything else is FATAL") {
  ModelData md; CHECK(md.daily_c_split == 0);
  CHECK(ModelData(ms5p("svp_ice_fix", 0)).daily_c_split == 0);   // key absent -> 0
  CHECK(ModelData(ms5p("daily_c_split", 0)).daily_c_split == 0);
  CHECK(ModelData(ms5p("daily_c_split", 1)).daily_c_split == 1);
  CHECK(ut::exit_status_of([] { ModelData m(ms5p("daily_c_split", 2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms5p("daily_c_split", -1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms5p("daily_c_split", 0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms5p("daily_c_split", true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms5p("daily_c_split", "1")); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(ms5p("daily_c_split", 1)); }) == 0);
  CHECK(dcsplit::valid_switch(Json::Value(0)));  CHECK(dcsplit::valid_switch(Json::Value(1.0)));
  CHECK_FALSE(dcsplit::valid_switch(Json::Value(1.5)));
  CHECK_FALSE(dcsplit::valid_switch(Json::Value(false)));
  CHECK_FALSE(dcsplit::valid_switch(Json::Value("0")));
}

TEST_CASE("output_spec requests: *_DS daily only, DS_FALLBACK monthly only, any of them with the switch at 0 is FATAL") {
  using dcsplit::request_error;
  for (int k = 0; k < dcsplit::NVAR; ++k) {
    std::string nm = dcsplit::VAR_NAMES[k]; CAPTURE(nm);
    CHECK(dcsplit::is_ds_name(nm));
    CHECK(request_error(nm, true, true, 0) != "");     // requested while off -> FATAL
    CHECK(request_error(nm, true, true, 1) == "");     // daily, switch on -> ok
    CHECK(request_error(nm, true, false, 1) != "");    // monthly/yearly only -> FATAL
    CHECK(request_error(nm, false, false, 0) == "");   // not requested -> ignored
  }
  CHECK(request_error("DS_FALLBACK", true, false, 0) != "");
  CHECK(request_error("DS_FALLBACK", true, false, 1) == "");
  CHECK(request_error("DS_FALLBACK", true, true, 1) != "");
  CHECK_FALSE(dcsplit::is_ds_name("GPP")); CHECK(request_error("GPP", true, true, 0) == "");
  CHECK(std::string(dcsplit::VAR_NAMES[dcsplit::NEE]) == "NEE_DS");
}

TEST_CASE("split: proportional to the weights, sum exact, negative/NaN weights count as 0") {
  const double w[5] = {1, 2, 3, 0, 4};
  double x[5];
  CHECK_FALSE(dcsplit::split(10.0, w, 5, x));
  CHECK(x[0] == doctest::Approx(1.0)); CHECK(x[1] == doctest::Approx(2.0)); CHECK(x[2] == doctest::Approx(3.0));
  CHECK(x[3] == 0.0); CHECK(x[4] == doctest::Approx(4.0));
  CHECK(std::fabs(sum(x, 5) - 10.0) <= TOL * 10.0);
  const double wb[4] = {-1.0, NAN, INFINITY, 2.0};
  CHECK_FALSE(dcsplit::split(3.0, wb, 4, x));
  CHECK(x[0] == 0.0); CHECK(x[1] == 0.0); CHECK(x[2] == 0.0); CHECK(x[3] == 3.0);
  CHECK(dcsplit::split(-2.0, w, 5, x) == false);          // negative totals split the same way
  CHECK(std::fabs(sum(x, 5) + 2.0) <= TOL * 2.0); CHECK(x[4] == doctest::Approx(-0.8));
}

TEST_CASE("zero flux: M == 0 gives all-zero days whatever the weights (no fallback)") {
  const double w[3] = {1, NAN, 0};
  double x[3] = {9, 9, 9};
  CHECK_FALSE(dcsplit::split(0.0, w, 3, x));
  CHECK(x[0] == 0.0); CHECK(x[1] == 0.0); CHECK(x[2] == 0.0);
  CHECK(dcsplit::split2(0.0, w, w, 3, x) == 0);
}

TEST_CASE("all-zero weights with M != 0: uniform M/n, fallback reported (split true, split2 levels)") {
  const double z[30] = {0}; const double g[30] = {0, 0, 1};
  double x[30];
  CHECK(dcsplit::split(6.0, z, 30, x));
  for (int d = 0; d < 30; ++d) CHECK(x[d] == doctest::Approx(0.2));
  CHECK(std::fabs(sum(x, 30) - 6.0) <= TOL * 6.0);
  CHECK(dcsplit::split2(6.0, g, z, 30, x) == 0);           // primary used
  CHECK(dcsplit::split2(6.0, z, g, 30, x) == 1);           // secondary (GPP weight)
  CHECK(x[2] == 6.0);
  CHECK(dcsplit::split2(6.0, z, z, 30, x) == 2);           // uniform
  CHECK(x[29] == doctest::Approx(0.2));
}

TEST_CASE("fallback counter: any fallback flag counts a month; the layer remap is counted separately") {
  long fb = 0, rm = 0;
  dcsplit::count_month(0, fb, rm);                                  CHECK(fb == 0); CHECK(rm == 0);
  dcsplit::count_month(dcsplit::FLAG_LAYER_REMAP, fb, rm);          CHECK(fb == 0); CHECK(rm == 1);
  dcsplit::count_month(dcsplit::FLAG_RG_GPP, fb, rm);               CHECK(fb == 1); CHECK(rm == 1);
  dcsplit::count_month(dcsplit::FLAG_GPP | dcsplit::FLAG_RG_UNI | dcsplit::FLAG_LAYER_REMAP, fb, rm);
  CHECK(fb == 2); CHECK(rm == 2);
  const int all[7] = {dcsplit::FLAG_GPP, dcsplit::FLAG_RM, dcsplit::FLAG_RG_GPP, dcsplit::FLAG_RG_UNI,
                      dcsplit::FLAG_RHSOM, dcsplit::FLAG_RHDWD, dcsplit::FLAG_LAYER_REMAP};
  int o = 0; for (int i = 0; i < 7; ++i) { CHECK((o & all[i]) == 0); o |= all[i]; }   // distinct bits
}

TEST_CASE("RG weight max(GPP_d - RM_d, 0); add_pft sums, derived NPP/INNPP, flags") {
  CHECK(dcsplit::rg_weight(3.0, 1.0) == 2.0); CHECK(dcsplit::rg_weight(1.0, 3.0) == 0.0);
  const int n = 4;
  const double wg[n] = {0, 1, 2, 1}, wr[n] = {1, 1, 1, 1};
  dcsplit::PftMonth m{8.0, 2.0, 1.0, 10.0, 7.0};             // gpp rm rg ingpp innpp
  double acc[dcsplit::NVAR][31]; dcsplit::clear(acc);
  int fl = dcsplit::add_pft(m, wg, wr, n, acc);
  CHECK(fl == 0);
  // GPP_d = 0,2,4,2; RM_d = 0.5 each; RG weight = 0,1.5,3.5,1.5 -> RG_d = 0,0.2307..,0.538..,0.2307..
  CHECK(acc[dcsplit::GPP][0] == 0.0); CHECK(acc[dcsplit::GPP][2] == doctest::Approx(4.0));
  CHECK(acc[dcsplit::RG][0] == 0.0); CHECK(acc[dcsplit::RG][2] == doctest::Approx(3.5 / 6.5));
  CHECK(std::fabs(sum(acc[dcsplit::GPP], n) - 8.0) <= TOL * 8);
  CHECK(std::fabs(sum(acc[dcsplit::RM], n) - 2.0) <= TOL * 2);
  CHECK(std::fabs(sum(acc[dcsplit::RG], n) - 1.0) <= TOL * 1);
  CHECK(std::fabs(sum(acc[dcsplit::INGPP], n) - 10.0) <= TOL * 10);
  CHECK(std::fabs(sum(acc[dcsplit::INNPP], n) - 7.0) <= TOL * 10);
  CHECK(std::fabs(sum(acc[dcsplit::NPP], n) - (8.0 - 2.0 - 1.0)) <= TOL * 10);
  for (int d = 0; d < n; ++d)
    CHECK(acc[dcsplit::NPP][d] == doctest::Approx(acc[dcsplit::GPP][d] - acc[dcsplit::RM][d] - acc[dcsplit::RG][d]));
  // winter-like month: GPP_d < RM_d every day -> RG on the GPP weight (flag 4)
  dcsplit::clear(acc);
  dcsplit::PftMonth w{0.4, 2.0, 0.1, 0.4, -1.7};
  fl = dcsplit::add_pft(w, wg, wr, n, acc);
  CHECK((fl & dcsplit::FLAG_RG_GPP) != 0); CHECK((fl & dcsplit::FLAG_RG_UNI) == 0);
  CHECK(acc[dcsplit::RG][0] == 0.0); CHECK(std::fabs(sum(acc[dcsplit::RG], n) - 0.1) <= TOL);
  // no GPP weight at all but GPP != 0 -> GPP and RG uniform (flags 1, 8)
  dcsplit::clear(acc);
  const double z[n] = {0, 0, 0, 0};
  fl = dcsplit::add_pft(w, z, wr, n, acc);
  CHECK((fl & dcsplit::FLAG_GPP) != 0); CHECK((fl & dcsplit::FLAG_RG_UNI) != 0);
  CHECK(acc[dcsplit::GPP][3] == doctest::Approx(0.1));
  // zero-flux PFT adds nothing
  dcsplit::clear(acc);
  dcsplit::PftMonth zz{0, 0, 0, 0, 0};
  CHECK(dcsplit::add_pft(zz, z, z, n, acc) == 0);
  for (int k = 0; k < dcsplit::NVAR; ++k) CHECK(sum(acc[k], n) == 0.0);
}

TEST_CASE("RH per layer, RECO = RHSOM + RHDWD + RM + RG, NEE = RECO - GPP") {
  const int n = 3;
  double acc[dcsplit::NVAR][31]; dcsplit::clear(acc);
  const double w1[n] = {1, 1, 2}, w0[n] = {0, 0, 0};
  CHECK_FALSE(dcsplit::add_rh(4.0, w1, n, acc[dcsplit::RHSOM]));
  CHECK(dcsplit::add_rh(3.0, w0, n, acc[dcsplit::RHSOM]));            // uniform fallback
  CHECK(acc[dcsplit::RHSOM][2] == doctest::Approx(2.0 + 1.0));
  CHECK(std::fabs(sum(acc[dcsplit::RHSOM], n) - 7.0) <= TOL * 7);
  acc[dcsplit::RHDWD][1] = 0.5; acc[dcsplit::RM][0] = 0.25; acc[dcsplit::RG][2] = 0.125; acc[dcsplit::GPP][2] = 4.0;
  dcsplit::finish(acc, n);
  CHECK(acc[dcsplit::RECO][2] == doctest::Approx(3.0 + 0.125));
  CHECK(acc[dcsplit::NEE][2] == doctest::Approx(3.125 - 4.0));
  CHECK(acc[dcsplit::NEE][0] == doctest::Approx((1.0 + 1.0) + 0.25));  // RHSOM_0 = 4/4 + 3/3
}

TEST_CASE("property: 500 random months (n = 28..31, mixed zero/negative/NaN weights) sum to the monthly value") {
  std::mt19937 rng(20261028); std::uniform_real_distribution<double> U(0.0, 1.0);
  double worst = 0;
  for (int t = 0; t < 500; ++t) {
    int n = 28 + t % 4; double w[31], wr[31], acc[dcsplit::NVAR][31]; dcsplit::clear(acc);
    for (int d = 0; d < n; ++d) {
      double u = U(rng);
      w[d] = u < 0.3 ? 0.0 : (u < 0.35 ? -u : (u < 0.37 ? NAN : U(rng) * 1e3));
      wr[d] = 0.5 + U(rng);
    }
    dcsplit::PftMonth m{U(rng) * 300, U(rng) * 80, U(rng) * 40, 0, 0};
    m.ingpp = m.gpp * (1 + U(rng)); m.innpp = m.ingpp - m.rm - U(rng) * 50;
    dcsplit::add_pft(m, w, wr, n, acc);
    double rh = U(rng) * 50; dcsplit::add_rh(rh, w, n, acc[dcsplit::RHSOM]);
    dcsplit::finish(acc, n);
    const double M[6] = {m.gpp, m.rm, m.rg, m.ingpp, m.innpp, m.gpp - m.rm - m.rg};
    const int K[6] = {dcsplit::GPP, dcsplit::RM, dcsplit::RG, dcsplit::INGPP, dcsplit::INNPP, dcsplit::NPP};
    for (int i = 0; i < 6; ++i) {
      double e = std::fabs(sum(acc[K[i]], n) - M[i]) / std::max(1.0, std::fabs(M[i]));
      worst = std::max(worst, e);
    }
    worst = std::max(worst, std::fabs(sum(acc[dcsplit::NEE], n) - (rh + m.rm + m.rg - m.gpp)) / 300.0);
    for (int k = 0; k < dcsplit::NVAR; ++k) for (int d = 0; d < n; ++d) REQUIRE(std::isfinite(acc[k][d]));
  }
  CHECK(worst <= 1e-12);
}

TEST_CASE("model weights use the model's own scalars (Vegetation_Bgc getGV/getTempFactor4GPP/getRaq10, Soil_Bgc getRhmoist/getRhq10)") {
  CohortData cd; BgcData bd; EnvData ed; Vegetation_Bgc vb;
  vb.setCohortData(&cd); vb.setBgcData(&bd); vb.setEnvData(&ed);
  vb.bgcpar.kc = 400; vb.bgcpar.ki = 75; vb.bgcpar.tmin = -1; vb.bgcpar.tmax = 30; vb.bgcpar.toptmax = 20;
  vb.bgcpar.raq10a0 = 2.35; vb.bgcpar.raq10a1 = -0.03; vb.bgcpar.raq10a2 = 0; vb.bgcpar.raq10a3 = 0;
  vb.prvtopt_used = 12; ed.m_atms.co2 = 400;
  const double ta = 8, par = 120, eet = 2, pet = 3;
  double gv = vb.getGV(eet, pet), ft = vb.getTempFactor4GPP(ta, 12.0);
  double ci = 400 * gv, expect = ci / (400 + ci) * (par / (75 + par)) * ft;
  CHECK(expect > 0);
  CHECK(vb.ds_gpp_weight(ta, par, eet, pet, 1.0) == doctest::Approx(expect).epsilon(1e-14));
  CHECK(vb.ds_gpp_weight(ta, par, eet, pet, 0.0) == 0.0);      // frozen day
  CHECK(vb.ds_gpp_weight(ta, 0.0, eet, pet, 1.0) == 0.0);      // dark day
  CHECK(vb.ds_gpp_weight(-5, par, eet, pet, 1.0) == 0.0);      // below tmin
  CHECK(vb.ds_gpp_weight(ta, par, eet, 0.0, 1.0) == doctest::Approx(  // PET floor as in prepareIntegration
        dcsplit::gpp_weight(400, vb.getGV(eet, 0.0001), 400, par, 75, ft, 1.0)));
  CHECK(vb.ds_rm_weight(ta) == doctest::Approx(vb.getRaq10(ta)).epsilon(1e-14));
  Soil_Bgc sb; sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed);
  sb.bgcpar.rhq10 = 2.0; sb.bgcpar.moistmin = 0; sb.bgcpar.moistmax = 1; sb.bgcpar.moistopt = 0.5; sb.bgcpar.rhmoistfrozen = 0.1;
  CHECK(sb.ds_rh_weight(5.0, 0.3) == doctest::Approx(sb.getRhmoist(0.3, 0, 1, 0.5) * sb.getRhq10(5.0)).epsilon(1e-14));
  CHECK(sb.ds_rh_weight(-3.0, 0.3) == doctest::Approx(0.1 * sb.getRhq10(-3.0)).epsilon(1e-14));
  CHECK(sb.ds_rh_weight(5.0, 1.2) == 0.0);                     // rhmoist clipped at 0 -> weight 0
}

TEST_CASE("layer helpers: midpoint remap and structure comparison") {
  const double z[3] = {0, 0.1, 0.3}, dz[3] = {0.1, 0.2, 0.5};
  CHECK(dcsplit::layer_containing(0.05, z, dz, 3) == 0);
  CHECK(dcsplit::layer_containing(0.1, z, dz, 3) == 1);
  CHECK(dcsplit::layer_containing(0.79, z, dz, 3) == 2);
  CHECK(dcsplit::layer_containing(5.0, z, dz, 3) == 2);
  CHECK(dcsplit::layer_containing(0.0, z, dz, 0) == -1);
  CHECK(dcsplit::same_layers(3, z, dz, 3, z, dz));
  const double dz2[3] = {0.1, 0.2, 0.4};
  CHECK_FALSE(dcsplit::same_layers(3, z, dz, 3, z, dz2));
  CHECK_FALSE(dcsplit::same_layers(3, z, dz, 2, z, dz));
}

}
