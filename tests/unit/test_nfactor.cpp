// Area 3: n-factor selection (nfactor_sameday on/off, first day of the month).
// Uses the behaviour-neutral seam include/NFactorSelect.h (commit "seam: n-factor ...").
// The day loop below mirrors Cohort::updateMonthly_Env: tdrv starts each month at
// TDRV_MONTH_START, the n-factor is selected, then tdrv = today's Tair.
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/NFactorSelect.h"

using namespace nfactor_sel;

namespace {
const double NFS = 1.15, NFW = 0.45;
std::vector<double> month_nfactors(int sameday, const std::vector<double>& tair) {
  std::vector<double> nf;
  double tdrv = TDRV_MONTH_START;
  for (size_t d = 0; d < tair.size(); ++d) {
    double tsel = selector_tair(sameday, tair[d], tdrv);
    nf.push_back(select_nfactor(tsel, NFS, NFW));
    tdrv = tair[d];
  }
  return nf;
}
}

TEST_SUITE("nfactor") {

TEST_CASE("select_nfactor: winter n-factor at and below 0 C, summer above") {
  CHECK(select_nfactor(-10.0, NFS, NFW) == NFW);
  CHECK(select_nfactor(0.0, NFS, NFW) == NFW);             // freezing includes 0 C
  CHECK(select_nfactor(1e-12, NFS, NFW) == NFS);
  CHECK(select_nfactor(6.95e-310, NFS, NFW) == NFS);       // the v5h stale-stack denormal -> summer
  CHECK(select_nfactor(25.0, NFS, NFW) == NFS);
  CHECK(select_nfactor(-0.0, NFS, NFW) == NFW);
}

TEST_CASE("selector_tair: sameday=1 uses today's Tair, otherwise the previous day's") {
  CHECK(selector_tair(1, -5.0, 3.0) == -5.0);
  CHECK(selector_tair(0, -5.0, 3.0) == 3.0);
  CHECK(selector_tair(2, -5.0, 3.0) == 3.0);               // only 1 switches it on
  CHECK(TDRV_MONTH_START > 0.0);                           // only the sign matters
}

TEST_CASE("legacy (sameday=0): one-day lag, day 1 of every month uses the SUMMER n-factor") {
  std::vector<double> t = {-20, -18, 5, -3, 2};
  std::vector<double> nf = month_nfactors(0, t);
  CHECK(nf[0] == NFS);   // first day of month: tdrv = TDRV_MONTH_START (> 0) even at -20 C
  CHECK(nf[1] == NFW);   // previous day -20
  CHECK(nf[2] == NFW);   // previous day -18 although today is +5
  CHECK(nf[3] == NFS);   // previous day +5 although today is -3
  CHECK(nf[4] == NFW);
}

TEST_CASE("sameday=1: n-factor follows today's Tair, including day 1 of the month") {
  std::vector<double> t = {-20, -18, 5, -3, 2};
  std::vector<double> nf = month_nfactors(1, t);
  CHECK(nf[0] == NFW); CHECK(nf[1] == NFW); CHECK(nf[2] == NFS); CHECK(nf[3] == NFW); CHECK(nf[4] == NFS);
}

TEST_CASE("constant-sign months: the two modes differ only on day 1 of a cold month") {
  std::vector<double> cold(31, -15.0), warm(30, 12.0);
  std::vector<double> c0 = month_nfactors(0, cold), c1 = month_nfactors(1, cold);
  std::vector<double> w0 = month_nfactors(0, warm), w1 = month_nfactors(1, warm);
  CHECK(c0[0] == NFS); CHECK(c1[0] == NFW);
  for (size_t d = 1; d < cold.size(); ++d) CHECK(c0[d] == c1[d]);
  for (size_t d = 0; d < warm.size(); ++d) CHECK(w0[d] == w1[d]);
}

} // TEST_SUITE nfactor
