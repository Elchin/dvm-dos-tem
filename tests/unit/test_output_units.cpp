// v5p-transp-units: output_spec.csv units labels of the accumulated water fluxes.
// TRANSPIRATION was labelled "mm/day" at every resolution, but EnvData::veg_endOfDay sums the
// daily d_v2a.tran (mm/day) into m_v2a.tran (mm/month) and EnvData::veg_endOfMonth sums that
// into y_v2a.tran (mm/year); Runner writes d_/m_/y_v2a.tran for daily/monthly/yearly output.
// The label is now "mm/time" (ModelData replaces "time" with day/month/year per file), like
// EVAP_CANOPY/EVAP_SOIL/EVAP_MOSS/SUBL_* (v5l). Metadata only, no numeric change.
// Refs: no published source; rationale: units attribute must describe the stored value
//   (CF Conventions v1.8 sec. 3.1 "Units" [CF18],
//   https://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html).
#include "doctest.h"
#include "test_helpers.h"

namespace {
// units column (3rd field) of a row in the repository's config/output_spec.csv
std::string spec_units(const std::string& var) {
  std::ifstream f("../../config/output_spec.csv");
  std::string line;
  while (std::getline(f, line)) {
    std::stringstream ss(line); std::string name, desc, units;
    std::getline(ss, name, ','); std::getline(ss, desc, ','); std::getline(ss, units, ',');
    if (name == var) return units;
  }
  return "<missing>";
}
}

TEST_SUITE("output_units") {

TEST_CASE("accumulated water fluxes carry a per-'time' unit in config/output_spec.csv") {
  const char* fluxes[] = {"TRANSPIRATION", "EVAP_CANOPY", "EVAP_MOSS", "EVAP_SOIL",
                          "SUBL_CANOPY", "SUBL_SNOW", "RAIN_INTERCEPT", "SNOW_INTERCEPT"};
  for (const char* v : fluxes) {
    CAPTURE(v);
    std::string u = spec_units(v);
    CHECK(u == "mm/time");
  }
  // EET/PET keep their historic "mm/m2/time" (time part correct; see tests/README.md)
  CHECK(spec_units("EET") == "mm/m2/time");
  CHECK(spec_units("PET") == "mm/m2/time");
}

TEST_CASE("transpiration is accumulated day -> month -> year (so the unit is per output step)") {
  EnvData ed;
  ed.veg_beginOfYear();
  double ysum = 0.0;
  for (int m = 0; m < 12; ++m) {
    const int dinm = DINM[m];
    ed.veg_beginOfMonth();
    double msum = 0.0;
    for (int d = 0; d < dinm; ++d) {
      ed.d_v2a.tran = 0.1 + 0.01 * d + 0.2 * m;   // mm/day
      msum += ed.d_v2a.tran;
      ed.veg_endOfDay(dinm);
    }
    CHECK(ed.m_v2a.tran == doctest::Approx(msum).epsilon(1e-12));                    // mm/month
    CHECK(ed.m_v2a.tran > 10.0 * ed.d_v2a.tran);
    ysum += msum;
    ed.veg_endOfMonth(m);
  }
  CHECK(ed.y_v2a.tran == doctest::Approx(ysum).epsilon(1e-12));                       // mm/year
}

} // TEST_SUITE output_units
