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
