// eq_nunlim_years (TEM-dev 2026-11-20): N-unlimited early eq spin-up switch, default 0 (legacy, one eq segment).
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/ModelData.h"
#include <json/json.h>
namespace {
Json::Value msn(Json::Value v) { Json::Value c; c["model_settings"]["eq_nunlim_years"] = v; return c; }
}
TEST_CASE("eq_nunlim_years: default 0 (ModelData() and ModelData(Json) without the key); values read; invalid values FATAL") {
  ModelData md; CHECK(md.eq_nunlim_years == 0);
  Json::Value c; c["model_settings"]["svp_ice_fix"] = 0; ModelData mj(c); CHECK(mj.eq_nunlim_years == 0);
  for (int k : {0, 1, 100, 200, 100000}) { CAPTURE(k); CHECK(ModelData(msn(k)).eq_nunlim_years == k); }
  CHECK(ut::exit_status_of([] { ModelData m(msn(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msn(100001)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msn(100.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msn(true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msn("100")); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(msn(100)); }) == 0);
}
TEST_CASE("eq_nunlim_segments: K 0 = one legacy segment; K > 0 = N-unlimited [0,K) + [K,n); invalid -> empty") {
  auto s0 = eq_nunlim_segments(0, 1000, false, 1);
  REQUIRE(s0.size() == 1); CHECK(s0[0].start == 0); CHECK(s0[0].end == 1000); CHECK_FALSE(s0[0].nfeed_off);
  auto s0f = eq_nunlim_segments(0, 1000, true, 5);   // K 0 ignores fire / interval (legacy path)
  REQUIRE(s0f.size() == 1); CHECK(s0f[0].end == 1000);
  for (int k : {1, 100, 200, 999}) {
    CAPTURE(k); auto s = eq_nunlim_segments(k, 1000, false, 1);
    REQUIRE(s.size() == 2); CHECK(s[0].start == 0); CHECK(s[0].end == k); CHECK(s[0].nfeed_off);
    CHECK(s[1].start == k); CHECK(s[1].end == 1000); CHECK_FALSE(s[1].nfeed_off);
    CHECK(s[0].end - s[0].start + s[1].end - s[1].start == 1000);   // year count unchanged
  }
  CHECK(eq_nunlim_segments(1000, 1000, false, 1).empty());   // K must be < eq years
  CHECK(eq_nunlim_segments(1001, 1000, false, 1).empty());
  CHECK(eq_nunlim_segments(-5, 1000, false, 1).empty());
  CHECK(eq_nunlim_segments(100, 1000, true, 1).empty());     // fire in eq: refused
  CHECK(eq_nunlim_segments(100, 1000, false, 0).empty());    // bad output interval
  CHECK(eq_nunlim_segments(150, 1000, false, 100).empty());  // K not a multiple of output_interval
  CHECK(eq_nunlim_segments(200, 1000, false, 100).size() == 2);
}
