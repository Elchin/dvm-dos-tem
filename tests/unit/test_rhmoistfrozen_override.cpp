// model_settings.rhmoistfrozen_override: per-CMT override of cmt_bgcsoil rhmoistfrozen (default absent/{} =
// parameter-file value, bit-identical). See include/RhmoistfrozenOverride.h and
// docs/WORK_NOTES_rhmoistfrozen_override.md. Refs: no published source (configuration switch).
#include "doctest.h"
#include "test_helpers.h"
#include "../../include/RhmoistfrozenOverride.h"
#include "../../include/CohortLookup.h"
#include "../../include/TEMUtilityFunctions.h"
#include <json/json.h>
#include <cstring>

namespace {
Json::Value rmsg(Json::Value v) { Json::Value c; c["model_settings"]["rhmoistfrozen_override"] = v; return c; }
Json::Value obj(std::initializer_list<std::pair<const char*, Json::Value>> kv) {
  Json::Value o(Json::objectValue); for (auto& p : kv) o[p.first] = p.second; return o;
}
bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }
}

TEST_CASE("rhmoistfrozen_override: absent, null and {} mean no override") {
  std::map<int, double> m; std::string err;
  CHECK(rmfo::parse(Json::Value(), m, err)); CHECK(m.empty());
  CHECK(rmfo::parse(Json::Value(Json::objectValue), m, err)); CHECK(m.empty());
  Json::Value c; c["model_settings"]["deep_carbon"] = 1;
  ModelData md0(c);
  CHECK(md0.rhmoistfrozen_override.empty());
  ModelData md1(rmsg(Json::Value(Json::objectValue)));
  CHECK(md1.rhmoistfrozen_override.empty());
}

TEST_CASE("rhmoistfrozen_override: valid keys and values") {
  std::map<int, double> m; std::string err;
  REQUIRE(rmfo::parse(obj({{"70", 0.0}, {"CMT66", 0.05}, {"7", 1}, {"05", 0.25}}), m, err));
  CHECK(m.size() == 4);
  CHECK(m.at(70) == 0.0); CHECK(m.at(66) == 0.05); CHECK(m.at(7) == 1.0); CHECK(m.at(5) == 0.25);
  double v = -1;
  CHECK(rmfo::lookup(m, temutil::cmtcode2num("CMT70"), v)); CHECK(v == 0.0);
  CHECK(rmfo::lookup(m, temutil::cmtcode2num("CMT05"), v)); CHECK(v == 0.25);
  CHECK_FALSE(rmfo::lookup(m, temutil::cmtcode2num("CMT73"), v));
  ModelData md(rmsg(obj({{"70", 0.0}})));
  CHECK(md.rhmoistfrozen_override.size() == 1);
  CHECK(md.rhmoistfrozen_override.at(70) == 0.0);
}

TEST_CASE("rhmoistfrozen_override: invalid input is rejected (FATAL in ModelData)") {
  std::map<int, double> m; std::string err;
  CHECK_FALSE(rmfo::parse(Json::Value(0.0), m, err));                       // not an object
  CHECK_FALSE(rmfo::parse(Json::Value("70"), m, err));
  Json::Value arr(Json::arrayValue); arr.append(0.0);
  CHECK_FALSE(rmfo::parse(arr, m, err));
  CHECK_FALSE(rmfo::parse(obj({{"x70", 0.0}}), m, err));                    // bad keys
  CHECK_FALSE(rmfo::parse(obj({{"700", 0.0}}), m, err));
  CHECK_FALSE(rmfo::parse(obj({{"", 0.0}}), m, err));
  CHECK_FALSE(rmfo::parse(obj({{"CMT", 0.0}}), m, err));
  CHECK_FALSE(rmfo::parse(obj({{"-7", 0.0}}), m, err));
  CHECK_FALSE(rmfo::parse(obj({{"70", "0.0"}}), m, err));                   // bad values
  CHECK_FALSE(rmfo::parse(obj({{"70", true}}), m, err));
  CHECK_FALSE(rmfo::parse(obj({{"70", -0.1}}), m, err));
  CHECK_FALSE(rmfo::parse(obj({{"70", 1.5}}), m, err));
  CHECK_FALSE(rmfo::parse(obj({{"70", Json::Value()}}), m, err));
  CHECK_FALSE(rmfo::parse(obj({{"70", 0.0}, {"CMT70", 0.1}}), m, err));     // duplicate CMT
  CHECK(err.find("twice") != std::string::npos);
}

TEST_CASE("rhmoistfrozen_override: Soil_Bgc::initializeParameter uses the override, else the file value") {
  static CohortLookup ch;
  ch.rhmoistfrozen = 1.0;
  Soil_Bgc sb; sb.setCohortLookup(&ch);
  sb.initializeParameter();
  CHECK(same_bits(sb.bgcpar.rhmoistfrozen, 1.0));          // default: file value, untouched
  sb.rmf_override_on = 1; sb.rmf_override_val = 0.0;
  sb.initializeParameter();
  CHECK(same_bits(sb.bgcpar.rhmoistfrozen, 0.0));
  sb.rmf_override_on = 0;
  sb.initializeParameter();
  CHECK(same_bits(sb.bgcpar.rhmoistfrozen, 1.0));
}
