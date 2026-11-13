// Reburn fire master switch, model_settings.fire (0 = off, default; 1 = on). Merge plan
// (dvmdostem_merge_plan.md section 2.3) items 0-3: parse/default/FATAL, the stage-flag truth
// table (stage_settings.X.fire, falling back to X.dsb), module_settings.dsb.fire .get() defaults
// (fire_process_wildfire default false), and the 18/21-line cmt_firepar.txt parser.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cstdio>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

namespace {
Json::Value fcfg(Json::Value fire, Json::Value eq_dsb, Json::Value eq_fire) {
  Json::Value c;
  if (!fire.isNull()) c["model_settings"]["fire"] = fire;
  if (!eq_dsb.isNull()) c["stage_settings"]["eq"]["dsb"] = eq_dsb;
  if (!eq_fire.isNull()) c["stage_settings"]["eq"]["fire"] = eq_fire;
  return c;
}
const Json::Value N;  // null: key absent

std::string firepar_block(int cmt, bool reburn) {
  std::string s = "// CMT" + std::string(cmt < 10 ? "0" : "") + std::to_string(cmt) + " // test block\n";
  for (int i = 0; i < 10; ++i) s += "0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 // pft row\n";   // fvcombust, fvslash
  for (int i = 0; i < 5; ++i) s += "0.2 // foslburn\n";
  s += "0.5 // vsmburn\n0.1 // r_retain_c\n0.3 // r_retain_n\n";
  if (reburn) s += "123 // fuel model\n0.15 // cbd\n2.5 // cbh\n";
  return s;
}
}

TEST_SUITE("fire_switch") {

TEST_CASE("fire master switch: default off; integers 0/1 only") {
  ModelData md; CHECK(md.fire == 0); CHECK_FALSE(md.fire_on_EQ); CHECK_FALSE(md.fire_process_wildfire);
  ModelData mj(fcfg(N, N, N)); CHECK(mj.fire == 0);
  CHECK_FALSE(mj.fire_on_EQ); CHECK_FALSE(mj.fire_on_SP); CHECK_FALSE(mj.fire_on_TR); CHECK_FALSE(mj.fire_on_SC);
  CHECK_FALSE(mj.fire_process_wildfire);       // fork default (Reburn upstream default was true)
  CHECK(mj.fire_ignition_tr == 0); CHECK(mj.fire_ignition_sc == 0);
  CHECK(mj.fire_temp_fm == -1); CHECK(mj.fire_dwd_combust_frac == 1.0); CHECK(mj.fire_windspeed == -1.0);
  CHECK(ModelData(fcfg(1, N, N)).fire == 1);
  CHECK(ut::exit_status_of([] { ModelData m(fcfg(2, N, N)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(fcfg(-1, N, N)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(fcfg(0.5, N, N)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(fcfg(true, N, N)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(fcfg("1", N, N)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(fcfg(1, N, N)); }) == 0);
}

TEST_CASE("stage truth table: fire_on = (fire == 1) && (stage.fire if present, else stage.dsb)") {
  // {fire, eq.dsb, eq.fire, expected}
  struct Row { Json::Value f, dsb, fr; bool on; };
  const Row rows[] = {
    {0, N, N, false},     {0, true, N, false},   {0, N, true, false},   {0, true, true, false},
    {1, N, N, false},     {1, false, N, false},  {1, true, N, true},     // legacy dsb key honoured
    {1, N, false, false}, {1, N, true, true},    {1, true, false, false}, {1, false, true, true},
    {N, true, true, false},                                                // master key absent = 0
  };
  for (const Row& r : rows) {
    ModelData m(fcfg(r.f, r.dsb, r.fr));
    CAPTURE(r.f.toStyledString()); CAPTURE(r.dsb.toStyledString()); CAPTURE(r.fr.toStyledString());
    CHECK(m.fire_on_EQ == r.on);
    CHECK(m.eq_dsb == r.on);                    // TEM.cpp set_dsbmodule(eq_dsb) uses the gated value
    CHECK_FALSE(m.fire_on_SP); CHECK_FALSE(m.sp_dsb);
  }
  // pre-run: dsb is forced off when fire == 0
  Json::Value c; c["stage_settings"]["pr"]["dsb"] = true;
  CHECK_FALSE(ModelData(c).pr_dsb);
  c["model_settings"]["fire"] = 1; CHECK(ModelData(c).pr_dsb);
}

TEST_CASE("module_settings.dsb.fire is read with defaults; ignition modes 0/1 only") {
  Json::Value c = fcfg(1, N, N);
  c["module_settings"]["dsb"]["fire"]["fire_process_wildfire"] = true;
  c["module_settings"]["dsb"]["fire"]["fire_windspeed"] = 120.0;
  ModelData m(c); CHECK(m.fire_process_wildfire); CHECK(m.fire_windspeed == 120.0); CHECK(m.fire_max_snow == 1.0);
  Json::Value b = fcfg(1, N, N); b["module_settings"]["dsb"]["fire"]["ignition_sc"] = 2;
  CHECK(ut::exit_status_of([b] { ModelData m2(b); }) == EXIT_FAILURE);
}

TEST_CASE("cmt_firepar.txt: legacy 18 lines (no fuel model) or Reburn 21 lines") {
  char tmpl[] = "/tmp/ut_firepar_XXXXXX"; std::string dir = std::string(mkdtemp(tmpl)) + "/";
  { std::ofstream f(dir + "cmt_firepar.txt"); f << firepar_block(5, false) << firepar_block(7, true); }
  CohortLookup a; a.cmtcode = "CMT05"; a.assignFirePar(dir);
  CHECK(a.cmt2fm == -1); CHECK(a.cbd == 0.0); CHECK(a.cbh == 0.0); CHECK(a.r_retain_n == doctest::Approx(0.3));
  CohortLookup b; b.cmtcode = "CMT07"; b.assignFirePar(dir);
  CHECK(b.cmt2fm == 123); CHECK(b.cbd == doctest::Approx(0.15)); CHECK(b.cbh == doctest::Approx(2.5));
  CHECK(b.r_retain_n == doctest::Approx(0.3));
  { std::ofstream f(dir + "cmt_firepar.txt"); std::string s = firepar_block(9, false); s += "1 // extra\n"; f << s; }
  CohortLookup c; c.cmtcode = "CMT09";
  CHECK(ut::exit_status_of([&c, &dir] { c.assignFirePar(dir); }) != 0);   // 19 lines: rejected
  std::remove((dir + "cmt_firepar.txt").c_str()); rmdir(dir.c_str());
}

}
