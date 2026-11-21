// restart-state-carry (Elchin 2026-11-21): state that a chained run keeps in memory across stages
// but a run restarted from a file lost. See docs/WORK_NOTES_restart_state_carry.md.
//   (1) bnf::YearMemory (Soil_Bgc::bnfmem): prev_npp[12], prev_et[12], prev_pft_npp[NUM_PFT], have_prev
//       -> restart vars bnf_prev_npp / bnf_prev_et (Y,X,prevtwelve), bnf_prev_pft_npp (Y,X,pft),
//       bnf_have_prev (Y,X). Old restarts without them keep the in-memory state (legacy).
//   (2) puddle_ice / magic_puddle written for every run (was runon_energy_mode >= 1 only).
// Bug fix: no published source; expected values are exact copies (bit-identical doubles).
#include "doctest.h"
#include "test_helpers.h"
#include "nc_fixtures.h"
#include <cmath>
#include <unistd.h>

namespace {
std::string rc_file(const char* tag) {
  return std::string("/tmp/ut_restart_carry_") + tag + "_" + std::to_string(getpid()) + ".nc";
}
// An "old" restart: Y/X dims and one variable only (no puddle or BNF vars).
std::string rc_old_file(const char* tag, int ny = 2, int nx = 3) {
  std::string fn = rc_file(tag);
  int nc, dy, dx, v;
  ut::ncchk(nc_create(fn.c_str(), NC_CLOBBER | NC_NETCDF4, &nc));
  ut::ncchk(nc_def_dim(nc, "Y", ny, &dy)); ut::ncchk(nc_def_dim(nc, "X", nx, &dx));
  int d2[2] = {dy, dx};
  ut::ncchk(nc_def_var(nc, "dsr", NC_INT, 2, d2, &v));
  ut::ncchk(nc_close(nc));
  return fn;
}
bool has_var(const std::string& fn, const char* name) {
  int nc, v; ut::ncchk(nc_open(fn.c_str(), NC_NOWRITE, &nc));
  bool ok = (nc_inq_varid(nc, name, &v) == NC_NOERR); nc_close(nc); return ok;
}
// One year of monthly NPP/ET (Jan-May ~0, summer peak) and PFT sums, recorded as the model does.
bnf::YearMemory one_year() {
  bnf::YearMemory y;
  const double npp[12] = {0, 0, 0, 0.1, 2.5, 21.7, 38.2, 30.4, 9.9, 0.7, 0, 0};
  for (int m = 0; m < 12; m++) {
    for (int ip = 0; ip < 3; ip++) { y.record_pft(ip, npp[m] * (0.5 + 0.2 * ip)); }
    y.record(m, npp[m], 1.3 * npp[m] + 0.05 * m);
  }
  return y;
}
}

TEST_SUITE("restart_carry") {

TEST_CASE("create_empty_file defines puddle_ice, magic_puddle and the four BNF vars; reinitValue clears them") {
  std::string fn = rc_file("empty");
  RestartData::create_empty_file(fn, 2, 3, 0);
  for (const char* n : {"puddle_ice", "magic_puddle", "bnf_prev_npp", "bnf_prev_et", "bnf_prev_pft_npp", "bnf_have_prev"}) {
    CAPTURE(n); CHECK(has_var(fn, n));
  }
  RestartData r; r.reinitValue();
  CHECK_FALSE(r.has_bnf_state); CHECK(r.bnf_have_prev == MISSING_I);
  CHECK(r.bnf_prev_npp[0] == MISSING_D); CHECK(r.bnf_prev_et[11] == MISSING_D); CHECK(r.bnf_prev_pft_npp[NUM_PFT - 1] == MISSING_D);
  CHECK(r.magic_puddle == MISSING_D); CHECK(r.puddle_ice_mm == MISSING_D);
  std::remove(fn.c_str());
}

TEST_CASE("BNF memory and puddle round trip through a restart file: values per pixel come back exactly") {
  std::string fn = rc_file("rt");
  RestartData::create_empty_file(fn, 2, 3, 0);
  Soil_Bgc a; a.bnfmem = one_year();
  RestartData w; w.reinitValue();
  a.bnf_memory_to_restart(w);
  w.magic_puddle = 3.178487; w.puddle_ice_mm = 0.25;
  w.write_px_bnf_vars(fn, 1, 2); w.write_px_puddle_vars(fn, 1, 2);
  RestartData w2; w2.reinitValue();            // second pixel, no previous year yet
  Soil_Bgc fresh; fresh.bnf_memory_to_restart(w2); w2.magic_puddle = 0.0; w2.puddle_ice_mm = 0.0;
  w2.write_px_bnf_vars(fn, 0, 0); w2.write_px_puddle_vars(fn, 0, 0);

  RestartData r; r.reinitValue(); r.read_px_bnf_vars(fn, 1, 2); r.read_px_puddle_vars(fn, 1, 2);
  REQUIRE(r.has_bnf_state);
  CHECK(r.bnf_have_prev == 1);
  for (int k = 0; k < 12; k++) { CHECK(r.bnf_prev_npp[k] == a.bnfmem.prev_npp[k]); CHECK(r.bnf_prev_et[k] == a.bnfmem.prev_et[k]); }
  for (int ip = 0; ip < NUM_PFT; ip++) { CHECK(r.bnf_prev_pft_npp[ip] == a.bnfmem.prev_pft_npp[ip]); }
  CHECK(r.magic_puddle == 3.178487); CHECK(r.puddle_ice_mm == 0.25);

  RestartData r2; r2.reinitValue(); r2.read_px_bnf_vars(fn, 0, 0);
  REQUIRE(r2.has_bnf_state); CHECK(r2.bnf_have_prev == 0); CHECK(r2.bnf_prev_npp[6] == 0.0);
  RestartData r3; r3.reinitValue(); r3.read_px_bnf_vars(fn, 0, 1); r3.read_px_puddle_vars(fn, 0, 1);  // never written
  CHECK(r3.has_bnf_state); CHECK(r3.bnf_have_prev == MISSING_I); CHECK(r3.bnf_prev_npp[0] == MISSING_D);
  CHECK(r3.magic_puddle == MISSING_D);
  std::remove(fn.c_str());
}

TEST_CASE("restored BNF memory: first-year BNF after a restart equals the chained (in-memory) BNF, month by month") {
  std::string fn = rc_file("bnf");
  RestartData::create_empty_file(fn, 1, 1, 0);
  Soil_Bgc chained; chained.bnfmem = one_year();
  RestartData w; w.reinitValue(); chained.bnf_memory_to_restart(w); w.write_px_bnf_vars(fn, 0, 0);
  RestartData r; r.reinitValue(); r.read_px_bnf_vars(fn, 0, 0);
  Soil_Bgc split;                                   // a fresh process
  CHECK_FALSE(split.bnfmem.have_prev);
  REQUIRE(split.bnf_memory_from_restart(r));
  CHECK(split.bnfmem.have_prev);
  bnf::Params p; p.scheme = 2; p.asym = 0.02;
  double sum = 0.0;
  for (int m = 0; m < 12; m++) {
    CAPTURE(m);
    double bc = bnf::month_bnf(p, chained.bnfmem, m, 1.0), bs = bnf::month_bnf(p, split.bnfmem, m, 1.0);
    CHECK(bs == bc); sum += bs;
    double c1, c2, s1, s2; bnf::month_parts(p, chained.bnfmem, m, 1.0, c1, c2); bnf::month_parts(p, split.bnfmem, m, 1.0, s1, s2);
    CHECK(s1 == c1); CHECK(s2 == c2);
  }
  CHECK(sum > 0.0);                                 // legacy split: 0 for the whole first year
  CHECK(bnf::month_bnf(p, split.bnfmem, 6, 1.0) > 0.0);
  bool act[NUM_PFT] = {true, true, true}; double sc[NUM_PFT], ss[NUM_PFT];
  CHECK(bnf::pft_shares(chained.bnfmem.prev_pft_npp, act, NUM_PFT, sc) == bnf::pft_shares(split.bnfmem.prev_pft_npp, act, NUM_PFT, ss));
  for (int ip = 0; ip < NUM_PFT; ip++) { CHECK(ss[ip] == sc[ip]); }
  // the memory then rolls on exactly as in the chained run
  for (int m = 0; m < 12; m++) { chained.bnfmem.record(m, 1.0 + m, 2.0); split.bnfmem.record(m, 1.0 + m, 2.0); }
  for (int k = 0; k < 12; k++) { CHECK(split.bnfmem.prev_npp[k] == chained.bnfmem.prev_npp[k]); }
  std::remove(fn.c_str());
}

TEST_CASE("old restart without the new vars: nothing read, BNF memory kept (legacy first-year 0 in a fresh process)") {
  std::string fn = rc_old_file("old");
  RestartData r; r.reinitValue(); r.bnf_have_prev = 1; r.bnf_prev_npp[3] = 9.0; r.magic_puddle = 4.0;
  r.read_px_bnf_vars(fn, 1, 2); r.read_px_puddle_vars(fn, 1, 2);
  CHECK_FALSE(r.has_bnf_state); CHECK(r.bnf_have_prev == MISSING_I); CHECK(r.bnf_prev_npp[3] == MISSING_D);
  CHECK(r.magic_puddle == MISSING_D); CHECK(r.puddle_ice_mm == MISSING_D);
  CHECK_FALSE(has_var(fn, "bnf_prev_npp")); CHECK_FALSE(has_var(fn, "magic_puddle"));   // reading adds nothing
  Soil_Bgc fresh;                                    // fresh process: stays without memory, as before the fix
  CHECK_FALSE(fresh.bnf_memory_from_restart(r)); CHECK_FALSE(fresh.bnfmem.have_prev);
  bnf::Params p; p.scheme = 2;
  CHECK(bnf::month_bnf(p, fresh.bnfmem, 6, 1.0) == 0.0);
  Soil_Bgc chained; chained.bnfmem = one_year();     // chained process: in-memory state untouched
  bnf::YearMemory before = chained.bnfmem;
  CHECK_FALSE(chained.bnf_memory_from_restart(r));
  CHECK(chained.bnfmem.have_prev);
  for (int k = 0; k < 12; k++) { CHECK(chained.bnfmem.prev_npp[k] == before.prev_npp[k]); }
  // a pixel never written (fill values) or a bad flag is treated like an old restart
  RestartData f; f.reinitValue(); f.has_bnf_state = true;
  CHECK_FALSE(fresh.bnf_memory_from_restart(f));
  chained.bnf_memory_to_restart(f); f.has_bnf_state = true; f.bnf_have_prev = 7;
  CHECK_FALSE(fresh.bnf_memory_from_restart(f));
  f.bnf_have_prev = 1; f.bnf_prev_et[5] = MISSING_D;
  CHECK_FALSE(fresh.bnf_memory_from_restart(f)); CHECK_FALSE(fresh.bnfmem.have_prev);
  std::remove(fn.c_str());
}

TEST_CASE("serial writers add the puddle and BNF vars to a restart made by an older binary") {
  std::string fn = rc_old_file("upgrade");
  Soil_Bgc a; a.bnfmem = one_year();
  RestartData w; w.reinitValue(); a.bnf_memory_to_restart(w); w.magic_puddle = 1.5; w.puddle_ice_mm = 0.0;
  w.write_px_bnf_vars(fn, 1, 1); w.write_px_puddle_vars(fn, 1, 1);
  RestartData r; r.reinitValue(); r.read_px_bnf_vars(fn, 1, 1); r.read_px_puddle_vars(fn, 1, 1);
  REQUIRE(r.has_bnf_state); CHECK(r.bnf_have_prev == 1);
  CHECK(r.bnf_prev_npp[6] == a.bnfmem.prev_npp[6]); CHECK(r.bnf_prev_pft_npp[2] == a.bnfmem.prev_pft_npp[2]);
  CHECK(r.magic_puddle == 1.5); CHECK(r.puddle_ice_mm == 0.0);
  std::remove(fn.c_str());
}

}
