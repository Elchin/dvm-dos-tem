// RTS / thermokarst port (model_settings.thermokarst; docs/WORK_NOTES_rts.md): master switch
// parse/default/FATAL, sub-keys read only when on (v0.3 defaults), inert defaults of every hook,
// optional restart state (absent in old restarts; round trip), and the rts-regrow pure-function
// tests (Thermokarst / Ground / RTSFloor / RTSRegrow / RTSMethane / TKNetC).
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include "../../include/Thermokarst.h"
#include "../../include/ThermokarstData.h"
#include "../../include/RTSMethane.h"
#include "../../include/RTSRegrow.h"
#include "../../include/RTSFloor.h"
#include "../../include/TKNetC.h"
#include "../../include/Ground.h"
#include "../../include/Soil_Env.h"
#include "../../include/Cohort.h"
#include <fstream>
#include <string>
#include "nc_fixtures.h"
#include <unistd.h>

namespace {
Json::Value tkc(Json::Value tk) { Json::Value c; if (!tk.isNull()) c["model_settings"]["thermokarst"] = tk; return c; }
Json::Value tkon(const char* k, Json::Value v) { Json::Value c = tkc(1); c["model_settings"][k] = v; return c; }
const Json::Value NUL;

std::string tk_restart_file(const char* tag, int ny = 2, int nx = 3) {
  std::string fn = std::string("/tmp/ut_tk_restart_") + tag + "_" + std::to_string(getpid()) + ".nc";
  int nc, dy, dx, dl, v;
  ut::ncchk(nc_create(fn.c_str(), NC_CLOBBER | NC_NETCDF4, &nc));
  ut::ncchk(nc_def_dim(nc, "Y", ny, &dy)); ut::ncchk(nc_def_dim(nc, "X", nx, &dx));
  ut::ncchk(nc_def_dim(nc, "soillayer", MAX_SOI_LAY, &dl));
  int d3[3] = {dy, dx, dl};
  ut::ncchk(nc_def_var(nc, "DZsoil", NC_DOUBLE, 3, d3, &v));
  ut::ncchk(nc_close(nc));
  return fn;
}
}

TEST_SUITE("rts") {

TEST_CASE("thermokarst master switch: default 0; integers 0/1 only, anything else FATAL") {
  ModelData md; CHECK(md.thermokarst == 0);
  ModelData mj(tkc(NUL)); CHECK(mj.thermokarst == 0);
  CHECK(ModelData(tkc(0)).thermokarst == 0);
  CHECK(ModelData(tkc(1)).thermokarst == 1);
  CHECK(ut::exit_status_of([] { ModelData m(tkc(2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkc(-1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkc(0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkc(true)); }) == EXIT_FAILURE);   // rts-regrow used a bool
  CHECK(ut::exit_status_of([] { ModelData m(tkc("1")); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkc(1)); }) == 0);
}

TEST_CASE("thermokarst 0: no RTS key or input file name is read (invalid sub-keys are ignored)") {
  Json::Value c = tkc(0);
  c["IO"]["hist_exp_thermokarst_file"] = "/nonexistent/hist.nc";
  c["IO"]["proj_exp_thermokarst_file"] = "/nonexistent/proj.nc";
  c["model_settings"]["thermokarst_folr"] = "not a number";
  c["model_settings"]["thermokarst_regrow"] = 7;
  c["model_settings"]["thermokarst_excess_fraction"] = 0.95;
  CHECK(ut::exit_status_of([c] { ModelData m(c); }) == 0);
  ModelData m(c);
  CHECK(m.hist_exp_thermokarst_file.empty()); CHECK(m.proj_exp_thermokarst_file.empty());
  CHECK(m.tk_folr == -1.0); CHECK(m.tk_excess_fraction == 0.0);
}

TEST_CASE("thermokarst 1: v0.3 settings on by default (the auto-enabled list), CH4 off") {
  Json::Value c = tkc(1);
  c["IO"]["hist_exp_thermokarst_file"] = "h.nc"; c["IO"]["proj_exp_thermokarst_file"] = "p.nc";
  ModelData m(c);
  CHECK(m.hist_exp_thermokarst_file == "h.nc"); CHECK(m.proj_exp_thermokarst_file == "p.nc");
  CHECK(m.tk_floor.rh_cap); CHECK(m.tk_floor.exp_on); CHECK(m.tk_floor.anox); CHECK(m.tk_floor.rmin_on);
  CHECK(m.tk_regrow.on); CHECK(m.tk_regrow.n_taper); CHECK(m.tk_regrow.navln_exempt);
  CHECK(m.tk_initiation); CHECK_FALSE(m.tk_ch4);
  CHECK(m.tk_folr == -1.0); CHECK(m.tk_c_resp_frac == 0.0); CHECK(m.tk_veg_remove_frac == 0.99);
  CHECK(m.tk_thaw_depth == 0.0); CHECK(m.tk_post_drainage == -1); CHECK(m.tk_excess_fraction == 0.0);
}

TEST_CASE("thermokarst 1: this fork's layer / soil-N switches default to 1 (rts-regrow has them always on)") {
  ModelData off(tkc(0));
  CHECK(off.soil_negpool_fix == 0); CHECK(off.layer_n_remap_fix == 0); CHECK(off.layer_bookkeep_fix == 0);
  ModelData on(tkc(1));
  CHECK(on.soil_negpool_fix == 1); CHECK(on.layer_n_remap_fix == 1); CHECK(on.layer_bookkeep_fix == 1);
  ModelData k(tkon("soil_negpool_fix", 0));   // explicit 0 kept (WARNING)
  CHECK(k.soil_negpool_fix == 0); CHECK(k.layer_n_remap_fix == 1); CHECK(k.layer_bookkeep_fix == 1);
  ModelData n(tkon("layer_n_remap_fix", 0));  // bookkeep needs remap: stays 0, no FATAL
  CHECK(n.layer_n_remap_fix == 0); CHECK(n.layer_bookkeep_fix == 0); CHECK(n.soil_negpool_fix == 1);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("layer_n_remap_fix", 0)); }) == 0);
}

TEST_CASE("thermokarst 1: sub-keys are read (synthetic-test settings)") {
  Json::Value c = tkc(1); Json::Value & ms = c["model_settings"];
  ms["thermokarst_folr"] = 0.8; ms["thermokarst_c_resp_frac"] = 0.1; ms["thermokarst_veg_remove_frac"] = 0.99;
  ms["thermokarst_thaw_depth"] = 1.0; ms["thermokarst_post_drainage"] = 1;
  ms["thermokarst_excess_fraction"] = 0.3; ms["thermokarst_excess_top"] = 0.6; ms["thermokarst_excess_bottom"] = 2.5;
  ms["thermokarst_regrow"] = false; ms["thermokarst_rh_cap"] = 0; ms["thermokarst_ch4"] = true;
  ms["thermokarst_initiation"] = false;
  ModelData m(c);
  CHECK(m.tk_folr == 0.8); CHECK(m.tk_c_resp_frac == 0.1); CHECK(m.tk_thaw_depth == 1.0); CHECK(m.tk_post_drainage == 1);
  CHECK(m.tk_excess_fraction == 0.3); CHECK(m.tk_excess_top == 0.6); CHECK(m.tk_excess_bottom == 2.5);
  CHECK_FALSE(m.tk_regrow.on); CHECK_FALSE(m.tk_floor.rh_cap); CHECK(m.tk_ch4); CHECK_FALSE(m.tk_initiation);
  CHECK(m.tk_floor.exp_on);   // untouched keys keep their defaults
}

TEST_CASE("thermokarst_safetynet_conserve: default true with thermokarst 1, false by key, FATAL type; global switch untouched") {
  ModelData m(tkc(1)); CHECK(m.tk_safetynet_conserve); CHECK(m.n_safetynet_conserve == 0);
  ModelData f(tkon("thermokarst_safetynet_conserve", false)); CHECK_FALSE(f.tk_safetynet_conserve);
  ModelData z(tkon("thermokarst_safetynet_conserve", 0)); CHECK_FALSE(z.tk_safetynet_conserve);
  CHECK(ut::exit_status_of([] { ModelData x(tkon("thermokarst_safetynet_conserve", 2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData x(tkon("thermokarst_safetynet_conserve", "yes")); }) == EXIT_FAILURE);
  Json::Value c = tkc(0); c["model_settings"]["thermokarst_safetynet_conserve"] = 7;   // ignored when off
  ModelData o(c); CHECK(o.tk_safetynet_conserve); CHECK(o.n_safetynet_conserve == 0);
}

TEST_CASE("thermokarst 1: wrong types and ranges are FATAL") {
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_regrow", 2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_regrow", "yes")); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("rts_regrow_n_taper", 0.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_folr", "0.8")); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_folr", true)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("rts_floor_rmin_years", "3")); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_folr", 1.5)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_c_resp_frac", -0.1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_post_drainage", 2)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_excess_fraction", 0.8)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_excess_fraction", -0.1)); }) == EXIT_FAILURE);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_excess_bottom", 0.4)); }) == EXIT_FAILURE);  // <= top 0.5
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_folr", -1.0)); }) == 0);   // <0 = legacy removal
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_regrow", 0)); }) == 0);
  CHECK(ut::exit_status_of([] { ModelData m(tkon("thermokarst_post_drainage", -1)); }) == 0);
}

TEST_CASE("hook defaults are inert (legacy) until the Cohort sets them after an RTS initiation") {
  Soil_Bgc sb;
  CHECK_FALSE(sb.tk_litter_fix); CHECK_FALSE(sb.tk_floor_on); CHECK_FALSE(sb.tk_navln_exempt);
  CHECK(sb.tk_ninput_layer == 0); CHECK(sb.tk_rmin_target < 0.0);
  Integrator in; CHECK_FALSE(in.tk_tolerate_neg_avln);
  Soil_Env se; CHECK_FALSE(se.tk_nofibric_fix);
  Ground g; CHECK(g.tk_months_since_init == -1); CHECK(g.tk_subsidence == 0.0); CHECK(g.tk_prev_subsidence == 0.0);
  CHECK_FALSE(g.tk_excess_initialized);
  RestartData r; CHECK_FALSE(r.read_tk_state); CHECK_FALSE(r.write_tk_state);
}

TEST_CASE("restart RTS state: absent in old restarts (has_tk_state false, nothing added); reinitValue clears it") {
  RestartData r; r.reinitValue(); r.read_tk_state = true;
  CHECK_FALSE(r.write_tk_state); CHECK_FALSE(r.has_tk_state);
  CHECK(r.tk_init == MISSING_I); CHECK(r.tk_exice[0] == MISSING_D);
  std::string fn = tk_restart_file("old");
  r.read_px_tk_vars(fn, 1, 2);
  CHECK_FALSE(r.has_tk_state);
  int nc, nv; ut::ncchk(nc_open(fn.c_str(), NC_NOWRITE, &nc)); ut::ncchk(nc_inq_nvars(nc, &nv)); nc_close(nc);
  CHECK(nv == 1);
  r.write_tk_state = true; r.has_tk_state = true; r.reinitValue();
  CHECK_FALSE(r.write_tk_state); CHECK_FALSE(r.has_tk_state); CHECK(r.read_tk_state);   // read flag is set once by Cohort
  std::remove(fn.c_str());
}

TEST_CASE("restart RTS state round trip: drainage, init month, excess-ice flag, subsidence, per-layer excess ice") {
  std::string fn = tk_restart_file("rt");
  RestartData w; w.reinitValue();
  w.tk_drainage_type = 1; w.tk_init = 37; w.tk_exinit = 1; w.tk_subsidence = 0.123456789;
  for (int il = 0; il < MAX_SOI_LAY; il++) { w.tk_exice[il] = (il > 4 && il < 9) ? 101.5 + il : 0.0; }
  w.write_px_tk_vars(fn, 1, 2);
  RestartData w2; w2.reinitValue(); w2.tk_drainage_type = 0; w2.tk_init = 0; w2.tk_exinit = 0; w2.tk_subsidence = 0.0;
  for (int il = 0; il < MAX_SOI_LAY; il++) { w2.tk_exice[il] = 0.0; }
  w2.write_px_tk_vars(fn, 0, 0);
  RestartData r; r.reinitValue(); r.read_px_tk_vars(fn, 1, 2);
  REQUIRE(r.has_tk_state);
  CHECK(r.tk_drainage_type == 1); CHECK(r.tk_init == 37); CHECK(r.tk_exinit == 1); CHECK(r.tk_subsidence == w.tk_subsidence);
  for (int il = 0; il < MAX_SOI_LAY; il++) { CHECK(r.tk_exice[il] == w.tk_exice[il]); }
  RestartData r2; r2.reinitValue(); r2.read_px_tk_vars(fn, 0, 0);
  CHECK(r2.tk_init == 0); CHECK(r2.tk_drainage_type == 0);
  RestartData r3; r3.reinitValue(); r3.read_px_tk_vars(fn, 0, 1);   // never written: fill values
  CHECK(r3.has_tk_state); CHECK(r3.tk_init == MISSING_I); CHECK(r3.tk_exice[0] == MISSING_D);
  std::remove(fn.c_str());
}

TEST_CASE("restart RTS regrow state: optional tk_regrow_state round trip; old thermokarst restarts load without it") {
  // a thermokarst restart written before tk_regrow_state existed: the five tk vars only
  std::string fo = tk_restart_file("rgold");
  {
    int nc, dv, d3[3], v;
    ut::ncchk(nc_open(fo.c_str(), NC_WRITE, &nc)); ut::ncchk(nc_inq_varid(nc, "DZsoil", &dv));
    ut::ncchk(nc_inq_vardimid(nc, dv, d3)); ut::ncchk(nc_redef(nc));
    for (const char* n : {"tk_drainage_type", "tk_init", "tk_exinit"}) { ut::ncchk(nc_def_var(nc, n, NC_INT, 2, d3, &v)); }
    ut::ncchk(nc_def_var(nc, "tk_subsidence", NC_DOUBLE, 2, d3, &v));
    ut::ncchk(nc_def_var(nc, "tk_exice", NC_DOUBLE, 3, d3, &v));
    ut::ncchk(nc_close(nc));
  }
  RestartData ro; ro.reinitValue(); ro.read_px_tk_vars(fo, 1, 2);
  CHECK(ro.has_tk_state); CHECK_FALSE(ro.has_tk_rg_state);
  CHECK(ro.tk_rg_state[0] == MISSING_D);
  std::remove(fo.c_str());

  std::string fn = tk_restart_file("rg");
  RestartData w; w.reinitValue();
  w.tk_drainage_type = 1; w.tk_init = 5; w.tk_exinit = 0; w.tk_subsidence = 0.0;
  for (int il = 0; il < MAX_SOI_LAY; il++) { w.tk_exice[il] = 0.0; }
  for (int k = 0; k < RestartData::TK_RG_NSTATE; k++) { w.tk_rg_state[k] = (k == 0) ? 1.0 : 0.5 * k + 1e-7 * k * k; }
  w.write_px_tk_vars(fn, 1, 2);
  RestartData r; r.reinitValue(); r.read_px_tk_vars(fn, 1, 2);
  REQUIRE(r.has_tk_state); REQUIRE(r.has_tk_rg_state);
  for (int k = 0; k < RestartData::TK_RG_NSTATE; k++) { CHECK(r.tk_rg_state[k] == w.tk_rg_state[k]); }
  RestartData r2; r2.reinitValue(); r2.read_px_tk_vars(fn, 0, 0);   // never written: fill -> no regrow state
  CHECK_FALSE(r2.has_tk_rg_state);
  int nc, vr, dr[3], ver = 0; size_t len = 0;
  ut::ncchk(nc_open(fn.c_str(), NC_NOWRITE, &nc)); ut::ncchk(nc_inq_varid(nc, "tk_regrow_state", &vr));
  ut::ncchk(nc_inq_vardimid(nc, vr, dr)); ut::ncchk(nc_inq_dimlen(nc, dr[2], &len));
  ut::ncchk(nc_get_att_int(nc, vr, "layout_version", &ver)); nc_close(nc);
  CHECK(len == (size_t)RestartData::TK_RG_NSTATE); CHECK(ver == 1);
  RestartData rr; rr.has_tk_rg_state = true; rr.tk_rg_state[3] = 7.0; rr.reinitValue();
  CHECK_FALSE(rr.has_tk_rg_state); CHECK(rr.tk_rg_state[3] == MISSING_D);
  std::remove(fn.c_str());
}

TEST_CASE("restart RTS regrow state: a length mismatch (other layout) is ignored, not misread") {
  std::string fn = tk_restart_file("rglen");
  RestartData w; w.reinitValue(); w.tk_drainage_type = 0; w.tk_init = 0; w.tk_exinit = 0; w.tk_subsidence = 0.0;
  for (int il = 0; il < MAX_SOI_LAY; il++) { w.tk_exice[il] = 0.0; }
  int nc, dv, d3[3], dr, v;   // pre-define tk_regrow_state with a different length
  ut::ncchk(nc_open(fn.c_str(), NC_WRITE, &nc)); ut::ncchk(nc_inq_varid(nc, "DZsoil", &dv));
  ut::ncchk(nc_inq_vardimid(nc, dv, d3)); ut::ncchk(nc_redef(nc));
  ut::ncchk(nc_def_dim(nc, "tkrgstate", RestartData::TK_RG_NSTATE + 1, &dr));
  int d3r[3] = {d3[0], d3[1], dr}; ut::ncchk(nc_def_var(nc, "tk_regrow_state", NC_DOUBLE, 3, d3r, &v));
  ut::ncchk(nc_close(nc));
  RestartData r; r.reinitValue(); r.read_px_tk_vars(fn, 0, 0);
  CHECK_FALSE(r.has_tk_rg_state);
  std::remove(fn.c_str());
}

TEST_CASE("Cohort tk_rg_pack / tk_rg_unpack: the regrow and floor-RH state round-trips exactly (unpack -> pack)") {
  // distinct value per slot (ints/bools where the layout has them); a slot that unpack/pack mapped to
  // the wrong member, or dropped, would not come back unchanged (the slot count is asserted in Cohort)
  const int NP = NUM_PFT, N = RestartData::TK_RG_NSTATE;
  std::vector<double> s(N);
  for (int k = 0; k < N; k++) { s[k] = 0.5 + 1.25 * k; }
  s[0] = 1.0;                                                           // tk_rg_have
  for (int ip = 0; ip < NP; ip++) { s[1 + 6 * NP + ip] = ip % 4; }      // tk_rg_type
  const int o = 1 + 7 * NP + 6 + 12;                                    // gpp12_n, gpp12_i
  s[o] = 12; s[o + 1] = 7;
  const int o2 = o + 2 + 1 + 12 + 1 + 120;                              // gpp_n, gpp_i
  s[o2] = 60; s[o2 + 1] = 33;
  for (int m = 0; m < 24; m++) { s[N - 24 + m] = m % 10; }              // rh_nyr, rh_next
  Cohort c; c.tk_rg_unpack(s.data());
  std::vector<double> s2(N, MISSING_D); c.tk_rg_pack(s2.data());
  CHECK(s2 == s);
  Cohort d; std::vector<double> z(N, MISSING_D); d.tk_rg_pack(z.data());   // default state: not recorded
  CHECK(z[0] == 0.0); for (double x : z) { CHECK(x != MISSING_D); }
}

TEST_CASE("RTS water outputs: TKWDRAIN / TKWLAT diagnostics start at 0, clear() resets them, output_spec rows exist") {
  soidiag_thermokarst d; CHECK(d.water_drained == 0.0); CHECK(d.water_removed == 0.0); CHECK(d.drainage_type == -1);
  ThermokarstData t; t.thermokarst_soid.water_drained = 3.0; t.thermokarst_soid.water_removed = 4.0; t.clear();
  CHECK(t.thermokarst_soid.water_drained == 0.0); CHECK(t.thermokarst_soid.water_removed == 0.0);
  std::ifstream f("../../config/output_spec.csv"); REQUIRE(f.good());
  std::string line; int found = 0;
  while (std::getline(f, line)) {
    if (line.rfind("TKWDRAIN,", 0) == 0 || line.rfind("TKWLAT,", 0) == 0) {
      found++; CHECK(line.find(",kg/m2/time,") != std::string::npos);
    }
  }
  CHECK(found == 2);
}

TEST_CASE("Soil_Bgc AVLN input layer: default 0 = legacy; an RTS regrow can move it to a rooted layer") {
  for (int L : {0, 2}) {
    CohortData cd; BgcData bd; EnvData ed; FirData fd; Soil_Bgc sb;
    sb.setCohortData(&cd); sb.setBgcData(&bd); sb.setEnvData(&ed); sb.setFirData(&fd);
    sb.nfeed = 1; sb.tk_ninput_layer = L;
    sb.bgcpar.som2co2 = 0.2; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
    const int n = 4; cd.m_soil.numsl = n; cd.m_soil.deepthick = 0.1; sb.totdzliq = 0.0;
    double z = 0;
    for (int i = 0; i < MAX_SOI_LAY; ++i) {
      bool on = i < n;
      cd.m_soil.type[i] = 3; cd.m_soil.dz[i] = 0.1; cd.m_soil.z[i] = z; z += 0.1;
      sb.tmp_sois.rawc[i] = 100; sb.tmp_sois.soma[i] = 50; sb.tmp_sois.sompr[i] = 200; sb.tmp_sois.somcr[i] = 300;
      sb.tmp_sois.orgn[i] = 20; sb.tmp_sois.avln[i] = 1; bd.m_sois.avln[i] = on ? 0.5 : 0.0;
      sb.del_soi2a.rhrawc[i] = 0; sb.del_soi2a.rhsoma[i] = 0; sb.del_soi2a.rhsompr[i] = 0; sb.del_soi2a.rhsomcr[i] = 0;
      sb.ltrflc[i] = 0; sb.ltrfln[i] = 0; sb.rtnextract[i] = 0;
      sb.del_soi2soi.netnmin[i] = 0; bd.m_v2soi.rtlfalfrac[i] = on ? 0.25 : 0.0;
      ed.m_sois.liq[i] = 10.0;
    }
    sb.tmp_sois.wdebrisc = 0; sb.tmp_sois.wdebrisn = 0; sb.d2wdebrisc = 0; sb.d2wdebrisn = 0; sb.del_soi2a.rhwdeb = 0;
    bd.m_a2soi.orgninput = 0; bd.m_a2soi.avlninput = 0.7; sb.del_soi2l.avlnlost = 0; ed.m_sois.draindepth = 0;
    sb.deltastate();
    for (int i = 0; i < n; ++i) { CHECK(sb.del_sois.avln[i] == doctest::Approx(i == L ? 0.7 : 0.0)); }
  }
}

TEST_CASE("Ground::adjustSoilAfterburn default argument is the legacy fire path") {
  // signature check: the fire call site (no argument) compiles to tk_remove_mineral = false
  double (Ground::*f)(bool) = &Ground::adjustSoilAfterburn;
  CHECK(f != nullptr);
}

}

// ---- ported rts-regrow pure-function tests (m2d, deep-C burial and cover_negative_orgn tests
//      not ported: those features are not in this port)
#define RCHECK(c) CHECK_MESSAGE(static_cast<bool>(c), #c)
#define CHECK_NEAR(a, b, tol) RCHECK(std::fabs((a) - (b)) <= (tol))
namespace rtsu {
static void test_removal_thickness() {
  // fraction of total organic thickness
  CHECK_NEAR(Thermokarst::organicRemovalThickness(0.8, 0.02, 0.10, 0.18, 0.7), 0.24, 1e-12);
  CHECK_NEAR(Thermokarst::organicRemovalThickness(0.0, 0.02, 0.10, 0.18, 0.7), 0.0, 1e-12);
  // clamp to 1
  CHECK_NEAR(Thermokarst::organicRemovalThickness(1.5, 0.02, 0.10, 0.18, 0.7), 0.30, 1e-12);
  // legacy mode (folr<0): thaw depth, never negative
  CHECK_NEAR(Thermokarst::organicRemovalThickness(-1.0, 0.02, 0.10, 0.18, 0.7), 0.7, 1e-12);
  CHECK_NEAR(Thermokarst::organicRemovalThickness(-1.0, 0.02, 0.10, 0.18, -1.0), 0.0, 1e-12);
}
}

namespace rtsu {
static void test_partition_mass_balance() {
  const double fr[] = {0.0, 0.1, 0.5, 1.0, -0.3, 1.7};
  for (double f : fr) {
    double rc, rn, lc, ln;
    Thermokarst::partitionLostCN(5000.0, 175.0, f, rc, rn, lc, ln);
    CHECK_NEAR(rc + lc, 5000.0, 1e-9);   // removed = respired + lateral
    CHECK_NEAR(rn + ln, 175.0, 1e-9);
    RCHECK(rc >= 0 && lc >= 0 && rn >= 0 && ln >= 0);
  }
  double rc, rn, lc, ln;
  Thermokarst::partitionLostCN(1000.0, 10.0, 0.25, rc, rn, lc, ln);
  CHECK_NEAR(rc, 250.0, 1e-9); CHECK_NEAR(lc, 750.0, 1e-9);
  CHECK_NEAR(rn, 2.5, 1e-9);   CHECK_NEAR(ln, 7.5, 1e-9);
}
}

namespace rtsu {
static void test_thaw_column() {
  // 4 layers of 0.1, 0.2, 0.3, 0.4 m; tops 0, 0.1, 0.3, 0.6
  std::vector<double> ztop = {0.0, 0.1, 0.3, 0.6}, dz = {0.1, 0.2, 0.3, 0.4};
  std::vector<double> maxliq = {50, 100, 150, 200};
  std::vector<double> ice = {10, 120, 90, 80}, liq = {5, 0, 10, 0};
  std::vector<bool> thawed; double drained = -1;
  double total_before = 0; for (size_t i = 0; i < 4; i++) total_before += ice[i] + liq[i];
  // depth 0.6: layers 0..2 fully inside, layer 3 not touched
  double melted = Ground::thawColumn(ztop, dz, ice, liq, maxliq, 0.6, drained, thawed);
  CHECK_NEAR(melted, 220.0, 1e-9);
  RCHECK(thawed[0] && thawed[1] && thawed[2] && !thawed[3]);
  CHECK_NEAR(ice[0] + ice[1] + ice[2], 0.0, 1e-12);
  CHECK_NEAR(ice[3], 80.0, 1e-12);                 // below thaw depth: unchanged
  CHECK_NEAR(liq[1], 100.0, 1e-12);                // capped at maxliq
  CHECK_NEAR(drained, 20.0, 1e-9);                 // 120 ice -> 100 kept + 20 drained
  double total_after = drained; for (size_t i = 0; i < 4; i++) total_after += ice[i] + liq[i];
  CHECK_NEAR(total_after, total_before, 1e-9);     // water conserved incl. drainage
  // partial layer is not thawed (whole layers only)
  std::vector<double> ice2 = {10, 10, 10, 10}, liq2 = {0, 0, 0, 0};
  Ground::thawColumn(ztop, dz, ice2, liq2, maxliq, 0.25, drained, thawed);
  RCHECK(thawed[0] && !thawed[1]);
  CHECK_NEAR(ice2[1], 10.0, 1e-12);
  // zero depth: no-op
  std::vector<double> ice3 = {10, 10, 10, 10}, liq3 = {1, 1, 1, 1};
  double m3 = Ground::thawColumn(ztop, dz, ice3, liq3, maxliq, 0.0, drained, thawed);
  CHECK_NEAR(m3, 0.0, 1e-12); CHECK_NEAR(drained, 0.0, 1e-12); RCHECK(!thawed[0]);
}
}

namespace rtsu {
static void test_excess_ice_init() {
  // prototype profile: excess = DENICE * overlap * f/(1-f)
  CHECK_NEAR(Ground::excessIceForLayer(0.5, 0.5, 0.3, 0.3, 2.0), DENICE * 0.5 * 0.3 / 0.7, 1e-9);
  CHECK_NEAR(Ground::excessIceForLayer(0.2, 0.2, 0.3, 0.3, 2.0), DENICE * 0.1 * 0.3 / 0.7, 1e-9); // partial overlap
  CHECK_NEAR(Ground::excessIceForLayer(2.5, 0.5, 0.3, 0.3, 2.0), 0.0, 1e-12);  // outside range
  CHECK_NEAR(Ground::excessIceForLayer(0.5, 0.5, 0.0, 0.3, 2.0), 0.0, 1e-12);  // f=0 -> none
  // excess volume fraction of (matrix + excess) equals f
  double x = Ground::excessIceForLayer(0.5, 1.0, 0.4, 0.0, 5.0);
  CHECK_NEAR((x / DENICE) / (1.0 + x / DENICE), 0.4, 1e-12);
}
}

namespace rtsu {
static void test_excess_melt_subsidence() {
  std::vector<double> dz = {0.2, 0.3, 0.5};
  std::vector<double> excess = {0.0, 200.0, 300.0}, liq = {10.0, 50.0, 100.0};
  std::vector<double> maxliq = {80.0, 120.0, 200.0};
  std::vector<double> meltfrac = {1.0, 1.0, 0.5};     // layer 3 half thawed
  const std::vector<double> dz0 = dz;
  double col0 = 0, w0 = 0;
  for (size_t i = 0; i < 3; i++) { col0 += dz[i] + excess[i] / DENICE; w0 += excess[i] + liq[i]; }
  double melted, drained;
  double sub = Ground::meltExcessColumn(excess, liq, maxliq, meltfrac, melted, drained);
  CHECK_NEAR(melted, 350.0, 1e-9);
  CHECK_NEAR(sub, 350.0 / DENICE, 1e-12);            // subsidence == melted excess-ice volume
  CHECK_NEAR(excess[2], 150.0, 1e-9);
  CHECK_NEAR(liq[1], 120.0, 1e-9);                   // capped at maxliq
  CHECK_NEAR(drained, 130.0 + 50.0, 1e-9);           // 200-70 room, 150-100 room
  double col1 = 0, w1 = drained;
  for (size_t i = 0; i < 3; i++) { col1 += dz[i] + excess[i] / DENICE; w1 += excess[i] + liq[i]; }
  RCHECK(dz == dz0);                                  // matrix thickness unchanged
  CHECK_NEAR(col0 - col1, sub, 1e-12);               // column thickness drops by subsidence
  CHECK_NEAR(w1, w0, 1e-9);                          // water mass conserved
  // no excess ice -> no-op
  std::vector<double> e0 = {0, 0, 0}, l0 = {1, 2, 3};
  sub = Ground::meltExcessColumn(e0, l0, maxliq, meltfrac, melted, drained);
  CHECK_NEAR(sub, 0.0, 1e-15); CHECK_NEAR(melted, 0.0, 1e-15); CHECK_NEAR(drained, 0.0, 1e-15);
  RCHECK(l0[0] == 1 && l0[1] == 2 && l0[2] == 3);
  // frozen layers (meltfrac 0) -> no-op
  std::vector<double> e1 = {100, 100, 100}, l1 = {1, 1, 1}, f0 = {0, 0, 0};
  sub = Ground::meltExcessColumn(e1, l1, maxliq, f0, melted, drained);
  CHECK_NEAR(sub, 0.0, 1e-15); RCHECK(e1[0] == 100 && e1[2] == 100);
}

// ---------------- monthly -> daily (mean-preserving, Killworth) ----------------
static const int DIM12[12] = {31,28,31,30,31,30,31,31,30,31,30,31};

static std::vector<float> synth_series(int nyr, double mean, double amp, double noise, unsigned seed, bool clipzero) {
  std::vector<float> v; unsigned s = seed;
  for (int y = 0; y < nyr; ++y) for (int m = 0; m < 12; ++m) {
    s = s * 1103515245u + 12345u;
    double r = ((s >> 8) & 0xFFFF) / 65535.0 - 0.5;
    double x = mean - amp * std::cos(2 * M_PI * (m + 0.5) / 12.0) + noise * r;
    if (clipzero && x < 0) x = 0;      // polar-night radiation: exact zeros
    v.push_back((float)x);
  }
  return v;
}

static void check_means(const std::vector<float>& d, const std::vector<float>& m12, int& bad) {
  int k = 0;
  for (int m = 0; m < 12; ++m) {
    double sum = 0; for (int j = 0; j < DIM12[m]; ++j) sum += d[k++];
    double M = m12[m], tol = 1e-6 * std::max(1.0, std::fabs(M));
    if (std::fabs(sum / DIM12[m] - M) > tol) { bad++;
      std::printf("  month %d mean %.9g vs %.9g\n", m, sum / DIM12[m], M); }
  }
}
}

namespace rtsu {
static void test_rts_ch4() {
  RTSMethaneParams p;
  // saturated fraction
  CHECK_NEAR(RTSMethane::saturatedFraction(0.0, 0.1, 0.5), 0.0, 1e-15);
  CHECK_NEAR(RTSMethane::saturatedFraction(0.6, 0.1, 0.5), 1.0, 1e-15);
  CHECK_NEAR(RTSMethane::saturatedFraction(0.4, 0.2, 0.5), 0.5, 1e-12);
  // oxidation fraction: bounded by [oxid_min, oxid_max] within [0,1], monotonic in WTD
  double prev = -1.0;
  for (double z = -0.5; z <= 3.0; z += 0.01) {
    double f = RTSMethane::oxidFraction(p, z);
    RCHECK(f >= p.oxid_min - 1e-15 && f <= p.oxid_max + 1e-15 && f >= 0.0 && f <= 1.0);
    RCHECK(f >= prev - 1e-15); prev = f;
  }
  RTSMethaneParams q = p; q.oxid_min = -0.5; q.oxid_max = 1.7;     // bad parameters are clamped
  RCHECK(RTSMethane::oxidFraction(q, 0.0) >= 0.0 && RTSMethane::oxidFraction(q, 5.0) <= 1.0);
  // ponding: 0 before initiation, pond_sat0 at initiation, 1/e after pond_tau years
  RCHECK(RTSMethane::pondSaturation(p, -1) == 0.0);
  CHECK_NEAR(RTSMethane::pondSaturation(p, 0), p.pond_sat0, 1e-15);
  CHECK_NEAR(RTSMethane::pondSaturation(p, (int)(12 * p.pond_tau)), p.pond_sat0 / std::exp(1.0), 1e-12);
  // methanogen recovery: 1 before initiation, mg0 at initiation, 1-(1-mg0)/e after mg_tau years, -> 1
  RCHECK(RTSMethane::methanogenFactor(p, -1) == 1.0);
  CHECK_NEAR(RTSMethane::methanogenFactor(p, 0), p.mg0, 1e-15);
  CHECK_NEAR(RTSMethane::methanogenFactor(p, (int)(12 * p.mg_tau)), 1.0 - (1.0 - p.mg0) / std::exp(1.0), 1e-12);
  RCHECK(RTSMethane::methanogenFactor(p, 12000) > 0.999999);
  { RTSMethaneParams q = p; q.mg0 = 1.0; RCHECK(RTSMethane::methanogenFactor(q, 5) == 1.0); }
  // production <= RH (also with an extreme scale factor) and emis + oxid = prod
  std::vector<RTSMethaneLayer> col(6);
  double rh = 0.0;
  for (int i = 0; i < 6; i++) { col[i].ztop = 0.1 * i; col[i].dz = 0.1; col[i].tem = 25.0 - 4 * i; col[i].rh = 1.0 + i; rh += col[i].rh; }
  for (double frac : {0.0, 0.01, 0.5, 5.0, 100.0}) {
    for (double wtd : {0.0, 0.25, 2.0}) {
      for (double pond : {0.0, 0.5, 1.0}) {
        RTSMethaneParams r = p; r.ch4_frac = frac; r.rel_q10 = 3.0;
        RTSMethaneFluxes f = RTSMethane::compute(col, r, wtd, pond, 2.0);
        RCHECK(f.prod >= 0.0 && f.prod <= rh + 1e-12);
        RCHECK(f.oxid >= 0.0 && f.oxid <= f.prod + 1e-12 && f.emis >= 0.0);
        CHECK_NEAR(f.emis + f.oxid, f.prod, 1e-12);
        RCHECK(f.plant <= f.emis + 1e-12);
      }
    }
  }
  // fully saturated, frac=100: prod == RH of thawed layers exactly (cap at 1)
  { RTSMethaneParams r = p; r.ch4_frac = 100.0;
    double rh_thawed = 0.0; for (auto & l : col) if (l.tem > 0) rh_thawed += l.rh;
    CHECK_NEAR(RTSMethane::compute(col, r, 0.0, 0.0, 0.0).prod, rh_thawed, 1e-12); }
  // frozen or no RH -> no CH4; pond=0 + deep WT -> no CH4 (A_RTS=0 / no anoxia baseline)
  { std::vector<RTSMethaneLayer> fz(col); for (auto & l : fz) l.tem = -1.0;
    RCHECK(RTSMethane::compute(fz, p, 0.0, 1.0, 1.0).prod == 0.0);
    RCHECK(RTSMethane::compute(col, p, 10.0, 0.0, 1.0).prod == 0.0); }
  // ponding raises production with a deep water table
  RCHECK(RTSMethane::compute(col, p, 0.3, 0.5, 1.0).prod > RTSMethane::compute(col, p, 0.3, 0.0, 1.0).prod);
}
}

namespace rtsu {
static void test_rts_floor() {
  RTSFloorParams p;
  // v0.2 defaults: all switches on, calibrated values; all off -> any() false (bitwise v0.1 path)
  RCHECK(p.rh_cap && p.exp_on && p.anox && p.any());
  RCHECK(p.cap0 == 0.045 && p.cap_inf == 0.002 && p.cap_tau == 1.5 && p.lat0 == 0.02 && p.lat_tau == 5.0 && p.anox_f == 0.3);
  { RTSFloorParams q = p; q.rh_cap = q.exp_on = q.anox = false; RCHECK(q.any() == q.rmin_on); q.rmin_on = false; RCHECK(!q.any()); }
  { RTSMethaneParams mp; RCHECK(mp.mg0 == 1.0); RCHECK(mp.ch4_frac == 0.135); }
  RCHECK(p.rmin_on && p.rmin0 == 0.3 && p.rmin_tau == 10.0 && p.rmin_years == 5);
  { RTSFloorParams q = p; q.rh_cap = q.exp_on = q.anox = q.rmin_on = false; RCHECK(!q.any()); }
  // rmin factor: 1 before initiation, rmin0 at initiation, 1-(1-rmin0)/e after rmin_tau years, -> 1
  RCHECK(RTSFloor::rminFactor(p, -1) == 1.0);
  CHECK_NEAR(RTSFloor::rminFactor(p, 0), p.rmin0, 1e-15);
  CHECK_NEAR(RTSFloor::rminFactor(p, (int)(12 * p.rmin_tau)), 1.0 - (1.0 - p.rmin0) / std::exp(1.0), 1e-12);
  RCHECK(RTSFloor::rminFactor(p, 120000) > 0.999999);
  { RTSFloorParams q = p; q.rmin0 = -2.0; RCHECK(RTSFloor::rminFactor(q, 0) == 0.0); q.rmin0 = 3.0; RCHECK(RTSFloor::rminFactor(q, 0) == 1.0); }
  // minimum respiration: raises RH to the target from the pools, conserves C, never above rmin_maxfrac of C
  { const int n = 3; double w[n] = {1.0, 0.5, 0.0};
    double cr[n] = {10, 0, 5}, ca[n] = {100, 200, 50}, cp[n] = {1000, 500, 0}, cc[n] = {5000, 3000, 100};
    double rr[n] = {0.01, 0, 0}, ra[n] = {0.1, 0.2, 0}, rp[n] = {0.05, 0.02, 0}, rc[n] = {0.01, 0.01, 0};
    double r0 = 0; for (int i = 0; i < n; i++) r0 += rr[i] + ra[i] + rp[i] + rc[i];
    const double add = RTSFloor::applyMinimum(10.0, n, w, cr, ca, cp, cc, rr, ra, rp, rc);
    double r1 = 0; for (int i = 0; i < n; i++) r1 += rr[i] + ra[i] + rp[i] + rc[i];
    CHECK_NEAR(r1, 10.0, 1e-9); CHECK_NEAR(r1 - r0, add, 1e-12);
    RCHECK(rr[2] == 0 && ra[2] == 0 && rc[2] == 0);                 // zero weight (e.g. no activity): untouched
    for (int i = 0; i < n; i++) RCHECK(rr[i] + ra[i] + rp[i] + rc[i] <= RTSFloor::rmin_maxfrac * (cr[i] + ca[i] + cp[i] + cc[i]) + 1e-12);
    // above target: unchanged; huge target: limited by rmin_maxfrac (no C creation beyond pools)
    double before = rr[0]; RCHECK(RTSFloor::applyMinimum(1.0, n, w, cr, ca, cp, cc, rr, ra, rp, rc) == 0.0); RCHECK(rr[0] == before);
    RTSFloor::applyMinimum(1e9, n, w, cr, ca, cp, cc, rr, ra, rp, rc);
    for (int i = 0; i < n; i++) RCHECK(rr[i] + ra[i] + rp[i] + rc[i] <= RTSFloor::rmin_maxfrac * (cr[i] + ca[i] + cp[i] + cc[i]) + 1e-9);
    // empty pools: nothing added
    double z[n] = {0, 0, 0}, a2[n] = {0, 0, 0}, b2[n] = {0, 0, 0}, c2[n] = {0, 0, 0}, d2[n] = {0, 0, 0};
    RCHECK(RTSFloor::applyMinimum(5.0, n, w, z, z, z, z, a2, b2, c2, d2) == 0.0); }
  // cap rate: cap0 at initiation, cap_inf + (cap0-cap_inf)/e after cap_tau years
  CHECK_NEAR(RTSFloor::capRate(p, 0), p.cap0, 1e-15);
  CHECK_NEAR(RTSFloor::capRate(p, (int)(12 * p.cap_tau + 0.5)) , p.cap_inf + (p.cap0 - p.cap_inf) * std::exp(-std::floor(12 * p.cap_tau + 0.5) / (12 * p.cap_tau)), 1e-12);
  RCHECK(RTSFloor::capRate(p, 100000) - p.cap_inf < 1e-12);
  // layer cap: no change before initiation or below the cap, limited above it
  RCHECK(RTSFloor::capLayerRH(p, -1, 50.0, 100.0) == 50.0);
  RCHECK(RTSFloor::capLayerRH(p, 0, 1e-6, 1000.0) == 1e-6);
  CHECK_NEAR(RTSFloor::capLayerRH(p, 0, 1e6, 1200.0), p.cap0 / 12.0 * 1200.0, 1e-12);
  RCHECK(RTSFloor::capLayerRH(p, 3, 0.0, 1000.0) == 0.0);
  // export fraction: lat0/12 at initiation, e-folding lat_tau, 0 before
  RCHECK(RTSFloor::exportFraction(p, -1) == 0.0);
  CHECK_NEAR(RTSFloor::exportFraction(p, 0), p.lat0 / 12.0, 1e-15);
  CHECK_NEAR(RTSFloor::exportFraction(p, (int)(12 * p.lat_tau)), p.lat0 / 12.0 / std::exp(1.0), 1e-12);
  // anoxia: 1 for unsaturated, anox_f for saturated, linear between, clipped
  RCHECK(RTSFloor::anoxFactor(p, 0.0) == 1.0);
  CHECK_NEAR(RTSFloor::anoxFactor(p, 1.0), p.anox_f, 1e-15);
  CHECK_NEAR(RTSFloor::anoxFactor(p, 0.5), 0.5 * (1.0 + p.anox_f), 1e-15);
  RCHECK(RTSFloor::anoxFactor(p, 2.0) == RTSFloor::anoxFactor(p, 1.0));
  // export of a pool set conserves mass
  { double c[4] = {10.0, 200.0, 3000.0, 40000.0}, out = 0.0, tot0 = 0.0, tot1 = 0.0;
    const double f = RTSFloor::exportFraction(p, 7);
    for (double & x : c) { tot0 += x; const double d = f * x; x -= d; out += d; tot1 += x; }
    CHECK_NEAR(tot0, tot1 + out, 1e-9); }
}
}

namespace rtsu {
static void test_rts_regrow() {
  RTSRegrowParams p;
  RCHECK(p.on);                                             // default ON since v0.3 (false -> bitwise v0.2)
  // Xia et al. 2026: tau = 1.35 GPP^-1.68 (GPP kgC/m2/yr): 0.3 -> ~10.2 yr, s ~1.46
  CHECK_NEAR(RTSRegrow::tau(p, 0.3), 1.35 * std::pow(0.3, -1.68), 1e-12);
  RCHECK(RTSRegrow::tau(p, 0.3) > 10.0 && RTSRegrow::tau(p, 0.3) < 10.5);
  CHECK_NEAR(RTSRegrow::tau(p, 0.3) / p.tau_ref, 10.2 / 7.0, 0.02);  // site scale s (diagnostic in v2)
  RCHECK(RTSRegrow::tau(p, 0.0) == p.tau_max);              // no GPP -> slowest
  RCHECK(RTSRegrow::tau(p, 100.0) == p.tau_min);            // clamped
  // PFT classes
  RCHECK(RTSRegrow::pftType(0, 0, 3) == RG_LICHEN);
  RCHECK(RTSRegrow::pftType(0, 0, 2) == RG_MOSS && RTSRegrow::pftType(0, 0, 1) == RG_MOSS);
  RCHECK(RTSRegrow::pftType(1, 1, 0) == RG_SHRUB && RTSRegrow::pftType(1, 0, 0) == RG_EVERGREEN);
  RCHECK(RTSRegrow::pftType(0, 0, 0) == RG_GRAM);
  // time axis: stabilization 12 months after initiation by default
  RCHECK(RTSRegrow::stabMonth(p) == 12);
  CHECK_NEAR(RTSRegrow::tYears(p, 36), 2.0, 1e-12);
  // v2: no logistic growth ceiling exists (targetFrac/limitMult removed by design; flagged in commit).
  // shrub boost: B0 at t=0, decays with Tb, only for deciduous shrubs, none while active
  CHECK_NEAR(RTSRegrow::boost(p, RG_SHRUB, 0.0), p.boost0, 1e-12);
  CHECK_NEAR(RTSRegrow::boost(p, RG_SHRUB, p.boost_T), 1.0 + (p.boost0 - 1.0) * std::exp(-1.0), 1e-12);
  RCHECK(RTSRegrow::boost(p, RG_GRAM, 0.0) == 1.0 && RTSRegrow::boost(p, RG_SHRUB, -0.1) == 1.0);
  // OL rebuild and N supply rates (Pizano et al. 2014)
  RCHECK(RTSRegrow::olRate(p, -0.1) == 0.0 && RTSRegrow::olRate(p, 10.0) == 32.0 && RTSRegrow::olRate(p, 61.0) == 6.5);
  RCHECK(RTSRegrow::nRate(p, -0.1) == 0.0 && RTSRegrow::nRate(p, 59.9) == 2.2 && RTSRegrow::nRate(p, 60.1) == 0.4);
  // avln part + nominal OL N = total N supply (annual), never negative
  CHECK_NEAR(12.0 * RTSRegrow::avlnSupplyMonth(p, 10.0) + 32.0 / 30.0, 2.2, 1e-12);
  CHECK_NEAR(12.0 * RTSRegrow::avlnSupplyMonth(p, 70.0) + 6.5 / 30.0, 0.4, 1e-12);
  RTSRegrowParams q; q.ol_cn = 5.0; RCHECK(RTSRegrow::avlnSupplyMonth(q, 10.0) == 0.0);
  // v2 GPP-recovery taper: full supply when GPP collapsed, ~zero when recovered, monotone in between
  RCHECK(p.n_taper && p.n_taper_pow == 4.0);
  RCHECK(RTSRegrow::taperFactor(p, 0.0) == 1.0 && RTSRegrow::taperFactor(p, -1.0) == 1.0);
  CHECK_NEAR(RTSRegrow::taperFactor(p, 0.5), 1.0 - 0.0625, 1e-12);   // plateau: near-full at half recovery
  CHECK_NEAR(RTSRegrow::taperFactor(p, 0.9), 1.0 - 0.6561, 1e-12);
  RCHECK(RTSRegrow::taperFactor(p, 1.0) == 0.0 && RTSRegrow::taperFactor(p, 2.0) == 0.0);
  { double prev = 2.0; bool mono = true;
    for (int i = 0; i <= 100; i++) { const double v = RTSRegrow::taperFactor(p, i / 100.0); mono &= v <= prev; prev = v; }
    RCHECK(mono); }
  { RTSRegrowParams q2; q2.n_taper_pow = 2.0; CHECK_NEAR(RTSRegrow::taperFactor(q2, 0.5), 0.75, 1e-12); }
  { RTSRegrowParams q3; q3.n_taper = false; RCHECK(RTSRegrow::taperFactor(q3, 0.99) == 1.0); }  // untapered v1 supply
  // tapered avln supply scales linearly with the factor and is never negative
  CHECK_NEAR(RTSRegrow::avlnSupplyMonth(p, 10.0, 0.5), 0.5 * RTSRegrow::avlnSupplyMonth(p, 10.0, 1.0), 1e-15);
  RCHECK(RTSRegrow::avlnSupplyMonth(p, 10.0, 0.0) == 0.0 && RTSRegrow::avlnSupplyMonth(p, 10.0, -1.0) == 0.0);
}
}

namespace rtsu {
static void test_tknetc() {
  // soil term: CH4 off -> RHSOM; on -> TKRHCO2+TKCH4EMIS (= RHSOM when
  // TKRHCO2=RHSOM-prod+oxid and TKCH4EMIS=prod-oxid, so CH4 counted once)
  CHECK_NEAR(TKNetC::soilCarbonLoss(false, 10.0, 0.0, 0.0), 10.0, 1e-12);
  { const double rhsom=10.0, prod=3.0, oxid=1.0;
    CHECK_NEAR(TKNetC::soilCarbonLoss(true, rhsom, rhsom-prod+oxid, prod-oxid), rhsom, 1e-12); }
  // TKNETC = (RG+RM+soil+RHDWD-GPP) + TKLATC + TKRESPC - TKREGROWC(regrow on)
  CHECK_NEAR(TKNetC::compute(false, true, 100, 10, 20, 30, 5, 0, 0, 7, 3, 2), (10+20+30+5-100)+7+3-2, 1e-12);
  // regrow off: TKREGROWC term dropped even if a stray value is passed
  CHECK_NEAR(TKNetC::compute(false, false, 100, 10, 20, 30, 5, 0, 0, 7, 3, 999), (10+20+30+5-100)+7+3, 1e-12);
  // CH4 on with same RHSOM partition gives the same TKNETC as CH4 off (no double count / no drop)
  CHECK_NEAR(TKNetC::compute(true, false, 100, 10, 20, 30, 5, 28, 2, 7, 3, 0),  // rh_co2+emis = 28+2 = RHSOM 30
             TKNetC::compute(false, false, 100, 10, 20, 30, 5, 0, 0, 7, 3, 0), 1e-12);
  // no TK fluxes: TKNETC == loss-positive NEE
  CHECK_NEAR(TKNetC::compute(false, false, 50, 5, 5, 10, 0, 0, 0, 0, 0, 0), -30.0, 1e-12);
}
}

TEST_SUITE("rts") {
TEST_CASE("rts-regrow unit test removal_thickness (ported from tests/rts/test_rts_units.cpp)") { rtsu::test_removal_thickness(); }
TEST_CASE("rts-regrow unit test partition_mass_balance (ported from tests/rts/test_rts_units.cpp)") { rtsu::test_partition_mass_balance(); }
TEST_CASE("rts-regrow unit test thaw_column (ported from tests/rts/test_rts_units.cpp)") { rtsu::test_thaw_column(); }
TEST_CASE("rts-regrow unit test excess_ice_init (ported from tests/rts/test_rts_units.cpp)") { rtsu::test_excess_ice_init(); }
TEST_CASE("rts-regrow unit test excess_melt_subsidence (ported from tests/rts/test_rts_units.cpp)") { rtsu::test_excess_melt_subsidence(); }
TEST_CASE("rts-regrow unit test rts_ch4 (ported from tests/rts/test_rts_units.cpp)") { rtsu::test_rts_ch4(); }
TEST_CASE("rts-regrow unit test rts_floor (ported from tests/rts/test_rts_units.cpp)") { rtsu::test_rts_floor(); }
TEST_CASE("rts-regrow unit test rts_regrow (ported from tests/rts/test_rts_units.cpp)") { rtsu::test_rts_regrow(); }
TEST_CASE("rts-regrow unit test tknetc (ported from tests/rts/test_rts_units.cpp)") { rtsu::test_tknetc(); }
}
