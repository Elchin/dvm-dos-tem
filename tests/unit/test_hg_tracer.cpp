// Hg tracer + MeHg index (model_settings.hg_tracer / hg_mehg; include/HgTracer.h).
// Refs: Schaefer et al. 2020 doi:10.1038/s41467-020-18398-5 (16% evasion, 1.5% export of decay-released Hg);
//       Thompson et al. 2025 doi:10.1021/acs.est.5c04510 (%MeHg targets). See docs_src/WORK_NOTES_hg_tracer.md.
#include "doctest.h"
#include "HgTracer.h"
#include <json/json.h>
#include <cmath>
#include <random>
#include <fstream>
#include <cstdio>
#include <vector>

namespace {
Json::Value ms0() { Json::Value m(Json::objectValue); return m; }

hgtr::Snap column(int n, double c0 = 1000.0) {
  hgtr::Snap s; s.nsl = n;
  for (int l = 0; l < n; l++) { s.soil[l] = c0 + 100.0 * l; s.dz[l] = (l < 3) ? 0.05 : 0.1; s.type[l] = (l < 3) ? (l == 0 ? 1 : 2) : 3; }
  s.veg[0][0] = 50; s.veg[0][1] = 200; s.veg[0][2] = 80; s.veg[3][0] = 10; s.dead[0] = 30; s.wd = 20;
  return s;
}
double tot(const hgtr::State& s) { return s.total(); }
}

TEST_CASE("hg_tracer parse: absent = off, defaults, errors") {
  hgtr::Params p;
  CHECK(hgtr::parse(ms0(), p) == "");
  CHECK(p.tracer == 0); CHECK(p.mehg == 0);
  CHECK(p.f_evade == doctest::Approx(0.02)); CHECK(p.f_export == doctest::Approx(0.005)); CHECK(p.spinup_closure == 0); CHECK(p.k_meth_ch4 == 3250.0); CHECK(p.k_meth_proxy == 87.0); CHECK(p.f_leach == 0.0);
  CHECK(p.spinup == 1); CHECK(p.spinup_cycle == 10); CHECK(p.transport == 1); CHECK(p.log_kd_min == doctest::Approx(3.8)); CHECK(p.dep_file.empty()); CHECK(hgtr::dep_at(p, 2000, 6).rate == 0.0);
  CHECK(p.r_org == doctest::Approx(0.35)); CHECK(p.r_min == doctest::Approx(1.47));
  Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_r_min"] = 1.6; m["other_switch"] = 3;
  CHECK(hgtr::parse(m, p) == ""); CHECK(p.tracer == 1); CHECK(p.r_min == doctest::Approx(1.6));
  Json::Value b = ms0(); b["hg_tracer"] = 2; CHECK(hgtr::parse(b, p) != "");
  Json::Value c = ms0(); c["hg_trcer"] = 1; CHECK(hgtr::parse(c, p) != "");          // typo guard
  Json::Value d = ms0(); d["hg_mehg"] = 1; CHECK(hgtr::parse(d, p) != "");           // needs the tracer
  Json::Value e = ms0(); e["hg_tracer"] = 1; e["hg_f_evade"] = 0.9; e["hg_f_export"] = 0.2; CHECK(hgtr::parse(e, p) != "");
  Json::Value f = ms0(); f["hg_tracer"] = true; CHECK(hgtr::parse(f, p) != "");     // bool is not an int here
}

TEST_CASE("hg_tracer initialize: organic / mineral ratios by layer type, tissue ratios") {
  hgtr::Params p; hgtr::State s; hgtr::Snap c = column(6);
  hgtr::initialize(p, s, c);
  CHECK(s.init);
  CHECK(s.soil[0] == doctest::Approx(c.soil[0] * p.r_org));
  CHECK(s.soil[2] == doctest::Approx(c.soil[2] * p.r_org));
  CHECK(s.soil[4] == doctest::Approx(c.soil[4] * p.r_min));
  CHECK(s.soil[6] == 0.0);
  CHECK(s.veg[0][1] == doctest::Approx(200 * p.r_stem));
  CHECK(s.balance_error() == 0.0);
}

TEST_CASE("hg_tracer fixed-geometry BGC step: uptake, decay partition, litter transfer, exact balance") {
  hgtr::Params p; hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  const double h0 = tot(s), hsoil1 = s.soil[1];
  hgtr::Snap b = a; hgtr::Flux f;
  f.npp[0][0] = 12.0; f.npp[0][2] = 5.0; f.npp[3][0] = -1.0;   // negative NPP: C loss without Hg
  b.veg[0][0] = a.veg[0][0] + 12.0 - 8.0;                        // 8 g leaf litter
  b.veg[0][2] = a.veg[0][2] + 5.0 - 3.0;                         // 3 g root litter
  b.veg[3][0] = a.veg[3][0] - 1.0;
  f.decomp[1] = 20.0; f.decomp_wd = 2.0;
  b.soil[0] += 8.0; b.soil[1] += 3.0 - 20.0; b.wd -= 2.0;
  hgtr::step(p, s, a, b, f);
  const double up = 12.0 * p.r_leaf + 5.0 * p.r_root;
  CHECK(s.m.uptake == doctest::Approx(up));
  const double rel = 20.0 * hsoil1 / a.soil[1] + 2.0 * p.r_leaf;
  CHECK(s.m.evasion == doctest::Approx(rel * p.f_evade));
  CHECK(s.m.export_ == doctest::Approx(rel * p.f_export));
  CHECK(std::fabs(tot(s) - (h0 + up - s.m.evasion - s.m.export_)) < 1e-10);
  CHECK(std::fabs(s.balance_error()) < 1e-10);
  CHECK(s.m.unexpl_closs == doctest::Approx(0.0)); CHECK(s.m.unexpl_cgain == doctest::Approx(0.0));
  // negative NPP left the Hg in the PFT 3 leaves (no donor)
  CHECK(s.veg[3][0] == doctest::Approx(10 * p.r_leaf));
}

TEST_CASE("hg_tracer geometry change: split / combine keep Hg with C (uniform ratio) and conserve") {
  hgtr::Params p; hgtr::State s; hgtr::Snap a = column(5);
  for (int l = 0; l < 5; l++) { a.type[l] = 3; }
  hgtr::initialize(p, s, a);
  const double h0 = tot(s);
  hgtr::Snap b = a; b.nsl = 6;   // split layer 1 into two halves
  b.soil[1] = a.soil[1] / 2; b.dz[1] = a.dz[1] / 2;
  b.soil[2] = a.soil[1] / 2; b.dz[2] = a.dz[1] / 2;
  for (int l = 2; l < 5; l++) { b.soil[l + 1] = a.soil[l]; b.dz[l + 1] = a.dz[l]; }
  hgtr::step(p, s, a, b, hgtr::Flux());
  for (int l = 0; l < 6; l++) { CHECK(s.soil[l] == doctest::Approx(b.soil[l] * p.r_min)); }
  CHECK(std::fabs(tot(s) - h0) < 1e-10);
  hgtr::Snap c = b; c.nsl = 4;   // combine 0+1+2
  c.soil[0] = b.soil[0] + b.soil[1] + b.soil[2]; c.dz[0] = b.dz[0] + b.dz[1] + b.dz[2];
  for (int l = 3; l < 6; l++) { c.soil[l - 2] = b.soil[l]; c.dz[l - 2] = b.dz[l]; }
  for (int l = 4; l < MAX_SOI_LAY; l++) { c.soil[l] = 0; c.dz[l] = 0; }
  hgtr::step(p, s, b, c, hgtr::Flux());
  CHECK(s.nsl == 4);
  for (int l = 0; l < 4; l++) { CHECK(s.soil[l] == doctest::Approx(c.soil[l] * p.r_min)); }
  for (int l = 4; l < MAX_SOI_LAY; l++) { CHECK(s.soil[l] == 0.0); }
  CHECK(std::fabs(tot(s) - h0) < 1e-10);
}

TEST_CASE("hg_tracer fire: burned top layers and veg emit their Hg, killed veg -> standing dead") {
  hgtr::Params p; hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  const double h0 = tot(s), htop = s.soil[0];
  hgtr::Snap b = a; hgtr::Flux f;
  // whole top layer burned and removed, veg PFT 0: leaves burned, stem -> standing dead
  b.nsl = 5; for (int l = 0; l < 5; l++) { b.soil[l] = a.soil[l + 1]; b.dz[l] = a.dz[l + 1]; b.type[l] = a.type[l + 1]; }
  b.soil[5] = 0; b.dz[5] = 0;
  b.veg[0][0] = 0; b.veg[0][1] = 0; b.dead[0] = a.dead[0] + 200.0;
  f.burn = a.soil[0] + 50.0;
  hgtr::step(p, s, a, b, f);
  CHECK(std::fabs(tot(s) - (h0 - s.m.fire)) < 1e-10);
  CHECK(std::fabs(s.balance_error()) < 1e-10);
  // the pot holds the top layer, leaves and stems; fire takes its C share of the pot
  const double pot = htop + 50 * p.r_leaf + 200 * p.r_stem, potc = a.soil[0] + 250.0;
  CHECK(s.m.fire == doctest::Approx(pot * f.burn / potc));
  for (int l = 0; l < 5; l++) { CHECK(s.soil[l] == doctest::Approx(b.soil[l] * (b.type[l] >= 3 ? p.r_min : p.r_org))); }
}

TEST_CASE("hg_tracer random steps (geometry changes, sinks, unexplained C): balance to roundoff") {
  std::mt19937 rng(42); std::uniform_real_distribution<double> U(0.0, 1.0);
  hgtr::Params p; p.dep_soil = 12.0; hgtr::State s; hgtr::Snap a = column(8, 3000.0); hgtr::initialize(p, s, a);
  for (int it = 0; it < 3000; it++) {
    hgtr::Snap b = a; hgtr::Flux f;
    for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) {
      if (a.veg[i][j] <= 0) continue;
      f.npp[i][j] = 10 * U(rng) - 1; b.veg[i][j] = std::fmax(0.0, a.veg[i][j] + f.npp[i][j] - 8 * U(rng));
    }
    b.dead[0] = std::fmax(0.0, a.dead[0] + 3 * U(rng) - 2); b.wd = std::fmax(0.0, a.wd + U(rng) - 0.5);
    f.dep_soil = (U(rng) < 0.5) ? 1.0 : 0.0;
    const int kind = int(U(rng) * 3);
    if (kind == 0) {   // fixed geometry with decay
      for (int l = 0; l < a.nsl; l++) { f.decomp[l] = 5 * U(rng); b.soil[l] = std::fmax(0.0, a.soil[l] - f.decomp[l] + 4 * U(rng)); }
    } else {           // random re-division with net change and sinks
      int n = 3 + int(U(rng) * 15); b.nsl = n;
      double c0 = 0; for (int l = 0; l < a.nsl; l++) c0 += a.soil[l];
      double c1 = c0 * (0.97 + 0.05 * U(rng));
      double w[MAX_SOI_LAY], ws = 0; for (int l = 0; l < n; l++) { w[l] = 0.2 + U(rng); ws += w[l]; }
      for (int l = 0; l < MAX_SOI_LAY; l++) { b.soil[l] = (l < n) ? c1 * w[l] / ws : 0; b.dz[l] = (l < n) ? 0.05 + 0.1 * U(rng) : 0; b.type[l] = (l < 2) ? 1 : 3; }
      if (kind == 2) { f.burn = std::fmax(0.0, c0 - c1) * U(rng); f.lateral = 3 * U(rng); f.pot_decomp = 2 * U(rng); }
    }
    hgtr::step(p, s, a, b, f);
    a = b;
    if (it % 12 == 11) { s.m.clear(); }
  }
  CHECK(std::fabs(s.balance_error()) < 1e-8 * std::fmax(1.0, s.total()));
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s.soil[l] >= 0.0); }
  CHECK(s.cum.fire > 0.0); CHECK(s.cum.dep > 0.0);
}

TEST_CASE("hg_tracer restart pack / unpack round trip") {
  hgtr::Params p; hgtr::State s; hgtr::Snap a = column(7); hgtr::initialize(p, s, a);
  s.mehg[2] = 0.5; s.cum.evasion = 3.0; s.cum.uptake = 4.0;
  double v[hgtr::NSTATE]; hgtr::pack(s, v);
  hgtr::State t; CHECK(hgtr::unpack(t, v, 7));
  CHECK(t.total() == s.total()); CHECK(t.mehg[2] == 0.5); CHECK(t.balance_error() == s.balance_error());
  v[0] = -9999.0; hgtr::State u; CHECK_FALSE(hgtr::unpack(u, v, 7)); CHECK_FALSE(u.init);
}

TEST_CASE("hg_mehg: steady state %MeHg = k_meth A / k_demeth at t_ref, capped by Hg, no change to Hg") {
  hgtr::Params p; p.tracer = 1; p.mehg = 1; p.k_meth_ch4 = 30.0; p.k_meth_proxy = 0.6; // explicit, independent of calibrated defaults
  hgtr::State s; hgtr::Snap a = column(4); hgtr::initialize(p, s, a);
  const double h = s.total();
  double A[MAX_SOI_LAY] = {}, T[MAX_SOI_LAY];
  for (int l = 0; l < MAX_SOI_LAY; l++) { A[l] = 1.0e-3 * (l + 1); T[l] = p.t_ref; }
  for (int m = 0; m < 60; m++) { hgtr::mehg_month(p, s, A, T, 4, true); }
  for (int l = 0; l < 4; l++) { CHECK(s.mehg[l] / s.soil[l] == doctest::Approx(p.k_meth_ch4 * A[l] / p.k_demeth).epsilon(1e-6)); }
  CHECK(s.total() == h);
  // colder: demethylation slower by q10, production unchanged (q10_meth = 1) -> higher fraction
  for (int l = 0; l < MAX_SOI_LAY; l++) { T[l] = p.t_ref - 10.0; }
  for (int m = 0; m < 200; m++) { hgtr::mehg_month(p, s, A, T, 4, false); }
  CHECK(s.mehg[0] / s.soil[0] == doctest::Approx(p.k_meth_proxy * A[0] / (p.k_demeth / 2.0)).epsilon(1e-6));
  for (int l = 0; l < MAX_SOI_LAY; l++) { A[l] = 1e3; }
  hgtr::mehg_month(p, s, A, T, 4, true);
  for (int l = 0; l < 4; l++) { CHECK(s.mehg[l] <= s.soil[l]); }
}

// ---------------- T1b: deposition series, leaching, spin-up, restart v2
TEST_CASE("hg_dep_file: CSV parse, annual interpolation, monthly rows, veg_scale, errors") {
  const char* fn = "/tmp/hg_dep_unittest.csv";
  { std::ofstream f(fn); f << "# test\nyear,month,hg2_wet,hg2_dry,hg0_soil,veg_scale\n1850,0,0.1,0.5,0.4,0.4\n1950,0,0.2,1.5,1.3,0.8\n2000,7,1,1,1,2\n"; }
  hgtr::Params p; Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_dep_file"] = fn; m["hg_dep_soil"] = 9.0;
  REQUIRE(hgtr::parse(m, p) == "");
  CHECK(hgtr::dep_at(p, 1800, 1).rate == doctest::Approx(1.0));            // clamped to the first row
  CHECK(hgtr::dep_at(p, 1900, 1).rate == doctest::Approx(2.0));            // (1.0 + 3.0) / 2
  CHECK(hgtr::dep_at(p, 1900, 1).veg_scale == doctest::Approx(0.6));
  CHECK(hgtr::dep_at(p, 2100, 1).rate == doctest::Approx(3.0));
  CHECK(hgtr::dep_at(p, 2000, 7).rate == doctest::Approx(3.0)); CHECK(hgtr::dep_at(p, 2000, 7).veg_scale == doctest::Approx(2.0));
  { std::ofstream f(fn); f << "year,foo\n1850,1\n"; }
  CHECK(hgtr::parse(m, p) != "");                                           // no deposition column
  { std::ofstream f(fn); f << "year,dep\n1850,-1\n"; }
  CHECK(hgtr::parse(m, p) != "");                                           // negative
  Json::Value b = ms0(); b["hg_tracer"] = 1; b["hg_dep_file"] = "/nonexistent/x.csv"; CHECK(hgtr::parse(b, p) != "");
  Json::Value c = ms0(); c["hg_tracer"] = 1; c["hg_dep_file"] = 3; CHECK(hgtr::parse(c, p) != "");
  Json::Value d = ms0(); d["hg_spinup"] = 1; CHECK(hgtr::parse(d, p) == "");    // T1c: ignored without the tracer
  Json::Value e = ms0(); e["hg_tracer"] = 1; e["hg_f_leach"] = 0.99; CHECK(hgtr::parse(e, p) != "");
  std::remove(fn);
}

TEST_CASE("hg_f_leach: decay-released Hg moves one layer down (bottom layer -> export), exact balance; veg_scale") {
  hgtr::Params p; p.f_leach = 0.3; hgtr::State s; hgtr::Snap a = column(4); hgtr::initialize(p, s, a);
  const double h1 = s.soil[1], h2 = s.soil[2], h3 = s.soil[3];
  hgtr::Snap b = a; hgtr::Flux f; f.decomp[1] = 10.0; f.decomp[3] = 5.0; b.soil[1] -= 10.0; b.soil[3] -= 5.0;
  f.npp[0][0] = 4.0; f.veg_scale = 0.5; b.veg[0][0] += 4.0;
  hgtr::step(p, s, a, b, f);
  const double r1 = 10.0 * h1 / a.soil[1], r3 = 5.0 * h3 / a.soil[3];
  CHECK(s.soil[2] == doctest::Approx(h2 + 0.3 * r1));
  CHECK(s.soil[1] == doctest::Approx(h1 - r1 * (p.f_evade + p.f_export + 0.3)));
  CHECK(s.m.export_ == doctest::Approx((r1 + r3) * p.f_export + 0.3 * r3));
  CHECK(s.m.uptake == doctest::Approx(4.0 * p.r_leaf * 0.5));
  CHECK(std::fabs(s.balance_error()) < 1e-10);
}

namespace {
// one synthetic year (12 months): NPP, litterfall to layer 0, decomposition in all layers, s2d-like
// transfer 0 -> 1, deposition; C in steady state (b == a for the soil totals over the month).
std::vector<hgtr::Rec> synth_year(const hgtr::Snap& c) {
  std::vector<hgtr::Rec> v;
  for (int mth = 0; mth < 12; mth++) {
    hgtr::Rec r; r.a = c; r.b = c;
    const double g = (mth >= 5 && mth <= 8) ? 1.0 : 0.1;
    r.f.npp[0][0] = 6.0 * g; r.f.npp[0][2] = 3.0 * g; r.f.dep_soil = 0.6;
    r.b.veg[0][0] = c.veg[0][0]; r.b.veg[0][2] = c.veg[0][2];          // litterfall = NPP
    r.f.decomp[0] = 6.0 * g; r.f.decomp[1] = 3.0 * g; r.f.decomp[2] = 0.5 * g; r.f.decomp[3] = 0.05 * g;
    r.b.soil[0] = c.soil[0] + 9.0 * g - 6.0 * g - 1.0 * g;            // litter in, decay, 1 g to layer 1
    r.b.soil[1] = c.soil[1] + 1.0 * g - 3.0 * g + 2.0 * g;             // (2 g from layer 2, see below)
    r.b.soil[2] = c.soil[2] - 0.5 * g - 2.0 * g + 2.5 * g;             // cryoturbation-like mixing
    r.b.soil[3] = c.soil[3] - 0.05 * g - 0.45 * g;                     // ... unbalanced on purpose
    v.push_back(r);
    hgtr::Rec back; back.a = r.b; back.b = c; v.push_back(back);       // month-end C returns to c (gap step)
  }
  return v;
}
}

TEST_CASE("hg_spinup: affine jump == explicit repetition of the recorded cycle; booked as spin; conserves") {
  hgtr::Params p; p.tracer = 1; p.spinup = 1; p.spinup_cycle = 1; p.f_leach = 0.1; p.spinup_years = 37;
  hgtr::Snap c = column(5);
  const std::vector<hgtr::Rec> yr = synth_year(c);
  hgtr::State s; hgtr::initialize(p, s, c);
  for (const hgtr::Rec& r : yr) hgtr::step(p, s, r.a, r.b, r.f);                 // year 1 (recorded)
  hgtr::State e = s;                                                              // explicit: 37 more years
  for (int y = 0; y < 37; y++) for (const hgtr::Rec& r : yr) hgtr::step(p, e, r.a, r.b, r.f);
  hgtr::State j = s;
  REQUIRE(hgtr::spinup_jump(p, j, yr, c.nsl));
  for (int l = 0; l < 5; l++) CHECK(j.soil[l] == doctest::Approx(e.soil[l]).epsilon(1e-9));
  CHECK(j.vegtot() == doctest::Approx(e.vegtot()).epsilon(1e-9));
  CHECK(j.cum.spin == doctest::Approx(e.total() - s.total()).epsilon(1e-9));
  CHECK(std::fabs(j.balance_error()) < 1e-8 * j.total());
  // long jump: converges to the periodic steady state (one more cycle changes nothing)
  p.spinup_years = 200000; hgtr::State k = s; REQUIRE(hgtr::spinup_jump(p, k, yr, c.nsl));
  hgtr::State k2 = k; for (const hgtr::Rec& r : yr) hgtr::step(p, k2, r.a, r.b, r.f);
  for (int l = 0; l < 4; l++) CHECK(k2.soil[l] == doctest::Approx(k.soil[l]).epsilon(1e-6));
  CHECK(std::fabs(k2.total() - k.total()) < 1e-3 * std::fabs(k.total() - s.total()));   // pools without turnover may still grow linearly
  // state not at the end geometry of the record: refused, state unchanged
  hgtr::State q = s; q.nsl = c.nsl + 1; const double q0 = q.total(); CHECK_FALSE(hgtr::spinup_jump(p, q, yr, c.nsl)); CHECK(q.total() == q0);
}

TEST_CASE("hg_spinup closure: non-periodic record (C drift, layer re-division) == explicit (cycle + closure) repetition; conserves; no runaway") {
  hgtr::Params p; p.tracer = 1; p.spinup = 1; p.spinup_cycle = 1; p.spinup_years = 15; p.dep_soil = 0.0; p.spinup_closure = 1;
  hgtr::Snap c0 = column(5);
  // record: 12 months of decomposition + litter, soil C grows 2 % in layers 1-2 over the year, then a re-division 5 -> 6 layers
  std::vector<hgtr::Rec> yr; hgtr::Snap cur = c0;
  for (int m = 0; m < 12; m++) {
    hgtr::Rec r; r.a = cur; r.b = cur;
    for (int l = 0; l < 5; l++) { r.f.decomp[l] = 2.0 + l; }
    r.f.npp[0][0] = 3.0;
    r.b.soil[1] += 1.0; r.b.soil[2] += 0.7;                       // net C gain (not periodic)
    r.b.soil[0] = cur.soil[0] - (2.0) + 3.0;                       // litter in, decay out
    for (int l = 3; l < 5; l++) { r.b.soil[l] = cur.soil[l] - (2.0 + l); }
    yr.push_back(r); cur = r.b;
  }
  hgtr::Rec dv; dv.a = cur; dv.b = cur; dv.b.nsl = 6;               // re-division: layer 4 split in two
  dv.b.soil[4] = cur.soil[4] * 0.5; dv.b.soil[5] = cur.soil[4] * 0.5; dv.b.dz[4] = cur.dz[4] * 0.5; dv.b.dz[5] = cur.dz[4] * 0.5; dv.b.type[5] = cur.type[4];
  yr.push_back(dv);
  hgtr::State s; hgtr::initialize(p, s, c0);
  for (const hgtr::Rec& r : yr) hgtr::step(p, s, r.a, r.b, r.f);
  REQUIRE(s.nsl == 6);
  const hgtr::Flux none;
  hgtr::State e = s;                                                // explicit: closure, 14 x (cycle + closure), cycle
  hgtr::step(p, e, yr.back().b, yr.front().a, none);
  for (int y = 0; y < 14; y++) { for (const hgtr::Rec& r : yr) hgtr::step(p, e, r.a, r.b, r.f); hgtr::step(p, e, yr.back().b, yr.front().a, none); }
  for (const hgtr::Rec& r : yr) hgtr::step(p, e, r.a, r.b, r.f);
  hgtr::State j = s; REQUIRE(hgtr::spinup_jump(p, j, yr, 0));
  CHECK(j.nsl == 6);
  for (int l = 0; l < 6; l++) CHECK(j.soil[l] == doctest::Approx(e.soil[l]).epsilon(1e-9));
  CHECK(j.vegtot() == doctest::Approx(e.vegtot()).epsilon(1e-9));
  CHECK(std::fabs(j.balance_error()) < 1e-8 * j.total());
  CHECK(hgtr::rec_c_drift(yr) > 0.0);
  // long jump: bounded periodic steady state (the C gain of the record is not booked again in every repetition)
  p.spinup_years = 200000; hgtr::State k = s; REQUIRE(hgtr::spinup_jump(p, k, yr, 0));
  p.spinup_years = 100000; hgtr::State k2 = s; REQUIRE(hgtr::spinup_jump(p, k2, yr, 0));
  for (int l = 0; l < 6; l++) CHECK(k.soil[l] == doctest::Approx(k2.soil[l]).epsilon(1e-3));
  CHECK(std::isfinite(k.total()));
}

TEST_CASE("hg restart layout v4 round trip (spin, spun) and v1 still read") {
  hgtr::Params p; hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  s.cum.spin = 12.5; s.spun = 1; s.mehg[1] = 0.25;
  double v[hgtr::NSTATE]; hgtr::pack(s, v); CHECK(v[0] == 4.0);   // T2: layout v4 (+ snow, + diss)
  hgtr::State t; REQUIRE(hgtr::unpack(t, v, 6));
  CHECK(t.spun == 1); CHECK(t.cum.spin == 12.5); CHECK(t.balance_error() == doctest::Approx(s.balance_error()));
  v[0] = 1.0; hgtr::State u; REQUIRE(hgtr::unpack(u, v, 6)); CHECK(u.spun == 0); CHECK(u.cum.spin == 0.0);
}

TEST_CASE("hg deposition goes to the top soil layer that holds C (moss layer without SOM C skipped)") {
  hgtr::Params p; hgtr::State s; hgtr::Snap a = column(5); a.soil[0] = 0.0; hgtr::initialize(p, s, a);
  hgtr::Snap b = a; hgtr::Flux f; f.dep_soil = 0.5;
  const double h1 = s.soil[1];
  hgtr::step(p, s, a, b, f);
  CHECK(s.soil[0] == 0.0); CHECK(s.soil[1] == doctest::Approx(h1 + 0.5));
  CHECK(std::fabs(s.balance_error()) < 1e-12);
}

namespace {
// transport flux: infiltration qin, total lateral drainage D spread evenly over the n layers (mm/month)
hgtr::Flux tr_flux(int n, double qin, double D, double liq = 40.0) {
  hgtr::Flux f; f.tr = 1; f.q_in = qin; f.q_out = D;
  for (int l = 0; l < n; l++) { f.liq[l] = liq; f.ice[l] = 0.0; f.por[l] = 0.5; f.drain_l[l] = D / n; }
  return f;
}
}

TEST_CASE("hg_transport: dissolved Hg moves down with the water, deepest thawed layer -> drainage export, exact balance") {
  hgtr::Params p; hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  for (int l = 0; l < 6; l++) s.mehg[l] = 0.01 * s.soil[l];
  const hgtr::State s0 = s;
  hgtr::Flux f = tr_flux(6, 30.0, 10.0);
  hgtr::step(p, s, a, a, f);
  CHECK(s.soil[0] < s0.soil[0]);                                                   // top layer loses
  CHECK(s.m.drain > 0.0); CHECK(s.m.export_ == doctest::Approx(s.m.drain));       // only the drainage leaves
  CHECK(s.m.down > 0.0);                                                           // organic -> mineral (layer 2 -> 3)
  CHECK(s.soiltot() == doctest::Approx(s0.soiltot() - s.m.drain).epsilon(1e-13));
  CHECK(std::fabs(s.balance_error()) < 1e-12 * s.total());
  for (int l = 0; l < 6; l++) { CHECK(s.mehg[l] <= s.soil[l]); CHECK(s.mehg[l] == doctest::Approx(0.01 * s.soil[l]).epsilon(1e-9)); }
  // analytic check of the top layer: phi = (q_v + d_0) / (liq + Kd M), q_v = (qin - D)(1 - z/Z) + (D - d_0)
  const double Z = 0.05 * 3 + 0.1 * 3, d0 = 10.0 / 6, qv0 = (30.0 - 10.0) * (1.0 - 0.05 / Z) + (10.0 - d0);
  const double phi0 = (qv0 + d0) / (40.0 + std::pow(10.0, p.log_kd_org) * a.soil[0] * 2.0e-3);
  CHECK(s.soil[0] == doctest::Approx(s0.soil[0] * (1.0 - phi0)).epsilon(1e-12));
  // C unchanged by the transport sub-step: no uptake / evasion / unexplained C
  CHECK(s.m.uptake == 0.0); CHECK(s.m.evasion == 0.0); CHECK(s.m.unexpl_closs == 0.0); CHECK(s.m.unexpl_cgain == 0.0);
}

TEST_CASE("hg_transport: frozen layer stops the water; off switch; no water no change; works after a geometry change") {
  hgtr::Params p; hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  const hgtr::State s0 = s;
  hgtr::Flux f = tr_flux(6, 30.0, 5.0); f.liq[2] = 1.0; f.ice[2] = 80.0;          // layer 2 frozen
  hgtr::step(p, s, a, a, f);
  for (int l = 2; l < 6; l++) CHECK(s.soil[l] == s0.soil[l]);                      // nothing reaches below the frost
  CHECK(s.soil[1] != s0.soil[1]); CHECK(s.m.drain > 0.0); CHECK(s.m.down == 0.0);
  CHECK(std::fabs(s.balance_error()) < 1e-12 * s.total());
  hgtr::Params q = p; q.transport = 0; hgtr::State t = s0; hgtr::step(q, t, a, a, tr_flux(6, 30.0, 5.0));
  for (int l = 0; l < 6; l++) CHECK(t.soil[l] == s0.soil[l]);
  hgtr::State w = s0; hgtr::Flux fd = tr_flux(6, 30.0, 5.0); fd.liq[0] = 0.0;   // dry top layer passes the water
  hgtr::step(p, w, a, a, fd); CHECK(w.m.drain > 0.0); CHECK(w.soil[5] != s0.soil[5]);
  hgtr::State nd = s0; hgtr::step(p, nd, a, a, tr_flux(6, 30.0, 0.0));        // no drainage: Hg moves down, nothing leaves
  CHECK(nd.m.drain == 0.0); CHECK(nd.soiltot() == doctest::Approx(s0.soiltot()).epsilon(1e-14)); CHECK(nd.soil[5] > s0.soil[5]);
  hgtr::State u = s0; hgtr::step(p, u, a, a, tr_flux(6, 0.0, 0.0));
  for (int l = 0; l < 6; l++) CHECK(u.soil[l] == s0.soil[l]);
  // geometry change in the same sub-step (split of layer 4): transport applies on the new geometry, conserves
  hgtr::Snap b = a; b.nsl = 7; b.soil[6] = b.soil[5] * 0.5; b.soil[5] *= 0.5; b.dz[5] = 0.05; b.dz[6] = 0.05; b.type[6] = 3;
  hgtr::State v = s0; hgtr::step(p, v, a, b, tr_flux(7, 30.0, 10.0));
  CHECK(v.nsl == 7); CHECK(v.m.drain > 0.0); CHECK(std::fabs(v.balance_error()) < 1e-12 * v.total());
}

TEST_CASE("hg_transport: spin-up jump with transport sub-steps == explicit repetition") {
  hgtr::Params p; p.tracer = 1; p.dep_soil = 0.0; p.spinup_years = 23; p.spinup_cycle = 1;
  hgtr::Snap c = column(5);
  std::vector<hgtr::Rec> yr;
  for (int m = 0; m < 12; m++) {
    hgtr::Rec r; r.a = c; r.b = c;
    for (int l = 0; l < 5; l++) { r.f.decomp[l] = 2.0 + l; r.b.soil[l] = c.soil[l]; }
    r.f.npp[0][0] = 3.0; r.b.veg[0][0] = c.veg[0][0]; r.f.dep_soil = 0.4;
    for (int l = 0; l < 5; l++) { r.b.soil[l] += 0.0; }
    yr.push_back(r);
    hgtr::Rec t; t.a = c; t.b = c; t.f = tr_flux(5, (m > 4 && m < 9) ? 40.0 : 0.0, (m > 4 && m < 9) ? 15.0 : 0.0); yr.push_back(t);
  }
  hgtr::State s; hgtr::initialize(p, s, c);
  for (const hgtr::Rec& r : yr) hgtr::step(p, s, r.a, r.b, r.f);
  hgtr::State e = s;
  for (int y = 0; y < 23; y++) for (const hgtr::Rec& r : yr) hgtr::step(p, e, r.a, r.b, r.f);
  hgtr::State j = s; REQUIRE(hgtr::spinup_jump(p, j, yr, c.nsl));
  for (int l = 0; l < 5; l++) CHECK(j.soil[l] == doctest::Approx(e.soil[l]).epsilon(1e-9));
  CHECK(std::fabs(j.balance_error()) < 1e-8 * j.total());
}

TEST_CASE("hg_spinup / hg_transport defaults are ignored without the tracer") {
  hgtr::Params p; Json::Value m(Json::objectValue); m["hg_spinup"] = 1; m["hg_transport"] = 1;
  CHECK(hgtr::parse(m, p) == ""); CHECK(p.tracer == 0);
  Json::Value b(Json::objectValue); b["hg_log_kd_org"] = 12; CHECK(hgtr::parse(b, p) != "");
}

TEST_CASE("hg_spinup closure 0 (default): depth-overlap Hg mass map; identity for same geometry; == explicit; conserves") {
  hgtr::Params p; p.tracer = 1; p.spinup = 1; p.spinup_cycle = 1; p.spinup_years = 12; p.dep_soil = 0.0;
  CHECK(p.spinup_closure == 0);
  hgtr::Snap c0 = column(5);
  // depth map: same geometry = identity; split layer 4 in two halves -> its Hg halves; merge back -> original
  hgtr::State s; hgtr::initialize(p, s, c0);
  hgtr::State i1 = s; hgtr::depth_map_(i1, c0, c0);
  for (int l = 0; l < 5; l++) CHECK(i1.soil[l] == s.soil[l]);
  hgtr::Snap c1 = c0; c1.nsl = 6; c1.dz[4] = c0.dz[4] * 0.5; c1.dz[5] = c0.dz[4] * 0.5; c1.soil[4] = c0.soil[4] * 0.5; c1.soil[5] = c0.soil[4] * 0.5; c1.type[5] = c0.type[4];
  hgtr::State i2 = s; hgtr::depth_map_(i2, c0, c1);
  CHECK(i2.nsl == 6); CHECK(i2.soil[4] == doctest::Approx(0.5 * s.soil[4])); CHECK(i2.soil[5] == doctest::Approx(0.5 * s.soil[4]));
  CHECK(i2.total() == doctest::Approx(s.total()).epsilon(1e-12));
  hgtr::State i3 = i2; hgtr::depth_map_(i3, c1, c0);
  for (int l = 0; l < 5; l++) CHECK(i3.soil[l] == doctest::Approx(s.soil[l]).epsilon(1e-12));
  // organic horizon thins (layers 0-2 lose 40 % of their thickness): no Hg crosses the organic-mineral interface
  hgtr::Snap c2 = c0; for (int l = 0; l < 3; l++) c2.dz[l] = c0.dz[l] * 0.6;
  hgtr::State i4 = s; hgtr::depth_map_(i4, c0, c2);
  for (int l = 0; l < 5; l++) CHECK(i4.soil[l] == doctest::Approx(s.soil[l]).epsilon(1e-12));
  // organic layers merge 3 -> 2: organic Hg stays organic, mineral unchanged
  hgtr::Snap c3 = c0; c3.nsl = 4; c3.dz[0] = c0.dz[0]; c3.dz[1] = c0.dz[1] + c0.dz[2]; c3.type[1] = c0.type[1];
  for (int l = 2; l < 4; l++) { c3.dz[l] = c0.dz[l + 1]; c3.type[l] = c0.type[l + 1]; c3.soil[l] = c0.soil[l + 1]; }
  hgtr::State i5 = s; hgtr::depth_map_(i5, c0, c3);
  CHECK(i5.soil[0] + i5.soil[1] == doctest::Approx(s.soil[0] + s.soil[1] + s.soil[2]).epsilon(1e-12));
  CHECK(i5.soil[2] == doctest::Approx(s.soil[3])); CHECK(i5.soil[3] == doctest::Approx(s.soil[4]));
  // non-periodic record (C drift + re-division 5 -> 6): jump == explicit closure-0 repetition
  std::vector<hgtr::Rec> yr; hgtr::Snap cur = c0;
  for (int m = 0; m < 12; m++) {
    hgtr::Rec r; r.a = cur; r.b = cur;
    for (int l = 0; l < 5; l++) { r.f.decomp[l] = 2.0 + l; }
    r.f.npp[0][0] = 3.0; r.b.soil[1] += 1.0; r.b.soil[2] += 0.7; r.b.soil[0] = cur.soil[0] + 1.0;
    for (int l = 3; l < 5; l++) { r.b.soil[l] = cur.soil[l] - (2.0 + l); }
    yr.push_back(r); cur = r.b;
  }
  hgtr::Rec dv; dv.a = cur; dv.b = cur; dv.b.nsl = 6; dv.b.soil[4] = cur.soil[4] * 0.5; dv.b.soil[5] = cur.soil[4] * 0.5;
  dv.b.dz[4] = cur.dz[4] * 0.5; dv.b.dz[5] = cur.dz[4] * 0.5; dv.b.type[5] = cur.type[4]; yr.push_back(dv);
  hgtr::State s1 = s; for (const hgtr::Rec& r : yr) hgtr::step(p, s1, r.a, r.b, r.f);
  hgtr::State e = s1; hgtr::depth_map_(e, yr.back().b, yr.front().a);
  for (int y = 0; y < 11; y++) { for (const hgtr::Rec& r : yr) hgtr::step(p, e, r.a, r.b, r.f); hgtr::depth_map_(e, yr.back().b, yr.front().a); }
  for (const hgtr::Rec& r : yr) hgtr::step(p, e, r.a, r.b, r.f);
  hgtr::State j = s1; REQUIRE(hgtr::spinup_jump(p, j, yr, 0));
  CHECK(j.nsl == 6);
  for (int l = 0; l < 6; l++) CHECK(j.soil[l] == doctest::Approx(e.soil[l]).epsilon(1e-9));
  CHECK(std::fabs(j.balance_error()) < 1e-8 * j.total());
  // closure 2 (raw) refuses a record whose layer count differs at the two ends
  p.spinup_closure = 2; hgtr::State q = s1; CHECK_FALSE(hgtr::spinup_jump(p, q, yr, 0));
}

// ---------------- T2 phase 1 (hg_t2): surface interface ----------------
namespace {
const char* T2DEP = "/tmp/hg_t2_dep_unittest.csv";
hgtr::Params t2params(bool dry0col = true) {
  { std::ofstream f(T2DEP);
    if (dry0col) f << "year,hg2_wet,hg2_dry,hg0_soil,hg0_dry,veg_scale\n1850,1.2,2.4,1.11,3.0,1\n2000,1.2,2.4,1.11,3.0,1\n";
    else f << "year,hg2_wet,hg2_dry,hg0_soil,veg_scale\n1850,1.2,2.4,1.11,1\n"; }
  hgtr::Params p; Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_t2"] = 1; m["hg_dep_file"] = T2DEP; m["hg_t2_preset"] = 0;   // component defaults (phase-5 preset off)
  REQUIRE(hgtr::parse(m, p) == "");
  return p;
}
// a month's T2 forcing (deposition per month, LAI, snow)
hgtr::Flux t2flux(double lai0, double lai3, double swe, double melt, double sublim, double cosz) {
  hgtr::Flux f; f.t2 = 1; f.dep_wet = 0.1; f.dep_dry2 = 0.2; f.dep_dry0 = 0.25; f.dep_soil0 = 0.0925;
  f.lai[0] = lai0; f.lai[3] = lai3; f.vegcov[0] = 0.6; f.vegcov[3] = 0.3;
  f.swe = swe; f.melt = melt; f.sublim = sublim; f.sw = 150.0; f.cosz = cosz; f.ndays = 30.0;
  return f;
}
}

TEST_CASE("hg_t2 parse: defaults, needs the tracer and deposition components; hg0_dry column or hg0_soil / 0.37") {
  hgtr::Params p; REQUIRE(hgtr::parse(ms0(), p) == "");
  CHECK(p.t2 == 0); CHECK(p.t2_foliar == 1); CHECK(p.t2_snow == 1); CHECK(p.t2_wetprec == 1); CHECK(p.t2_photored == 0);
  CHECK(p.fol_fmax == 0.25); CHECK(p.fol_lai == 1.25); CHECK(p.fol_scale == 1.0); CHECK(p.snow_kred == 0.090);
  Json::Value a = ms0(); a["hg_t2"] = 1; CHECK(hgtr::parse(a, p) != "");                       // needs the tracer
  Json::Value b = ms0(); b["hg_tracer"] = 1; b["hg_t2"] = 1; CHECK(hgtr::parse(b, p) != "");   // needs components
  { std::ofstream f(T2DEP); f << "year,dep\n1850,3\n"; }
  Json::Value c = ms0(); c["hg_tracer"] = 1; c["hg_t2"] = 1; c["hg_dep_file"] = T2DEP; CHECK(hgtr::parse(c, p) != "");
  Json::Value d = ms0(); d["hg_tracer"] = 1; d["hg_fol_fmax"] = 1.5; CHECK(hgtr::parse(d, p) != "");
  hgtr::Params q = t2params(true);
  hgtr::DepVal dv = hgtr::dep_at(q, 1900, 3);
  CHECK(dv.rate == doctest::Approx(1.2 + 2.4 + 1.11));   // T1 rate unchanged: hg0_dry not added
  CHECK(dv.wet == doctest::Approx(1.2)); CHECK(dv.dry2 == doctest::Approx(2.4)); CHECK(dv.dry0 == doctest::Approx(3.0)); CHECK(dv.soil0 == doctest::Approx(1.11));
  hgtr::Params r = t2params(false);
  CHECK(hgtr::dep_at(r, 1900, 3).dry0 == doctest::Approx(1.11 / 0.37));
  std::remove(T2DEP);
}

TEST_CASE("hg_t2 f_fix(LAI) and 24-h mean cos zenith") {
  hgtr::Params p;
  CHECK(hgtr::f_fix(p, 0.0) == 0.0); CHECK(hgtr::f_fix(p, 0.625) == doctest::Approx(0.125));
  CHECK(hgtr::f_fix(p, 1.25) == doctest::Approx(0.25)); CHECK(hgtr::f_fix(p, 4.0) == doctest::Approx(0.25));
  // numerical 24-h integration of max(0, cos Z)
  auto num = [](double lat, double doy) {
    const double pi = 3.14159265358979323846, phi = lat * pi / 180, dl = 23.45 * pi / 180 * std::sin(2 * pi * (284 + doy) / 365);
    double sum = 0; const int n = 200000;
    for (int i = 0; i < n; i++) { const double h = -pi + 2 * pi * (i + 0.5) / n; sum += std::fmax(0.0, std::sin(phi) * std::sin(dl) + std::cos(phi) * std::cos(dl) * std::cos(h)); }
    return sum / n;
  };
  for (double lat : {0.0, 45.0, 68.6, 75.0, 85.0, -70.0}) for (double doy : {15.5, 105.0, 172.0, 288.0, 355.0}) {
    CHECK(hgtr::cosz24(lat, doy) == doctest::Approx(num(lat, doy)).epsilon(1e-5));
  }
  CHECK(hgtr::cosz24(80.0, 355.0) == 0.0);                 // polar night
  CHECK(hgtr::cosz24(80.0, 172.0) > 0.3);                  // polar day
  CHECK(hgtr::mid_doy(0) == 16.0); CHECK(hgtr::days_in(1) == 28.0);
}

TEST_CASE("hg_t2 off: identical (bit-exact) to the T1 step, also with T2 fields set in the flux") {
  std::mt19937 rng(7); std::uniform_real_distribution<double> U(0.0, 1.0);
  hgtr::Params p; p.dep_soil = 12.0; hgtr::State s1, s2; hgtr::Snap a = column(8, 3000.0); hgtr::initialize(p, s1, a); s2 = s1;
  hgtr::Params p2 = p; p2.t2 = 0;
  for (int it = 0; it < 600; it++) {
    hgtr::Snap b = a; hgtr::Flux f;
    for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) {
      if (a.veg[i][j] <= 0) continue;
      f.npp[i][j] = 10 * U(rng) - 1; b.veg[i][j] = std::fmax(0.0, a.veg[i][j] + f.npp[i][j] - 8 * U(rng));
    }
    for (int l = 0; l < a.nsl; l++) { f.decomp[l] = 5 * U(rng); b.soil[l] = std::fmax(0.0, a.soil[l] - f.decomp[l] + 4 * U(rng)); }
    f.dep_soil = 1.0;
    hgtr::Flux g = f; g.t2 = 1; g.lai[0] = 2.0; g.vegcov[0] = 0.5; g.swe = 50; g.melt = 10; g.dep_wet = 1; g.dep_dry0 = 1;   // ignored: p2.t2 = 0
    hgtr::step(p, s1, a, b, f); hgtr::step(p2, s2, a, b, g);
    a = b;
  }
  for (int l = 0; l < MAX_SOI_LAY; l++) CHECK(s1.soil[l] == s2.soil[l]);
  for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) CHECK(s1.veg[i][j] == s2.veg[i][j]);
  CHECK(s1.total() == s2.total()); CHECK(s2.snow == 0.0); CHECK(s1.cum.evasion == s2.cum.evasion);
}

TEST_CASE("hg_t2 foliar uptake: leaf Hg from LAI x canopy x dry deposition (snow-free), stem/root keep NPP x ratio; exact balance") {
  hgtr::Params p = t2params(); hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  const double h0 = s.total(), l00 = s.veg[0][0], l30 = s.veg[3][0];
  hgtr::Snap b = a; hgtr::Flux f = t2flux(2.0, 0.625, 0.0, 0.0, 0.0, 0.3);
  f.npp[0][0] = 10.0; f.npp[0][1] = 4.0; b.veg[0][0] += 10.0; b.veg[0][1] += 4.0;
  hgtr::step(p, s, a, b, f);
  const double u0 = 0.25 * 0.6 * 0.45, u3 = 0.125 * 0.3 * 0.45;   // f_fix x vegcov x (dry0 + dry2)
  CHECK(s.m.foliar == doctest::Approx(u0 + u3));
  CHECK(s.veg[0][0] == doctest::Approx(l00 + u0)); CHECK(s.veg[3][0] == doctest::Approx(l30 + u3));
  CHECK(s.m.uptake == doctest::Approx(u0 + u3 + 4.0 * p.r_stem));   // no leaf NPP x r_leaf
  const double fc = 0.25 * 0.6 + 0.125 * 0.3;
  CHECK(s.m.dep == doctest::Approx(0.1 + 0.2 * (1 - fc) + 0.0925));
  CHECK(s.m.snowdep == 0.0); CHECK(s.snow == 0.0);
  CHECK(std::fabs(s.total() - (h0 + s.m.uptake + s.m.dep - s.m.evasion - s.m.export_)) < 1e-14 * s.total());
  CHECK(std::fabs(s.balance_error()) < 1e-14 * s.total());
  // full snow cover: no foliar uptake; canopy cap: sum of fixing fractions <= 1
  hgtr::State s2; hgtr::initialize(p, s2, a); hgtr::Flux g = t2flux(2.0, 2.0, 50.0, 0.0, 0.0, 0.1);
  hgtr::step(p, s2, a, a, g); CHECK(s2.m.foliar == 0.0);
  hgtr::Params pc = p; pc.fol_scale = 10.0; hgtr::State s3; hgtr::initialize(pc, s3, a); hgtr::Flux h = t2flux(2.0, 2.0, 0.0, 0.0, 0.0, 0.1);
  hgtr::step(pc, s3, a, a, h); CHECK(s3.m.foliar == doctest::Approx(0.45)); CHECK(s3.m.dep == doctest::Approx(0.1 + 0.0925));
  CHECK(std::fabs(s3.balance_error()) < 1e-14 * s3.total());
}

TEST_CASE("hg_t2 snowpack: accumulation, photic photoreduction, melt / sublimation release, exact balance") {
  hgtr::Params p = t2params(); p.snow_sublim_loss = 1;  /* ORCHIDEE form: sublimation share lost */ hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  const double soil0 = s.soiltot();
  // winter month, deep snow, polar night: everything Hg(II) goes to the snow, nothing reduced or released
  hgtr::Flux w = t2flux(0.5, 0.0, 120.0, 0.0, 0.0, 0.0);
  hgtr::step(p, s, a, a, w);
  CHECK(s.snow == doctest::Approx(0.1 + 0.2)); CHECK(s.m.snowdep == doctest::Approx(0.3)); CHECK(s.m.foliar == 0.0);
  CHECK(s.soiltot() - soil0 == doctest::Approx(0.0925));   // only the soil Hg(0) sink
  // spring: photoreduction of the photic fraction 0.1 / (120 / 250) = 0.208 at k 0.09 /h x 24 x cos Z x 30 d
  s.m.clear(); const double sn = s.snow;
  hgtr::Flux sp = t2flux(0.5, 0.0, 120.0, 0.0, 0.0, 0.2);
  hgtr::step(p, s, a, a, sp);
  const double fph = 0.1 / (120.0 / 250.0), red = (sn + 0.3) * fph * (1 - std::exp(-0.09 * 24 * 0.2 * 30));
  CHECK(s.m.snowred == doctest::Approx(red)); CHECK(s.snow == doctest::Approx(sn + 0.3 - red));
  // melt month: release = (melt + sublim) / (melt + sublim + swe); melt share to soil, sublimation share to air
  s.m.clear(); const double sn2 = s.snow, so2 = s.soiltot();
  hgtr::Flux m = t2flux(0.5, 0.0, 40.0, 60.0, 20.0, 0.0);
  hgtr::step(p, s, a, a, m);
  const double fs = 1.0, add = 0.3 * fs, tot = sn2 + add, rel = 80.0 / 120.0;
  CHECK(s.m.snowmelt == doctest::Approx(tot * rel * 0.75)); CHECK(s.m.snowred == doctest::Approx(tot * rel * 0.25));
  CHECK(s.snow == doctest::Approx(tot * (1 - rel)));
  CHECK(s.soiltot() - so2 == doctest::Approx(s.m.snowmelt + 0.0925));
  // snow gone: everything left is released
  s.m.clear(); hgtr::Flux z = t2flux(0.5, 0.0, 0.0, 5.0, 0.0, 0.3); hgtr::step(p, s, a, a, z); CHECK(s.snow == 0.0);
  CHECK(std::fabs(s.balance_error()) < 1e-14 * s.total());
}

TEST_CASE("hg_t2 photoreduction of ground Hg(II) (hg_t2_photored)") {
  hgtr::Params p = t2params(); p.t2_photored = 1; hgtr::State s; hgtr::Snap a = column(5); hgtr::initialize(p, s, a);
  hgtr::Flux f = t2flux(0.0, 0.0, 0.0, 0.0, 0.0, 0.3); f.sw = 100.0;
  hgtr::step(p, s, a, a, f);
  const double fpr = 0.67 * std::exp(-0.016 * 100.0);
  CHECK(s.m.photored == doctest::Approx(0.3 * fpr)); CHECK(s.m.evasion == doctest::Approx(0.3 * fpr));
  CHECK(std::fabs(s.balance_error()) < 1e-14 * s.total());
}

TEST_CASE("hg_t2 random steps (geometry changes, fire, snow seasons): balance to roundoff, pools >= 0") {
  std::mt19937 rng(11); std::uniform_real_distribution<double> U(0.0, 1.0);
  hgtr::Params p = t2params(); p.t2_photored = 1; hgtr::State s; hgtr::Snap a = column(8, 3000.0); hgtr::initialize(p, s, a);
  for (int it = 0; it < 3000; it++) {
    const int mth = it % 12;
    hgtr::Snap b = a; hgtr::Flux f = t2flux(3 * U(rng), U(rng), (mth < 5 || mth > 9) ? 200 * U(rng) : 0.0, 30 * U(rng), 5 * U(rng), 0.4 * U(rng));
    if (U(rng) < 0.3) { f.t2 = 0; f.dep_soil = 1.0; }   // non-BGC sub-steps
    for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) {
      if (a.veg[i][j] <= 0) continue;
      f.npp[i][j] = 10 * U(rng) - 1; b.veg[i][j] = std::fmax(0.0, a.veg[i][j] + f.npp[i][j] - 8 * U(rng));
    }
    const int kind = int(U(rng) * 3);
    if (kind == 0) {
      for (int l = 0; l < a.nsl; l++) { f.decomp[l] = 5 * U(rng); b.soil[l] = std::fmax(0.0, a.soil[l] - f.decomp[l] + 4 * U(rng)); }
    } else {
      int n = 3 + int(U(rng) * 15); b.nsl = n;
      double c0 = 0; for (int l = 0; l < a.nsl; l++) c0 += a.soil[l];
      double c1 = c0 * (0.97 + 0.05 * U(rng));
      double w[MAX_SOI_LAY], ws = 0; for (int l = 0; l < n; l++) { w[l] = 0.2 + U(rng); ws += w[l]; }
      for (int l = 0; l < MAX_SOI_LAY; l++) { b.soil[l] = (l < n) ? c1 * w[l] / ws : 0; b.dz[l] = (l < n) ? 0.05 + 0.1 * U(rng) : 0; b.type[l] = (l < 2) ? 1 : 3; }
      if (kind == 2) { f.burn = std::fmax(0.0, c0 - c1) * U(rng); f.lateral = 3 * U(rng); f.pot_decomp = 2 * U(rng); }
    }
    hgtr::step(p, s, a, b, f);
    a = b;
    if (it % 12 == 11) { s.m.clear(); }
  }
  CHECK(std::fabs(s.balance_error()) < 1e-8 * std::fmax(1.0, s.total()));
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s.soil[l] >= 0.0); }
  CHECK(s.snow >= 0.0); CHECK(s.cum.foliar > 0.0); CHECK(s.cum.snowred > 0.0); CHECK(s.cum.snowmelt > 0.0);
  std::remove(T2DEP);
}

TEST_CASE("hg_t2 spin-up jump with the snowpack == explicit repetition; restart v3 keeps the snow, v2 reads as 0") {
  hgtr::Params p = t2params(); p.spinup = 1; p.spinup_cycle = 1; p.spinup_years = 23;
  hgtr::Snap c = column(5);
  std::vector<hgtr::Rec> yr = synth_year(c);
  for (size_t k = 0; k < yr.size(); k += 2) {   // BGC records: T2 surface forcing with a snow season
    const int mth = (int)(k / 2);
    const bool snow = (mth < 4 || mth > 9);
    yr[k].f = [&]() { hgtr::Flux f = yr[k].f; hgtr::Flux g = t2flux(mth >= 5 && mth <= 8 ? 1.0 : 0.2, 0.5, snow ? 80.0 : 0.0,
                       mth == 4 ? 80.0 : 0.0, snow ? 2.0 : 0.0, 0.05 + 0.03 * mth);
                       for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) g.npp[i][j] = f.npp[i][j];
                       for (int l = 0; l < MAX_SOI_LAY; l++) g.decomp[l] = f.decomp[l];
                       return g; }();
  }
  hgtr::State s; hgtr::initialize(p, s, c);
  for (const hgtr::Rec& r : yr) hgtr::step(p, s, r.a, r.b, r.f);
  hgtr::State e = s; for (int y = 0; y < 23; y++) for (const hgtr::Rec& r : yr) hgtr::step(p, e, r.a, r.b, r.f);
  hgtr::State j = s; REQUIRE(hgtr::spinup_jump(p, j, yr, c.nsl));
  for (int l = 0; l < 5; l++) CHECK(j.soil[l] == doctest::Approx(e.soil[l]).epsilon(1e-9));
  CHECK(j.snow == doctest::Approx(e.snow).epsilon(1e-9)); CHECK(j.vegtot() == doctest::Approx(e.vegtot()).epsilon(1e-9));
  CHECK(std::fabs(j.balance_error()) < 1e-8 * j.total());
  // restart
  s.snow = 0.75; double v[hgtr::NSTATE]; hgtr::pack(s, v); CHECK(v[0] == 4.0); CHECK(v[hgtr::NSTATE - 1 - MAX_SOI_LAY] == 0.75);
  hgtr::State t; REQUIRE(hgtr::unpack(t, v, 5)); CHECK(t.snow == 0.75); CHECK(t.total() == s.total());
  v[0] = 2.0; hgtr::State u; REQUIRE(hgtr::unpack(u, v, 5)); CHECK(u.snow == 0.0);
  std::remove(T2DEP);
}

// ---------------- T2 phase 2 (hg_t2_diss): dissolved pool, daily export, mixing ----------------
namespace {
hgtr::Params t2dparams() { hgtr::Params p = t2params(); p.t2_diss = 1; return p; }
// transport sub-step with nd days of constant water
hgtr::Flux dflux(int nl, int nd, double qin, double qover, double drain_l2, double liq, int frozen_from = 99) {
  hgtr::Flux f; f.tr = 1; f.nd = nd;
  for (int l = 0; l < nl; l++) { f.liq[l] = liq; f.ice[l] = (l >= frozen_from) ? 100.0 : 0.0; f.por[l] = 0.5; }
  for (int d = 0; d < nd; d++) {
    f.dq_in[d] = qin; f.dq_over[d] = qover;
    for (int l = 0; l < nl; l++) { f.d_liq[d][l] = liq; f.d_ice[d][l] = (l >= frozen_from) ? 100.0 : 0.0; }
    if (nl > 2) f.d_drain[d][2] = drain_l2;
  }
  return f;
}
}

TEST_CASE("hg_t2_diss parse: defaults, needs hg_t2; D7 default keeps sublimated snow Hg in the pack") {
  hgtr::Params p; REQUIRE(hgtr::parse(ms0(), p) == "");
  CHECK(p.t2_diss == 0); CHECK(p.t2_log_kd_org == 3.5); CHECK(p.tau_sorb == 0.25); CHECK(p.rho_moss == 30.0); CHECK(p.v_adv == 0.2); CHECK(p.z_runoff == 0.045); CHECK(p.d_cryo == 1e-3); CHECK(p.d_bio == 1e-4);
  CHECK(p.snow_sublim_loss == 0);
  Json::Value a = ms0(); a["hg_tracer"] = 1; a["hg_t2_diss"] = 1; CHECK(hgtr::parse(a, p) != "");
  hgtr::Params q = t2params(); hgtr::State s; hgtr::Snap c = column(5); hgtr::initialize(q, s, c);
  s.snow = 1.0; s.h_init += 1.0; const double so = s.soiltot();
  hgtr::Flux m = t2flux(0.5, 0.0, 40.0, 60.0, 20.0, 0.0); m.dep_wet = 0; m.dep_dry2 = 0; m.dep_dry0 = 0; m.dep_soil0 = 0;
  hgtr::step(q, s, c, c, m);
  CHECK(s.m.snowmelt == doctest::Approx(60.0 / 100.0)); CHECK(s.m.snowred == 0.0); CHECK(s.snow == doctest::Approx(0.4));
  CHECK(s.soiltot() - so == doctest::Approx(0.6)); CHECK(std::fabs(s.balance_error()) < 1e-14 * s.total());
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_diss: decay-released Hg goes to the dissolved pool, no f_export; geometry change re-sorbs; exact balance") {
  hgtr::Params p = t2dparams(); hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  const double h1 = s.soil[1];
  hgtr::Snap b = a; hgtr::Flux f; f.decomp[1] = 20.0; b.soil[1] -= 20.0;
  hgtr::step(p, s, a, b, f);
  const double rel = 20.0 * h1 / a.soil[1];
  CHECK(s.m.export_ == 0.0); CHECK(s.m.evasion == doctest::Approx(rel * p.f_evade));
  CHECK(s.diss[1] == doctest::Approx(rel * (1 - p.f_evade)));
  CHECK(std::fabs(s.balance_error()) < 1e-14 * s.total());
  hgtr::Snap g = b; g.nsl = 7; for (int l = 0; l < 7; l++) { g.soil[l] = (l < 6) ? b.soil[l] * 6.0 / 7.0 : 0.0; g.dz[l] = 0.06; g.type[l] = (l < 3) ? 1 : 3; }
  double c0 = 0, c1 = 0; for (int l = 0; l < 6; l++) c0 += b.soil[l]; for (int l = 0; l < 7; l++) c1 += g.soil[l]; g.soil[0] += c0 - c1;
  hgtr::step(p, s, b, g, hgtr::Flux());
  CHECK(s.disstot() == 0.0); CHECK(std::fabs(s.balance_error()) < 1e-13 * s.total());
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_diss sorption: no water flow -> dissolved fraction relaxes to liq / (liq + Kd M); no export") {
  hgtr::Params p = t2dparams(); p.d_cryo = 0; p.d_bio = 0; hgtr::State s; hgtr::Snap a = column(5); hgtr::initialize(p, s, a);
  const double h0 = s.total();
  for (int m = 0; m < 6; m++) { hgtr::Flux f = dflux(5, 30, 0.0, 0.0, 0.0, 30.0); hgtr::step(p, s, a, a, f); }
  const double kdo = std::pow(10.0, p.t2_log_kd_org), kdm = std::pow(10.0, p.log_kd_min);
  for (int l = 0; l < 5; l++) {
    const bool mn = hgtr::is_mineral(a.type[l]);
    const double M = mn ? 0.5 * a.dz[l] * 2650.0 : a.soil[l] * 2e-3, phi = 30.0 / (30.0 + (mn ? kdm : kdo) * M);
    CHECK(s.diss[l] / (s.diss[l] + s.soil[l]) == doctest::Approx(phi).epsilon(1e-6));
  }
  CHECK(s.cum.export_ == 0.0); CHECK(std::fabs(s.total() - h0) < 1e-13 * h0);
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_diss transport: overland from the top zone, lateral drainage, vertical advection; frozen layer stops it; exact balance") {
  hgtr::Params p = t2dparams(); p.d_cryo = 0; p.d_bio = 0; hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  for (int l = 0; l < 6; l++) { s.diss[l] = 1.0; }
  const double h0 = s.total();
  hgtr::Flux f = dflux(6, 1, 10.0, 5.0, 2.0, 20.0, 4); p.tau_sorb = 1e5;   // ~no sorption in 1 day
  hgtr::step(p, s, a, a, f);
  CHECK(s.m.runoff > 0.0); CHECK(s.m.drain > 0.0);
  CHECK(s.m.export_ == doctest::Approx(s.m.runoff + s.m.drain));
  CHECK(s.diss[4] == doctest::Approx(1.0).epsilon(1e-4)); CHECK(s.diss[5] == doctest::Approx(1.0).epsilon(1e-4));   // frozen: untouched (sorption ~0)
  CHECK(s.diss[3] > 1.0);   // receives from above, no outflow into the frozen layer
  CHECK(std::fabs(s.total() - (h0 - s.m.export_)) < 1e-13 * h0);
  // all frozen: nothing moves
  hgtr::State t = s; t.m.clear(); hgtr::Flux z = dflux(6, 3, 10.0, 5.0, 2.0, 20.0, 0); hgtr::step(p, t, a, a, z); CHECK(t.m.export_ == 0.0);
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_diss moss sorbent (dz x rho_moss) and joint removal fractions capped at 0.9") {
  hgtr::Params p = t2dparams(); p.d_cryo = 0; p.d_bio = 0; hgtr::State s; hgtr::Snap a = column(5);
  a.type[0] = 0; a.soil[0] = 0.0; a.dz[0] = 0.02; hgtr::initialize(p, s, a);
  s.soil[0] = 5.0; s.h_init = s.total(); s.cum.clear();
  const double h0 = s.total();
  for (int m = 0; m < 3; m++) { hgtr::Flux f = dflux(5, 30, 0.0, 0.0, 0.0, 30.0); hgtr::step(p, s, a, a, f); }
  const double phi = 30.0 / (30.0 + std::pow(10.0, p.t2_log_kd_org) * 0.02 * p.rho_moss);
  CHECK(s.diss[0] / (s.diss[0] + s.soil[0]) == doctest::Approx(phi).epsilon(1e-6));
  CHECK(std::fabs(s.total() - h0) < 1e-13 * h0);
  // huge overland + drainage + infiltration, no sorption: one day removes at most 90 % of each layer's diss
  hgtr::State t; hgtr::Snap c = column(5); hgtr::initialize(p, t, c); p.tau_sorb = 1e5;
  for (int l = 0; l < 5; l++) { t.diss[l] = 1.0; }
  const double t0 = t.total();
  hgtr::Flux g = dflux(5, 1, 1.0e4, 1.0e4, 1.0e3, 1.0, 99); hgtr::step(p, t, c, c, g);
  for (int l = 0; l < 5; l++) { CHECK(t.diss[l] >= 0.1 - 1e-4); }   // (tau 1e5 d: ~1e-5 sorption)
  CHECK(std::fabs(t.total() - (t0 - t.m.export_)) < 1e-13 * t0);
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_diss mixing: uniform concentration -> no flux; gradient -> downward, conservative, pools >= 0") {
  hgtr::Params p = t2dparams(); p.d_cryo = 1.0; hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  for (int l = 0; l < 6; l++) { s.soil[l] = 10.0 * a.dz[l]; s.diss[l] = 0.0; }
  s.h_init = s.total(); s.cum.clear();
  hgtr::Flux f = dflux(6, 0, 0, 0, 0, 20.0);
  hgtr::State u = s; hgtr::step(p, u, a, a, f);
  for (int l = 0; l < 6; l++) CHECK(u.soil[l] == doctest::Approx(s.soil[l]).epsilon(1e-12));
  s.soil[0] = 5.0; s.h_init = s.total();
  const double h0 = s.total(); hgtr::step(p, s, a, a, f);
  CHECK(s.soil[0] < 5.0); CHECK(s.soil[1] > 0.5); CHECK(std::fabs(s.total() - h0) < 1e-13 * h0);
  for (int l = 0; l < 6; l++) CHECK(s.soil[l] >= 0.0);
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_diss random steps (geometry changes, fire, snow, daily transport): balance to roundoff, pools >= 0") {
  std::mt19937 rng(23); std::uniform_real_distribution<double> U(0.0, 1.0);
  hgtr::Params p = t2dparams(); hgtr::State s; hgtr::Snap a = column(8, 3000.0); hgtr::initialize(p, s, a);
  for (int it = 0; it < 2000; it++) {
    const int mth = it % 12;
    hgtr::Snap b = a; hgtr::Flux f = t2flux(3 * U(rng), U(rng), (mth < 5 || mth > 9) ? 200 * U(rng) : 0.0, 30 * U(rng), 5 * U(rng), 0.4 * U(rng));
    for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) {
      if (a.veg[i][j] <= 0) continue;
      f.npp[i][j] = 10 * U(rng) - 1; b.veg[i][j] = std::fmax(0.0, a.veg[i][j] + f.npp[i][j] - 8 * U(rng));
    }
    if (U(rng) < 0.6) {
      for (int l = 0; l < a.nsl; l++) { f.decomp[l] = 5 * U(rng); b.soil[l] = std::fmax(0.0, a.soil[l] - f.decomp[l] + 4 * U(rng)); }
    } else {
      int n = 3 + int(U(rng) * 15); b.nsl = n;
      double c0 = 0; for (int l = 0; l < a.nsl; l++) c0 += a.soil[l];
      double c1 = c0 * (0.97 + 0.05 * U(rng));
      double w[MAX_SOI_LAY], ws = 0; for (int l = 0; l < n; l++) { w[l] = 0.2 + U(rng); ws += w[l]; }
      for (int l = 0; l < MAX_SOI_LAY; l++) { b.soil[l] = (l < n) ? c1 * w[l] / ws : 0; b.dz[l] = (l < n) ? 0.05 + 0.1 * U(rng) : 0; b.type[l] = (l < 2) ? 1 : 3; }
      if (U(rng) < 0.5) { f.burn = std::fmax(0.0, c0 - c1) * U(rng); f.lateral = 3 * U(rng); f.pot_decomp = 2 * U(rng); }
    }
    hgtr::step(p, s, a, b, f);
    hgtr::Flux t = dflux(b.nsl, 30, 20 * U(rng), 10 * U(rng), 3 * U(rng), 5 + 40 * U(rng), (mth < 5 || mth > 9) ? int(U(rng) * 3) : 99);
    hgtr::step(p, s, b, b, t);
    a = b;
    if (it % 12 == 11) { s.m.clear(); }
  }
  CHECK(std::fabs(s.balance_error()) < 1e-8 * std::fmax(1.0, s.total()));
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s.soil[l] >= 0.0); CHECK(s.diss[l] >= 0.0); }
  CHECK(s.cum.export_ > 0.0);
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_diss spin-up jump == explicit repetition (dissolved pool in the affine vector); restart v4 / v3") {
  hgtr::Params p = t2dparams(); p.spinup = 1; p.spinup_cycle = 1; p.spinup_years = 19;
  hgtr::Snap c = column(5);
  std::vector<hgtr::Rec> yr = synth_year(c);
  std::vector<hgtr::Rec> rec;
  for (size_t k = 0; k < yr.size(); k++) {
    rec.push_back(yr[k]);
    if (k % 2 == 1) {   // after each month: a transport sub-step with daily water
      const int mth = (int)(k / 2);
      hgtr::Rec t; t.a = c; t.b = c; t.f = dflux(5, 30, (mth >= 4 && mth <= 8) ? 3.0 : 0.0, mth == 4 ? 4.0 : 0.0, 0.5, 25.0, (mth < 4 || mth > 9) ? 1 : 99);
      rec.push_back(t);
    }
  }
  hgtr::State s; hgtr::initialize(p, s, c);
  for (const hgtr::Rec& r : rec) hgtr::step(p, s, r.a, r.b, r.f);
  hgtr::State e = s; for (int y = 0; y < 19; y++) for (const hgtr::Rec& r : rec) hgtr::step(p, e, r.a, r.b, r.f);
  hgtr::State j = s; REQUIRE(hgtr::spinup_jump(p, j, rec, c.nsl));
  for (int l = 0; l < 5; l++) { CHECK(j.soil[l] == doctest::Approx(e.soil[l]).epsilon(1e-9)); CHECK(j.diss[l] == doctest::Approx(e.diss[l]).epsilon(1e-9)); }
  CHECK(std::fabs(j.balance_error()) < 1e-8 * j.total());
  double v[hgtr::NSTATE]; hgtr::pack(e, v); CHECK(v[0] == 4.0);
  hgtr::State t; REQUIRE(hgtr::unpack(t, v, 5)); for (int l = 0; l < 5; l++) CHECK(t.diss[l] == e.diss[l]); CHECK(t.total() == e.total());
  v[0] = 3.0; hgtr::State u; REQUIRE(hgtr::unpack(u, v, 5)); CHECK(u.disstot() == 0.0); CHECK(u.snow == e.snow);
  std::remove(T2DEP);
}

// ---- T2 phase 3: process soil re-emission, ground Hg(0) sink, stomatal (GPP) foliar uptake, daily melt routing
TEST_CASE("hg_t2 phase 3 parse: defaults off, ranges, dependencies") {
  hgtr::Params p = t2params();
  CHECK(p.t2_evasion == 0); CHECK(p.t2_meltroute == 0); CHECK(p.t2_foliar == 1);
  CHECK(p.fd == 0.053); CHECK(p.ev_z == 0.05); CHECK(p.v_ground == 8.5e-5); CHECK(p.gem == 1.3); CHECK(p.gs_m == 4.0); CHECK(p.gs_eff == 0.75); CHECK(p.gs_ingpp == 1);
  Json::Value a = ms0(); a["hg_tracer"] = 1; a["hg_t2_evasion"] = 1; CHECK(hgtr::parse(a, p) != "");          // needs hg_t2
  Json::Value b = ms0(); b["hg_tracer"] = 1; b["hg_t2"] = 1; b["hg_dep_file"] = T2DEP; b["hg_t2_preset"] = 0; b["hg_t2_meltroute"] = 1;
  CHECK(hgtr::parse(b, p) != "");                                                                          // needs hg_t2_diss
  b["hg_t2_diss"] = 1; CHECK(hgtr::parse(b, p) == ""); CHECK(p.t2_meltroute == 1);
  Json::Value c = ms0(); c["hg_tracer"] = 1; c["hg_t2"] = 1; c["hg_dep_file"] = T2DEP; c["hg_t2_preset"] = 0; c["hg_t2_foliar"] = 2; c["hg_gs_eff"] = 0.5;
  CHECK(hgtr::parse(c, p) == ""); CHECK(p.t2_foliar == 2); CHECK(p.gs_eff == 0.5);
  c["hg_t2_foliar"] = 3; CHECK(hgtr::parse(c, p) != "");
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_evasion: re-emission fd x (1 - f_snow) x [T > 0] x exp(-z / ev_z) of decay-released Hg; ground sink v x GEM x s") {
  hgtr::Params p = t2params(); p.t2_evasion = 1; p.t2_foliar = 0;
  hgtr::State s; hgtr::Snap a = column(4); hgtr::initialize(p, s, a);
  hgtr::Snap b = a; hgtr::Flux f = t2flux(0, 0, 0, 0, 0, 0.3);
  f.decomp[0] = 10.0; f.decomp[3] = 10.0; b.soil[0] -= 10.0; b.soil[3] -= 10.0;
  f.ts[0] = 5.0; f.ts[3] = 5.0; f.ta = 5.0;
  const double h0 = s.soil[0], h3 = s.soil[3];
  hgtr::step(p, s, a, b, f);
  const double r0 = 10.0 * h0 / a.soil[0], r3 = 10.0 * h3 / a.soil[3];
  const double z3 = 0.05 * 3 + 0.05;   // layer 3 (dz 0.1) centred at 0.20 m
  const double e = r0 * p.fd * std::exp(-0.025 / p.ev_z) + r3 * p.fd * std::exp(-z3 / p.ev_z);
  CHECK(s.m.soilev == doctest::Approx(e).epsilon(1e-12));
  CHECK(s.m.sink0 == doctest::Approx(p.v_ground * p.gem * 30.0 * 86.4).epsilon(1e-12));
  CHECK(std::fabs(s.balance_error()) < 1e-13 * s.total());
  // frozen layers or full snow cover: no re-emission
  hgtr::State t; hgtr::initialize(p, t, a); hgtr::Flux g = f; g.ts[0] = -1.0; g.ts[3] = -0.1; hgtr::step(p, t, a, b, g); CHECK(t.m.soilev == 0.0);
  hgtr::State u; hgtr::initialize(p, u, a); hgtr::Flux k = f; k.swe = 100.0; hgtr::step(p, u, a, b, k); CHECK(u.m.soilev == 0.0);
  // off: the T1 fraction hg_f_evade
  hgtr::Params q = t2params(); hgtr::State v; hgtr::initialize(q, v, a); hgtr::step(q, v, a, b, f);
  CHECK(v.m.soilev == 0.0); CHECK(v.m.evasion == doctest::Approx((r0 + r3) * q.f_evade + v.m.snowred + v.m.photored).epsilon(1e-12));
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_evasion with a geometry change: decomposition-weighted re-emission of the pot; exact balance") {
  hgtr::Params p = t2params(); p.t2_evasion = 1;
  hgtr::State s; hgtr::Snap a = column(4); hgtr::initialize(p, s, a);
  hgtr::Snap b = a; b.nsl = 5; for (int l = 0; l < 5; l++) { b.soil[l] = (a.soil[0] + a.soil[1] + a.soil[2] + a.soil[3] - 40.0) / 5.0; b.dz[l] = 0.06; b.type[l] = (l < 2) ? 1 : 3; }
  hgtr::Flux f = t2flux(1, 1, 0, 0, 0, 0.3); f.decomp[0] = 20.0; f.ts[0] = 4.0;
  hgtr::step(p, s, a, b, f);
  CHECK(s.m.soilev > 0.0); CHECK(std::fabs(s.balance_error()) < 1e-13 * s.total());
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_foliar 2: stomatal Hg(0) uptake = k x GPP to the leaves (+ LAI interception of dry Hg(II)); linear in GPP") {
  hgtr::Params p = t2params(); p.t2_foliar = 2;
  hgtr::State s; hgtr::Snap a = column(4); hgtr::initialize(p, s, a);
  hgtr::Flux f = t2flux(0.0, 0.0, 0, 0, 0, 0.3); f.gpp[0] = 60.0; f.gpp[3] = 20.0; f.ta = 10.0; f.co2 = 400.0;
  const double k = hgtr::stomatal_k(p, f);
  const double kexp = 0.75 * 1.3 * 1.6 * 4.0 / (12.0 * 400e-6 * 3.3365 * (101325.0 / (8.314 * 283.15))) * 1e-3;
  CHECK(k == doctest::Approx(kexp).epsilon(1e-12));
  const double l0 = s.veg[0][0], l3 = s.veg[3][0];
  hgtr::step(p, s, a, a, f);
  CHECK(s.m.stomatal == doctest::Approx(80.0 * k).epsilon(1e-12));
  CHECK(s.veg[0][0] - l0 == doctest::Approx(60.0 * k).epsilon(1e-9)); CHECK(s.veg[3][0] - l3 == doctest::Approx(20.0 * k).epsilon(1e-9));
  CHECK(s.m.foliar == doctest::Approx(s.m.stomatal));   // LAI 0: no Hg(II) interception
  CHECK(std::fabs(s.balance_error()) < 1e-13 * s.total());
  // with LAI: + dry Hg(II) share; dry Hg(0) of the LAI scheme is not used
  hgtr::State t; hgtr::initialize(p, t, a); hgtr::Flux g = t2flux(2.0, 2.0, 0, 0, 0, 0.3); g.gpp[0] = 60.0; g.gpp[3] = 20.0; g.ta = 10.0; g.co2 = 400.0;
  hgtr::step(p, t, a, a, g);
  CHECK(t.m.foliar - t.m.stomatal == doctest::Approx(0.25 * (0.6 + 0.3) * 0.2).epsilon(1e-12));
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_meltroute: melt Hg delivered with the day's melt, overland share exported; pending pool empty after the month") {
  hgtr::Params p = t2dparams(); p.t2_meltroute = 1; p.d_cryo = 0; p.d_bio = 0;
  hgtr::State s; hgtr::Snap a = column(5); hgtr::initialize(p, s, a);
  s.snow = 2.0; s.h_init = s.total(); s.cum.clear();
  hgtr::Flux f = t2flux(0, 0, 20.0, 30.0, 0.0, 0.0);   // melt 30 of SWE 20 -> release 0.6 of the pack
  hgtr::step(p, s, a, a, f);
  CHECK(s.mpend > 0.0); CHECK(s.m.snowmelt == doctest::Approx(s.mpend)); CHECK(std::fabs(s.balance_error()) < 1e-13 * s.total());
  const double mp = s.mpend;
  hgtr::Flux t = dflux(5, 10, 2.0, 6.0, 0.0, 20.0); for (int d = 0; d < 10; d++) t.d_melt[d] = (d < 5) ? 3.0 : 0.0;
  p.tau_sorb = 1e5; hgtr::State u = s; u.m.clear(); hgtr::step(p, u, a, a, t);
  CHECK(u.mpend == 0.0);
  CHECK(u.m.meltrun == doctest::Approx(mp * 6.0 / 8.0).epsilon(1e-12));
  CHECK(u.m.runoff >= u.m.meltrun);
  CHECK(std::fabs(u.balance_error()) < 1e-13 * u.total());
  // pack / unpack at a month end keeps the total (mpend folded into the top dissolved pool)
  double v[hgtr::NSTATE]; hgtr::State w = s; hgtr::pack(w, v); hgtr::State x; REQUIRE(hgtr::unpack(x, v, 5));
  CHECK(x.total() == doctest::Approx(w.total()).epsilon(1e-14));
  std::remove(T2DEP);
}

TEST_CASE("hg_t2 phase 3 all on, random steps with daily transport: balance to roundoff, pools >= 0") {
  std::mt19937 rng(31); std::uniform_real_distribution<double> U(0.0, 1.0);
  hgtr::Params p = t2dparams(); p.t2_evasion = 1; p.t2_foliar = 2; p.t2_meltroute = 1;
  hgtr::State s; hgtr::Snap a = column(8, 3000.0); hgtr::initialize(p, s, a);
  for (int it = 0; it < 1500; it++) {
    const int mth = it % 12;
    hgtr::Snap b = a; hgtr::Flux f = t2flux(3 * U(rng), U(rng), (mth < 5 || mth > 9) ? 200 * U(rng) : 0.0, 30 * U(rng), 5 * U(rng), 0.4 * U(rng));
    for (int i = 0; i < NUM_PFT; i++) { f.gpp[i] = 50 * U(rng); }
    f.ta = 30 * U(rng) - 15; f.co2 = 300 + 200 * U(rng);
    for (int l = 0; l < a.nsl; l++) { f.ts[l] = 20 * U(rng) - 10; f.decomp[l] = 5 * U(rng); b.soil[l] = std::fmax(0.0, a.soil[l] - f.decomp[l] + 4 * U(rng)); }
    if (U(rng) < 0.3) {   // geometry change
      int n = 3 + int(U(rng) * 12); b.nsl = n; double c0 = 0; for (int l = 0; l < a.nsl; l++) c0 += a.soil[l];
      for (int l = 0; l < MAX_SOI_LAY; l++) { b.soil[l] = (l < n) ? c0 / n : 0; b.dz[l] = (l < n) ? 0.05 + 0.1 * U(rng) : 0; b.type[l] = (l == 0) ? 0 : ((l < 3) ? 1 : 3); }
    }
    hgtr::step(p, s, a, b, f);
    hgtr::Flux t = dflux(b.nsl, 30, 20 * U(rng), 10 * U(rng), 3 * U(rng), 5 + 40 * U(rng), (mth < 5 || mth > 9) ? int(U(rng) * 3) : 99);
    for (int d = 0; d < 30; d++) t.d_melt[d] = (U(rng) < 0.3) ? 5 * U(rng) : 0.0;
    hgtr::step(p, s, b, b, t);
    CHECK(s.mpend == 0.0);
    a = b;
    if (it % 12 == 11) { s.m.clear(); }
  }
  CHECK(std::fabs(s.balance_error()) < 1e-8 * std::fmax(1.0, s.total()));
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s.soil[l] >= 0.0); CHECK(s.diss[l] >= 0.0); }
  CHECK(s.cum.export_ > 0.0); CHECK(s.cum.evasion > 0.0);
  std::remove(T2DEP);
}

// ---- T2 phase 4: fire ash retention, RTS particulate export, canopy wash-off
namespace {
// whole top layer + PFT-0 leaves and stems removed (fire or slump); geometry change
void remove_top(const hgtr::Snap& a, hgtr::Snap& b) {
  b = a; b.nsl = a.nsl - 1;
  for (int l = 0; l < b.nsl; l++) { b.soil[l] = a.soil[l + 1]; b.dz[l] = a.dz[l + 1]; b.type[l] = a.type[l + 1]; }
  b.soil[b.nsl] = 0; b.dz[b.nsl] = 0; b.veg[0][0] = 0; b.veg[0][1] = 0;
}
hgtr::State scaled(const hgtr::State& s, double k) {
  hgtr::State t = s;
  for (int i = 0; i < NUM_PFT; i++) { t.dead[i] *= k; for (int j = 0; j < NUM_PFT_PART; j++) t.veg[i][j] *= k; }
  t.wd *= k; t.snow *= k; t.mpend *= k;
  for (int l = 0; l < MAX_SOI_LAY; l++) { t.soil[l] *= k; t.diss[l] *= k; t.mehg[l] *= k; }
  t.h_init = t.total(); t.cum.clear(); t.m.clear(); t.y.clear();
  return t;
}
}

TEST_CASE("hg_t2 phase 4 parse: defaults off, ranges, dependencies") {
  hgtr::Params p = t2params();
  CHECK(p.fire_ash == 0.0); CHECK(p.t2_rts == 0); CHECK(p.wash == 0.0); CHECK(p.wash_pref == 10.0);
  Json::Value a = ms0(); a["hg_tracer"] = 1; a["hg_fire_ash"] = 0.1; CHECK(hgtr::parse(a, p) != "");      // needs hg_t2
  Json::Value b = ms0(); b["hg_tracer"] = 1; b["hg_t2_rts"] = 1; CHECK(hgtr::parse(b, p) != "");
  Json::Value c = ms0(); c["hg_tracer"] = 1; c["hg_t2"] = 1; c["hg_dep_file"] = T2DEP; c["hg_t2_preset"] = 0; c["hg_fire_ash"] = 0.1; c["hg_t2_rts"] = 1; c["hg_wash"] = 0.5; c["hg_wash_pref"] = 20.0;
  CHECK(hgtr::parse(c, p) == ""); CHECK(p.fire_ash == 0.1); CHECK(p.t2_rts == 1); CHECK(p.wash == 0.5); CHECK(p.wash_pref == 20.0);
  c["hg_t2_foliar"] = 0; CHECK(hgtr::parse(c, p) != "");                                                   // wash-off needs a canopy scheme
  c["hg_t2_foliar"] = 1; c["hg_fire_ash"] = 1.5; CHECK(hgtr::parse(c, p) != "");
  std::remove(T2DEP);
}

TEST_CASE("hg_fire_ash: the burned soil share of the fire Hg stays in the ash on the new top layer; veg Hg all emitted") {
  hgtr::Params p0 = t2params(); hgtr::Params p = p0; p.fire_ash = 0.1;
  hgtr::Snap a = column(6), b; remove_top(a, b); b.dead[0] = a.dead[0] + 200.0;
  hgtr::Flux f; f.burn = a.soil[0] + 50.0;   // fire step: f.t2 = 0
  hgtr::State s0; hgtr::initialize(p0, s0, a); hgtr::State s; hgtr::initialize(p, s, a);
  const double htop = s.soil[0], pot = htop + 50 * p.r_leaf + 200 * p.r_stem;
  hgtr::step(p0, s0, a, b, f); hgtr::step(p, s, a, b, f);
  const double F = s0.m.fire;
  CHECK(F == doctest::Approx(pot * f.burn / (a.soil[0] + 250.0)));
  CHECK(s.m.ash == doctest::Approx(0.1 * F * htop / pot).epsilon(1e-12));
  CHECK(s.m.fire == doctest::Approx(F - s.m.ash).epsilon(1e-12));
  CHECK(s.soil[0] - s0.soil[0] == doctest::Approx(s.m.ash).epsilon(1e-9));   // to the new top C layer
  CHECK(std::fabs(s.balance_error()) < 1e-13 * s.total());
  // vegetation-only fire: nothing retained
  hgtr::Snap bv = a; bv.veg[0][0] = 0; hgtr::Flux fv; fv.burn = 50.0;
  hgtr::State v; hgtr::initialize(p, v, a); hgtr::step(p, v, a, bv, fv); CHECK(v.m.ash == 0.0); CHECK(v.m.fire > 0.0);
  // linear in Hg (affine spin-up): doubling the state doubles the fluxes
  hgtr::State d; hgtr::initialize(p, d, a); d = scaled(d, 2.0); hgtr::step(p, d, a, b, f);
  CHECK(d.m.ash == doctest::Approx(2.0 * s.m.ash).epsilon(1e-12)); CHECK(d.m.fire == doctest::Approx(2.0 * s.m.fire).epsilon(1e-12));
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_rts: the slumped material's Hg leaves as particulate export, respired-C Hg not returned to the column") {
  hgtr::Params p0 = t2params(); hgtr::Params p = p0; p.t2_rts = 1;
  hgtr::Snap a = column(6), b; remove_top(a, b);
  const double X = a.soil[0] + 250.0; hgtr::Flux f; f.lateral = 0.9 * X; f.pot_decomp = 0.1 * X;
  hgtr::State s0; hgtr::initialize(p0, s0, a); hgtr::State s; hgtr::initialize(p, s, a);
  const double h0 = s.total(), pot = s.soil[0] + 50 * p.r_leaf + 200 * p.r_stem;
  hgtr::step(p0, s0, a, b, f); hgtr::step(p, s, a, b, f);
  CHECK(s0.m.partic == 0.0);
  CHECK(s0.m.export_ == doctest::Approx(pot * (0.9 + 0.1 * p.f_export)).epsilon(1e-12));   // T1: respired-C Hg mostly back to the column
  CHECK(s.m.evasion == doctest::Approx(pot * 0.1 * p.f_evade).epsilon(1e-12));
  CHECK(s.m.partic == doctest::Approx(pot * (0.9 + 0.1 * (1.0 - p.f_evade))).epsilon(1e-12));
  CHECK(s.m.export_ == doctest::Approx(s.m.partic).epsilon(1e-14));
  CHECK(s.total() == doctest::Approx(h0 - pot).epsilon(1e-12));
  CHECK(std::fabs(s.balance_error()) < 1e-13 * s.total());
  hgtr::State d; hgtr::initialize(p, d, a); d = scaled(d, 3.0); hgtr::step(p, d, a, b, f);
  CHECK(d.m.partic == doctest::Approx(3.0 * s.m.partic).epsilon(1e-12));
  std::remove(T2DEP);
}

TEST_CASE("hg_wash: share wash x P / (P + P_ref) of the intercepted dry Hg(II) goes to the ground (both canopy schemes)") {
  for (int sch = 1; sch <= 2; sch++) {
    hgtr::Params p0 = t2params(); p0.t2_foliar = sch; hgtr::Params p = p0; p.wash = 0.5;
    hgtr::Snap a = column(4);
    hgtr::Flux f = t2flux(2.0, 2.0, 0, 0, 0, 0.3); f.gpp[0] = 60.0; f.ta = 10.0; f.co2 = 400.0; f.rain = 30.0;
    hgtr::State s0; hgtr::initialize(p0, s0, a); hgtr::State s; hgtr::initialize(p, s, a);
    hgtr::step(p0, s0, a, a, f); hgtr::step(p, s, a, a, f);
    const double w = 0.25 * (0.6 + 0.3) * 0.2 * 0.5 * 30.0 / 40.0;   // intercepted dry Hg(II) x fw
    CHECK(s0.m.wash == 0.0);
    CHECK(s.m.wash == doctest::Approx(w).epsilon(1e-12));
    CHECK(s0.m.foliar - s.m.foliar == doctest::Approx(w).epsilon(1e-12));
    CHECK(s.m.dep - s0.m.dep == doctest::Approx(w).epsilon(1e-12));
    CHECK(s.m.uptake + s.m.dep == doctest::Approx(s0.m.uptake + s0.m.dep).epsilon(1e-13));
    CHECK(std::fabs(s.balance_error()) < 1e-13 * s.total());
    hgtr::State t; hgtr::initialize(p, t, a); hgtr::Flux g = f; g.rain = 0.0; hgtr::step(p, t, a, a, g); CHECK(t.m.wash == 0.0);
  }
  std::remove(T2DEP);
}

TEST_CASE("hg_t2 phase 4 all on (with phase 3), random steps with fire / RTS and daily transport: balance to roundoff, pools >= 0") {
  std::mt19937 rng(47); std::uniform_real_distribution<double> U(0.0, 1.0);
  hgtr::Params p = t2dparams(); p.t2_evasion = 1; p.t2_foliar = 2; p.t2_meltroute = 1; p.fire_ash = 0.2; p.t2_rts = 1; p.wash = 0.5;
  hgtr::State s; hgtr::Snap a = column(8, 3000.0); hgtr::initialize(p, s, a);
  for (int it = 0; it < 1500; it++) {
    const int mth = it % 12;
    hgtr::Snap b = a; hgtr::Flux f = t2flux(3 * U(rng), U(rng), (mth < 5 || mth > 9) ? 200 * U(rng) : 0.0, 30 * U(rng), 5 * U(rng), 0.4 * U(rng));
    for (int i = 0; i < NUM_PFT; i++) { f.gpp[i] = 50 * U(rng); }
    f.ta = 30 * U(rng) - 15; f.co2 = 300 + 200 * U(rng); f.rain = 60 * U(rng);
    for (int l = 0; l < a.nsl; l++) { f.ts[l] = 20 * U(rng) - 10; f.decomp[l] = 5 * U(rng); b.soil[l] = std::fmax(0.0, a.soil[l] - f.decomp[l] + 4 * U(rng)); }
    hgtr::step(p, s, a, b, f);
    a = b;
    if (U(rng) < 0.15) {   // disturbance sub-step: fire or RTS with a geometry change
      hgtr::Snap c = a; hgtr::Flux g;
      int n = 3 + int(U(rng) * 12); c.nsl = n; double c0 = 0; for (int l = 0; l < a.nsl; l++) c0 += a.soil[l];
      const double c1 = c0 * (0.85 + 0.1 * U(rng));
      for (int l = 0; l < MAX_SOI_LAY; l++) { c.soil[l] = (l < n) ? c1 / n : 0; c.dz[l] = (l < n) ? 0.05 + 0.1 * U(rng) : 0; c.type[l] = (l == 0) ? 0 : ((l < 3) ? 1 : 3); }
      if (U(rng) < 0.5) { g.burn = (c0 - c1) * U(rng); } else { g.lateral = 0.9 * (c0 - c1) * U(rng); g.pot_decomp = 0.1 * (c0 - c1) * U(rng); }
      hgtr::step(p, s, a, c, g);
      a = c;
    }
    hgtr::Flux t = dflux(a.nsl, 30, 20 * U(rng), 10 * U(rng), 3 * U(rng), 5 + 40 * U(rng), (mth < 5 || mth > 9) ? int(U(rng) * 3) : 99);
    for (int d = 0; d < 30; d++) t.d_melt[d] = (U(rng) < 0.3) ? 5 * U(rng) : 0.0;
    hgtr::step(p, s, a, a, t);
    if (it % 12 == 11) { s.m.clear(); }
  }
  CHECK(std::fabs(s.balance_error()) < 1e-8 * std::fmax(1.0, s.total()));
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s.soil[l] >= 0.0); CHECK(s.diss[l] >= 0.0); }
  CHECK(s.cum.fire > 0.0); CHECK(s.cum.export_ > 0.0);
  std::remove(T2DEP);
}

// ---- T2 phase 5b: humic Kd offset, export path (flow-path retention), ground Hg(0) sink under snow, hg_t2 preset, T2 copies
TEST_CASE("hg_t2 phase 5 parse: preset with hg_t2 = 1 (explicit wins), hg_t2_preset 0, T2-only copies, T1 untouched, dependencies") {
  hgtr::Params p = t2params();   // writes T2DEP
  CHECK(p.t2_kdsplit == 0); CHECK(p.t2_xpath == 0); CHECK(p.t2_snowsink == 0); CHECK(p.f_conn == 1.0); CHECK(p.t2_dlog_kd_hum == 0.0);
  CHECK(p.r_stem == 0.02); CHECK(p.r_root == 0.07); CHECK(p.k_meth_proxy == 87.0); CHECK(p.t2_diss == 0);   // preset 0: component defaults
  Json::Value t1 = ms0(); t1["hg_tracer"] = 1; t1["hg_mehg"] = 1; t1["hg_r_stem"] = 0.01; t1["hg_k_meth_proxy"] = 50.0;
  hgtr::Params q; REQUIRE(hgtr::parse(t1, q) == "");   // T1 path: T1 names, no preset
  CHECK(q.r_stem == 0.01); CHECK(q.k_meth_proxy == 50.0); CHECK(q.t2_diss == 0); CHECK(q.t2_evasion == 0); CHECK(q.v_ground == 8.5e-5); CHECK(q.t2_log_kd_org == 3.5);
  Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_t2"] = 1; m["hg_dep_file"] = T2DEP;
  hgtr::Params r; REQUIRE(hgtr::parse(m, r) == "");   // preset
  CHECK(r.t2_diss == 1); CHECK(r.t2_evasion == 1); CHECK(r.t2_foliar == 2); CHECK(r.t2_meltroute == 1); CHECK(r.t2_rts == 1);
  CHECK(r.fire_ash == 0.1); CHECK(r.wash == 0.5);
  hgtr::Params pp; hgtr::t2_preset_(pp);
  CHECK(r.t2_log_kd_org == pp.t2_log_kd_org); CHECK(r.v_ground == pp.v_ground); CHECK(r.fd == pp.fd); CHECK(r.ev_z == pp.ev_z);
  CHECK(r.t2_kdsplit == pp.t2_kdsplit); CHECK(r.t2_xpath == pp.t2_xpath); CHECK(r.t2_snowsink == pp.t2_snowsink);
  CHECK(r.r_stem == pp.t2_r_stem); CHECK(r.r_root == pp.t2_r_root); CHECK(r.k_meth_proxy == pp.t2_k_meth_proxy);
  Json::Value e = m; e["hg_v_ground"] = 2e-4; e["hg_t2_rts"] = 0; e["hg_t2_r_stem"] = 0.01;
  REQUIRE(hgtr::parse(e, r) == ""); CHECK(r.v_ground == 2e-4); CHECK(r.t2_rts == 0); CHECK(r.r_stem == 0.01); CHECK(r.t2_diss == 1);   // explicit wins
  Json::Value g = m; g["hg_t2_diss"] = 0; REQUIRE(hgtr::parse(g, r) == ""); CHECK(r.t2_meltroute == 0); CHECK(r.t2_xpath == 0); CHECK(r.t2_kdsplit == 0);
  Json::Value g2 = m; g2["hg_t2_evasion"] = 0; REQUIRE(hgtr::parse(g2, r) == ""); CHECK(r.t2_snowsink == 0);
  Json::Value h = m; h["hg_r_stem"] = 0.01; CHECK(hgtr::parse(h, r) != "");          // T1 name with hg_t2 = 1
  Json::Value h2 = m; h2["hg_k_meth_proxy"] = 70.0; CHECK(hgtr::parse(h2, r) != "");
  Json::Value k = m; k["hg_t2_preset"] = 0; k["hg_t2_kdsplit"] = 1; CHECK(hgtr::parse(k, r) != "");   // needs hg_t2_diss
  k["hg_t2_diss"] = 1; k["hg_t2_dlog_kd_hum"] = 0.5; REQUIRE(hgtr::parse(k, r) == ""); CHECK(r.t2_kdsplit == 1); CHECK(r.t2_dlog_kd_hum == 0.5);
  Json::Value x = m; x["hg_t2_preset"] = 0; x["hg_t2_xpath"] = 1; CHECK(hgtr::parse(x, r) != "");
  Json::Value y = m; y["hg_t2_preset"] = 0; y["hg_t2_snowsink"] = 1; CHECK(hgtr::parse(y, r) != "");   // needs hg_t2_evasion
  y["hg_t2_evasion"] = 1; y["hg_v_snow"] = 3e-4; REQUIRE(hgtr::parse(y, r) == ""); CHECK(r.v_snow == 3e-4);
  Json::Value z = m; z["hg_f_conn"] = 1.5; CHECK(hgtr::parse(z, r) != "");
  Json::Value w = ms0(); w["hg_tracer"] = 1; w["hg_t2_xpath"] = 1; CHECK(hgtr::parse(w, r) != "");    // needs hg_t2
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_snowsink: ground sink ((1 - f_snow) v_ground + f_snow v_snow) x GEM x s; v_snow = v_ground == off") {
  hgtr::Params p = t2params(); p.t2_evasion = 1; p.t2_foliar = 0;
  hgtr::Snap a = column(4); hgtr::Flux f = t2flux(0, 0, 5.0, 0, 0, 0.3);   // SWE 5 of swe_full 10: f_snow 0.5
  hgtr::Params q = p; q.t2_snowsink = 1; q.v_snow = 3.0e-4;
  hgtr::State s; hgtr::initialize(q, s, a); hgtr::step(q, s, a, a, f);
  CHECK(s.m.sink0 == doctest::Approx((0.5 * p.v_ground + 0.5 * 3.0e-4) * p.gem * 30.0 * 86.4).epsilon(1e-12));
  CHECK(std::fabs(s.balance_error()) < 1e-13 * s.total());
  hgtr::Params r = q; r.v_snow = r.v_ground;
  hgtr::State u, v; hgtr::initialize(p, u, a); hgtr::initialize(r, v, a); hgtr::step(p, u, a, a, f); hgtr::step(r, v, a, a, f);
  CHECK(v.m.sink0 == doctest::Approx(u.m.sink0).epsilon(1e-15)); CHECK(v.total() == doctest::Approx(u.total()).epsilon(1e-15));
  hgtr::Flux g = t2flux(0, 0, 0.0, 0, 0, 0.3); hgtr::State w; hgtr::initialize(q, w, a); hgtr::step(q, w, a, a, g);   // no snow: v_ground only
  CHECK(w.m.sink0 == doctest::Approx(p.v_ground * p.gem * 30.0 * 86.4).epsilon(1e-12));
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_kdsplit: humic layers sorb with Kd x 10^dlog; dlog 0 == off (bit-identical)") {
  hgtr::Params p = t2dparams(); p.d_cryo = 0; p.d_bio = 0;
  hgtr::Snap a = column(5); hgtr::Flux t = dflux(5, 10, 2.0, 0.0, 1.0, 20.0);
  auto run = [&](const hgtr::Params& q) { hgtr::State s; hgtr::initialize(q, s, a); for (int l = 0; l < 5; l++) { s.diss[l] = 0.5 * s.soil[l]; s.soil[l] *= 0.5; }
    s.h_init = s.total(); hgtr::step(q, s, a, a, t); return s; };
  hgtr::State s0 = run(p); hgtr::Params q = p; q.t2_kdsplit = 1; hgtr::State s1 = run(q);
  for (int l = 0; l < 5; l++) { CHECK(s1.diss[l] == s0.diss[l]); CHECK(s1.soil[l] == s0.soil[l]); }
  CHECK(s1.m.export_ == s0.m.export_);
  q.t2_dlog_kd_hum = 1.0; hgtr::State s2 = run(q);
  CHECK(s2.diss[1] < s0.diss[1]); CHECK(s2.diss[2] < s0.diss[2]);    // humic (type 2) layers: more sorbed
  CHECK(s2.m.drain < s0.m.drain);                                    // the drained layer 2 is humic
  CHECK(std::fabs(s2.balance_error()) < 1e-13 * s2.total());
  q.t2_dlog_kd_hum = -1.0; hgtr::State s3 = run(q); CHECK(s3.diss[2] > s0.diss[2]); CHECK(std::fabs(s3.balance_error()) < 1e-13 * s3.total());
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_xpath: share f_conn of the organic-layer overland / lateral export reaches the stream, the rest goes to the first mineral layer") {
  hgtr::Params p = t2dparams(); p.d_cryo = 0; p.d_bio = 0; p.tau_sorb = 1e5;   // no sorption during the month: transport only
  hgtr::Snap a = column(5); hgtr::Flux t = dflux(5, 10, 2.0, 6.0, 1.0, 20.0);
  auto run = [&](const hgtr::Params& q) { hgtr::State s; hgtr::initialize(q, s, a); for (int l = 0; l < 5; l++) { s.diss[l] = 0.5 * s.soil[l]; s.soil[l] *= 0.5; }
    s.h_init = s.total(); hgtr::step(q, s, a, a, t); return s; };
  hgtr::State s0 = run(p); hgtr::Params q = p; q.t2_xpath = 1; q.f_conn = 1.0; hgtr::State s1 = run(q);
  for (int l = 0; l < 5; l++) { CHECK(s1.diss[l] == s0.diss[l]); } CHECK(s1.m.export_ == s0.m.export_); CHECK(s1.m.down == s0.m.down);
  q.f_conn = 0.4; hgtr::State s2 = run(q);
  CHECK(s2.m.export_ < s0.m.export_); CHECK(s2.m.down > s0.m.down); CHECK(s2.diss[3] > s0.diss[3]);
  CHECK(std::fabs(s2.balance_error()) < 1e-13 * s2.total());
  q.f_conn = 0.0; hgtr::State s3 = run(q); CHECK(s3.m.runoff == 0.0); CHECK(s3.m.drain == 0.0);   // organic layers only: column(5) mineral layers 3-4 not drained
  CHECK(std::fabs(s3.balance_error()) < 1e-13 * s3.total());
  // melt overland share retained too
  hgtr::Params r = p; r.t2_meltroute = 1; r.t2_xpath = 1; r.f_conn = 0.25;
  hgtr::State u; hgtr::initialize(r, u, a); u.mpend = 1.0; u.h_init = u.total(); hgtr::Flux m = t; for (int d = 0; d < 10; d++) m.d_melt[d] = 3.0;
  hgtr::step(r, u, a, a, m); CHECK(u.m.meltrun == doctest::Approx(0.25 * 6.0 / 8.0).epsilon(1e-12)); CHECK(std::fabs(u.balance_error()) < 1e-13 * u.total());
  std::remove(T2DEP);
}

TEST_CASE("hg_t2 phase 5b all on (with phases 3-4), random steps: balance to roundoff, pools >= 0") {
  std::mt19937 rng(59); std::uniform_real_distribution<double> U(0.0, 1.0);
  hgtr::Params p = t2dparams(); p.t2_evasion = 1; p.t2_foliar = 2; p.t2_meltroute = 1; p.fire_ash = 0.2; p.t2_rts = 1; p.wash = 0.5;
  p.t2_kdsplit = 1; p.t2_dlog_kd_hum = 0.7; p.t2_xpath = 1; p.f_conn = 0.3; p.t2_snowsink = 1; p.v_snow = 4e-4;
  hgtr::State s; hgtr::Snap a = column(8, 3000.0); hgtr::initialize(p, s, a);
  for (int it = 0; it < 1500; it++) {
    const int mth = it % 12;
    hgtr::Snap b = a; hgtr::Flux f = t2flux(3 * U(rng), U(rng), (mth < 5 || mth > 9) ? 200 * U(rng) : 0.0, 30 * U(rng), 5 * U(rng), 0.4 * U(rng));
    for (int i = 0; i < NUM_PFT; i++) { f.gpp[i] = 50 * U(rng); }
    f.ta = 30 * U(rng) - 15; f.co2 = 300 + 200 * U(rng); f.rain = 60 * U(rng);
    for (int l = 0; l < a.nsl; l++) { f.ts[l] = 20 * U(rng) - 10; f.decomp[l] = 5 * U(rng); b.soil[l] = std::fmax(0.0, a.soil[l] - f.decomp[l] + 4 * U(rng)); }
    if (U(rng) < 0.3) {
      int n = 3 + int(U(rng) * 12); b.nsl = n; double c0 = 0; for (int l = 0; l < a.nsl; l++) c0 += a.soil[l];
      for (int l = 0; l < MAX_SOI_LAY; l++) { b.soil[l] = (l < n) ? c0 / n : 0; b.dz[l] = (l < n) ? 0.05 + 0.1 * U(rng) : 0; b.type[l] = (l == 0) ? 0 : ((l < 2) ? 1 : ((l < 3) ? 2 : 3)); }
    }
    hgtr::step(p, s, a, b, f);
    hgtr::Flux t = dflux(b.nsl, 30, 20 * U(rng), 10 * U(rng), 3 * U(rng), 5 + 40 * U(rng), (mth < 5 || mth > 9) ? int(U(rng) * 3) : 99);
    for (int d = 0; d < 30; d++) t.d_melt[d] = (U(rng) < 0.3) ? 5 * U(rng) : 0.0;
    hgtr::step(p, s, b, b, t);
    a = b;
    if (it % 12 == 11) { s.m.clear(); }
  }
  CHECK(std::fabs(s.balance_error()) < 1e-8 * std::fmax(1.0, s.total()));
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s.soil[l] >= 0.0); CHECK(s.diss[l] >= 0.0); }
  CHECK(s.cum.export_ > 0.0); CHECK(s.cum.down > 0.0);
  std::remove(T2DEP);
}

// ---- T2 phase 5c: flow-path retention routed below the mixing zone (hg_t2_deeproute, hg_f_deep, hg_z_deep)
TEST_CASE("hg_t2_deeproute parse: default off, needs hg_t2_xpath, follows an explicit xpath off under the preset, ranges") {
  hgtr::Params p = t2params();   // writes T2DEP
  CHECK(p.t2_deeproute == 0); CHECK(p.f_deep == 1.0); CHECK(p.z_deep == 1.0);
  Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_t2"] = 1; m["hg_dep_file"] = T2DEP;
  hgtr::Params r; hgtr::Params pp; hgtr::t2_preset_(pp);
  REQUIRE(hgtr::parse(m, r) == ""); CHECK(r.t2_deeproute == pp.t2_deeproute); CHECK(r.z_deep == pp.z_deep); CHECK(r.f_deep == pp.f_deep);
  Json::Value a = m; a["hg_t2_preset"] = 0; a["hg_t2_diss"] = 1; a["hg_t2_deeproute"] = 1; CHECK(hgtr::parse(a, r) != "");   // needs xpath
  a["hg_t2_xpath"] = 1; a["hg_z_deep"] = 1.5; a["hg_f_deep"] = 0.5; REQUIRE(hgtr::parse(a, r) == "");
  CHECK(r.t2_deeproute == 1); CHECK(r.z_deep == 1.5); CHECK(r.f_deep == 0.5);
  Json::Value b = m; b["hg_t2_xpath"] = 0; REQUIRE(hgtr::parse(b, r) == ""); CHECK(r.t2_deeproute == 0);   // follows xpath off
  Json::Value c = m; c["hg_f_deep"] = 1.5; CHECK(hgtr::parse(c, r) != "");
  Json::Value d = m; d["hg_z_deep"] = -0.1; CHECK(hgtr::parse(d, r) != "");
  Json::Value w = ms0(); w["hg_tracer"] = 1; w["hg_t2_deeproute"] = 1; CHECK(hgtr::parse(w, r) != "");    // needs hg_t2
  Json::Value t1 = ms0(); t1["hg_tracer"] = 1; t1["hg_mehg"] = 1; hgtr::Params q; REQUIRE(hgtr::parse(t1, q) == ""); CHECK(q.t2_deeproute == 0);
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_deeproute: retained Hg goes to the layer at hg_z_deep (f_deep share); off / f_deep 0 / z_deep above the mineral == phase 5b") {
  hgtr::Params p = t2dparams(); p.d_cryo = 0; p.d_bio = 0; p.tau_sorb = 1e5; p.t2_xpath = 1; p.f_conn = 0.4;
  hgtr::Snap a = column(6); hgtr::Flux t = dflux(6, 10, 2.0, 6.0, 1.0, 20.0);   // layer tops 0, .05, .10 organic; .15, .25, .35 mineral
  auto run = [&](const hgtr::Params& q) { hgtr::State s; hgtr::initialize(q, s, a); for (int l = 0; l < 6; l++) { s.diss[l] = 0.5 * s.soil[l]; s.soil[l] *= 0.5; }
    s.h_init = s.total(); hgtr::step(q, s, a, a, t); return s; };
  hgtr::State s0 = run(p);
  for (double fd : {0.0, 1.0}) for (double zd : {0.0, 0.3}) {
    if (fd > 0.0 && zd > 0.0) continue;
    hgtr::Params q = p; q.t2_deeproute = 1; q.f_deep = fd; q.z_deep = zd; hgtr::State s1 = run(q);   // f_deep 0 or receiving layer = first mineral
    for (int l = 0; l < 6; l++) { CHECK(s1.diss[l] == s0.diss[l]); CHECK(s1.soil[l] == s0.soil[l]); }
    CHECK(s1.m.export_ == s0.m.export_); CHECK(s1.m.down == s0.m.down);
  }
  hgtr::Params q = p; q.t2_deeproute = 1; q.f_deep = 1.0; q.z_deep = 0.3; hgtr::State s2 = run(q);   // layer 4 (top 0.25)
  for (int l = 0; l < 3; l++) { CHECK(s2.diss[l] == s0.diss[l]); }   // organic layers unchanged
  CHECK(s2.m.export_ == s0.m.export_); CHECK(s2.m.down == s0.m.down);
  CHECK(s2.diss[3] < s0.diss[3]); CHECK(s2.diss[4] > s0.diss[4]);
  CHECK(std::fabs(s2.balance_error()) < 1e-13 * s2.total());
  q.z_deep = 50.0; hgtr::State s3 = run(q); CHECK(s3.diss[5] > s0.diss[5]); CHECK(std::fabs(s3.balance_error()) < 1e-13 * s3.total());   // deepest layer
  q.z_deep = 0.3; q.f_deep = 0.5; hgtr::State s4 = run(q); CHECK(s4.diss[4] > s0.diss[4]); CHECK(s4.diss[4] < s2.diss[4]);
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_deeproute below the thawed zone: routed Hg stays in the frozen layer, not mixed back into the organic layers") {
  hgtr::Params p = t2dparams(); p.tau_sorb = 1e5; p.t2_xpath = 1; p.f_conn = 0.4; p.d_cryo = 1e-2; p.d_bio = 1e-3;   // strong mixing
  hgtr::Snap a = column(6); hgtr::Flux t = dflux(6, 10, 2.0, 6.0, 1.0, 20.0, 5);   // layer 5 (top 0.35) frozen
  auto run = [&](const hgtr::Params& q) { hgtr::State s; hgtr::initialize(q, s, a); for (int l = 0; l < 6; l++) { s.diss[l] = 0.5 * s.soil[l]; s.soil[l] *= 0.5; }
    for (int l = 3; l < 6; l++) { s.soil[l] *= 0.01; s.diss[l] *= 0.01; }   // Hg-poor mineral soil: mixing goes downward
    s.h_init = s.total(); hgtr::step(q, s, a, a, t); return s; };
  hgtr::Params off = p; off.t2_xpath = 0; hgtr::State s0 = run(off);
  hgtr::State s1 = run(p);                                                   // phase 5b: to the first mineral layer
  hgtr::Params q = p; q.t2_deeproute = 1; q.z_deep = 0.4; hgtr::State s2 = run(q);   // to the frozen layer 5
  auto org = [](const hgtr::State& s) { double x = 0; for (int l = 0; l < 3; l++) x += s.soil[l] + s.diss[l]; return x; };
  auto lay = [](const hgtr::State& s, int l) { return s.soil[l] + s.diss[l]; };
  CHECK(org(s2) < org(s1));                                   // the phase-5b retention partly mixes back up
  CHECK(org(s2) == doctest::Approx(org(s0)).epsilon(1e-12));  // deep route: organic layers as with no retention
  CHECK(lay(s2, 5) - lay(s0, 5) == doctest::Approx(s2.m.down - s0.m.down).epsilon(1e-10));   // all retention in layer 5
  CHECK(s2.m.export_ < s0.m.export_);
  CHECK(std::fabs(s2.balance_error()) < 1e-13 * s2.total());
  std::remove(T2DEP);
}

TEST_CASE("hg_t2 phase 5c all on (deep route), random steps with geometry changes: balance to roundoff, pools >= 0") {
  std::mt19937 rng(61); std::uniform_real_distribution<double> U(0.0, 1.0);
  hgtr::Params p = t2dparams(); p.t2_evasion = 1; p.t2_foliar = 2; p.t2_meltroute = 1; p.fire_ash = 0.2; p.t2_rts = 1; p.wash = 0.5;
  p.t2_kdsplit = 1; p.t2_dlog_kd_hum = 0.4; p.t2_xpath = 1; p.f_conn = 0.3; p.t2_snowsink = 1; p.v_snow = 4e-4;
  p.t2_deeproute = 1; p.f_deep = 0.8; p.z_deep = 0.5;
  hgtr::State s; hgtr::Snap a = column(8, 3000.0); hgtr::initialize(p, s, a);
  for (int it = 0; it < 1500; it++) {
    const int mth = it % 12;
    hgtr::Snap b = a; hgtr::Flux f = t2flux(3 * U(rng), U(rng), (mth < 5 || mth > 9) ? 200 * U(rng) : 0.0, 30 * U(rng), 5 * U(rng), 0.4 * U(rng));
    for (int i = 0; i < NUM_PFT; i++) { f.gpp[i] = 50 * U(rng); }
    f.ta = 30 * U(rng) - 15; f.co2 = 300 + 200 * U(rng); f.rain = 60 * U(rng);
    for (int l = 0; l < a.nsl; l++) { f.ts[l] = 20 * U(rng) - 10; f.decomp[l] = 5 * U(rng); b.soil[l] = std::fmax(0.0, a.soil[l] - f.decomp[l] + 4 * U(rng)); }
    if (U(rng) < 0.3) {
      int n = 3 + int(U(rng) * 12); b.nsl = n; double c0 = 0; for (int l = 0; l < a.nsl; l++) c0 += a.soil[l];
      for (int l = 0; l < MAX_SOI_LAY; l++) { b.soil[l] = (l < n) ? c0 / n : 0; b.dz[l] = (l < n) ? 0.05 + 0.1 * U(rng) : 0; b.type[l] = (l == 0) ? 0 : ((l < 2) ? 1 : ((l < 3) ? 2 : 3)); }
    }
    hgtr::step(p, s, a, b, f);
    hgtr::Flux t = dflux(b.nsl, 30, 20 * U(rng), 10 * U(rng), 3 * U(rng), 5 + 40 * U(rng), (mth < 5 || mth > 9) ? int(U(rng) * 3) : 3 + int(U(rng) * 6));
    for (int d = 0; d < 30; d++) t.d_melt[d] = (U(rng) < 0.3) ? 5 * U(rng) : 0.0;
    hgtr::step(p, s, b, b, t);
    a = b;
    if (it % 12 == 11) { s.m.clear(); }
  }
  CHECK(std::fabs(s.balance_error()) < 1e-8 * std::fmax(1.0, s.total()));
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s.soil[l] >= 0.0); CHECK(s.diss[l] >= 0.0); }
  CHECK(s.cum.export_ > 0.0); CHECK(s.cum.down > 0.0);
  std::remove(T2DEP);
}

// ---- T2 phase 5d: T2-only copy of hg_k_meth_ch4
TEST_CASE("hg_t2_k_meth_ch4: T2 copy (default 3250 = T1 value), T1 name with hg_t2 = 1 is an error, T1 path untouched") {
  hgtr::Params p = t2params();   // writes T2DEP
  CHECK(p.t2_k_meth_ch4 == 3250.0); CHECK(p.k_meth_ch4 == 3250.0);
  Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_t2"] = 1; m["hg_dep_file"] = T2DEP;
  hgtr::Params r; REQUIRE(hgtr::parse(m, r) == ""); CHECK(r.k_meth_ch4 == 3250.0);
  Json::Value a = m; a["hg_t2_k_meth_ch4"] = 2000.0; REQUIRE(hgtr::parse(a, r) == ""); CHECK(r.k_meth_ch4 == 2000.0); CHECK(r.t2_k_meth_ch4 == 2000.0);
  Json::Value b = m; b["hg_k_meth_ch4"] = 3250.0; CHECK(hgtr::parse(b, r) != "");   // T1 name with hg_t2 = 1
  Json::Value c = m; c["hg_t2_k_meth_ch4"] = -1.0; CHECK(hgtr::parse(c, r) != "");
  Json::Value t1 = ms0(); t1["hg_tracer"] = 1; t1["hg_mehg"] = 1; t1["hg_k_meth_ch4"] = 1000.0; t1["hg_t2_k_meth_ch4"] = 5000.0;
  hgtr::Params q; REQUIRE(hgtr::parse(t1, q) == ""); CHECK(q.k_meth_ch4 == 1000.0);   // hg_t2 = 0: T1 name, the copy unused
  // same MeHg as T1 with the default copy: CH4-driven step, hg_t2 = 1 parse vs T1 value
  hgtr::Params s1; hgtr::Params s0; s1.tracer = s0.tracer = 1; s1.mehg = s0.mehg = 1; s1.k_meth_ch4 = r.k_meth_ch4; s0.k_meth_ch4 = 3250.0;
  CHECK(s1.k_meth_ch4 == s0.k_meth_ch4);
  std::remove(T2DEP);
}

// ---- T2 phase 6b: deep spread of the routed retention (hg_t2_deepspread) and thaw release (hg_t2_thawrel)
TEST_CASE("hg_t2 phase 6b/6c parse: component defaults off, preset on (6c), Schaefer fractions, dependencies, ranges") {
  hgtr::Params p = t2params();   // writes T2DEP
  CHECK(p.t2_deepspread == 0); CHECK(p.t2_thawrel == 0); CHECK(p.f_thaw_evade == 0.16); CHECK(p.f_thaw_export == 0.015); CHECK(p.thaw_ref_year == 2024);
  Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_t2"] = 1; m["hg_dep_file"] = T2DEP;
  hgtr::Params r; REQUIRE(hgtr::parse(m, r) == ""); CHECK(r.t2_deepspread == 1); CHECK(r.t2_thawrel == 1);   // preset (6c): on
  { hgtr::Params pp; hgtr::t2_preset_(pp); CHECK(pp.t2_deepspread == 1); CHECK(pp.t2_thawrel == 1); }
  Json::Value z = m; z["hg_t2_preset"] = 0; REQUIRE(hgtr::parse(z, r) == ""); CHECK(r.t2_deepspread == 0); CHECK(r.t2_thawrel == 0);   // preset 0: off
  Json::Value z2 = m; z2["hg_t2_deepspread"] = 0; z2["hg_t2_thawrel"] = 0; REQUIRE(hgtr::parse(z2, r) == ""); CHECK(r.t2_deepspread == 0); CHECK(r.t2_thawrel == 0);   // explicit wins
  Json::Value z3 = m; z3["hg_t2_diss"] = 0; REQUIRE(hgtr::parse(z3, r) == ""); CHECK(r.t2_thawrel == 0);   // preset: follows diss off
  Json::Value a = m; a["hg_t2_deepspread"] = 1; a["hg_t2_thawrel"] = 1; a["hg_f_thaw_evade"] = 0.2; a["hg_f_thaw_export"] = 0.05; a["hg_t2_thaw_ref_year"] = 2010;
  REQUIRE(hgtr::parse(a, r) == ""); CHECK(r.t2_deepspread == 1); CHECK(r.t2_thawrel == 1); CHECK(r.f_thaw_evade == 0.2); CHECK(r.f_thaw_export == 0.05); CHECK(r.thaw_ref_year == 2010);
  Json::Value b = m; b["hg_t2_deeproute"] = 0; b["hg_t2_deepspread"] = 1; CHECK(hgtr::parse(b, r) != "");   // needs deeproute
  Json::Value b2 = m; b2["hg_t2_xpath"] = 0; REQUIRE(hgtr::parse(b2, r) == ""); CHECK(r.t2_deepspread == 0);   // preset: follows
  Json::Value c = m; c["hg_t2_diss"] = 0; c["hg_t2_thawrel"] = 1; CHECK(hgtr::parse(c, r) != "");   // needs diss
  Json::Value c2 = m; c2["hg_t2_thawrel"] = 1; c2["hg_transport"] = 0; CHECK(hgtr::parse(c2, r) != "");   // needs transport
  Json::Value d = m; d["hg_f_thaw_evade"] = 0.9; d["hg_f_thaw_export"] = 0.2; CHECK(hgtr::parse(d, r) != "");   // sum > 1
  Json::Value e = m; e["hg_f_thaw_evade"] = -0.1; CHECK(hgtr::parse(e, r) != "");
  Json::Value w = ms0(); w["hg_tracer"] = 1; w["hg_t2_thawrel"] = 1; CHECK(hgtr::parse(w, r) != "");   // needs hg_t2
  std::remove(T2DEP);
}

TEST_CASE("hg_t2 phase 6b thaw_track: daily max thaw depth -> zmax (always), zpf (years <= ref); off: no change") {
  hgtr::Params p = t2dparams(); hgtr::Snap a = column(6);   // tops 0 .05 .10 | .15 .25 .35, bottom .45
  hgtr::State s; hgtr::initialize(p, s, a);
  hgtr::Flux t = dflux(6, 10, 2.0, 0.0, 0.0, 20.0, 4); t.year = 2000;
  hgtr::thaw_track(p, s, a, t); CHECK(s.zmax == -1.0); CHECK(s.zpf == -1.0);   // switches off
  p.t2_thawrel = 1; hgtr::thaw_track(p, s, a, t); CHECK(s.zmax == doctest::Approx(0.25)); CHECK(s.zpf == doctest::Approx(0.25));
  hgtr::Flux u = dflux(6, 10, 2.0, 0.0, 0.0, 20.0, 5); u.year = 2030;   // deeper thaw after the reference period
  hgtr::thaw_track(p, s, a, u); CHECK(s.zmax == doctest::Approx(0.35)); CHECK(s.zpf == doctest::Approx(0.25));
  hgtr::Flux v = dflux(6, 10, 2.0, 0.0, 0.0, 20.0, 1); v.year = 2000; for (int l = 0; l < 6; l++) { v.d_ice[3][l] = 0.0; }   // one thawed day
  hgtr::State s2; hgtr::initialize(p, s2, a); hgtr::thaw_track(p, s2, a, v); CHECK(s2.zmax == doctest::Approx(0.45));   // daily maximum
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_deepspread: routed share to the frozen layers below zmax by soil C; zmax unset / no layer == deeproute") {
  hgtr::Params p = t2dparams(); p.d_cryo = 0; p.d_bio = 0; p.tau_sorb = 1e5; p.t2_xpath = 1; p.f_conn = 0.4;
  p.t2_deeproute = 1; p.f_deep = 1.0; p.z_deep = 0.3;
  hgtr::Snap a = column(6); hgtr::Flux t = dflux(6, 10, 2.0, 6.0, 1.0, 20.0, 4);   // layers 4-5 frozen
  auto run = [&](const hgtr::Params& q, double zmax) { hgtr::State s; hgtr::initialize(q, s, a); for (int l = 0; l < 6; l++) { s.diss[l] = 0.5 * s.soil[l]; s.soil[l] *= 0.5; }
    s.zmax = zmax; s.h_init = s.total(); hgtr::step(q, s, a, a, t); return s; };
  hgtr::State s0 = run(p, 0.25);                         // deeproute: all to layer 4
  hgtr::Params q = p; q.t2_deepspread = 1;
  hgtr::State s1 = run(q, -1.0), s2 = run(q, 10.0);     // unset / no layer below zmax: as deeproute (bit-identical)
  for (int l = 0; l < 6; l++) { CHECK(s1.diss[l] == s0.diss[l]); CHECK(s1.soil[l] == s0.soil[l]); CHECK(s2.diss[l] == s0.diss[l]); }
  hgtr::State s3 = run(q, 0.25);                         // layers 4 (C 1400) and 5 (C 1500)
  const double d4 = s3.diss[4] + s3.soil[4] - s0.diss[4] - s0.soil[4], d5 = s3.diss[5] + s3.soil[5] - s0.diss[5] - s0.soil[5];
  CHECK(d5 > 0.0); CHECK(d4 < 0.0);
  const double moved = -d4;   // routed Hg that went to layer 5 instead of 4: share 1500 / 2900
  hgtr::Params q0 = p; q0.t2_xpath = 1; q0.t2_deeproute = 0; hgtr::State sx = run(q0, 0.25);   // retention to lm (layer 3)
  const double routed = (s0.diss[4] + s0.soil[4]) - (sx.diss[4] + sx.soil[4]);
  CHECK(routed > 0.0); CHECK(d5 == doctest::Approx(routed * 1500.0 / 2900.0).epsilon(1e-6)); CHECK(moved == doctest::Approx(d5).epsilon(1e-9));
  for (int l = 0; l < 4; l++) { CHECK(s3.diss[l] == s0.diss[l]); }
  CHECK(s3.m.export_ == s0.m.export_); CHECK(s3.m.down == s0.m.down);
  CHECK(std::fabs(s3.balance_error()) < 1e-13 * s3.total());
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_thawrel: decay-released Hg below zpf evades / is exported at the Schaefer fractions; zpf unset or below the column == off") {
  hgtr::Params p = t2dparams(); p.t2_evasion = 1;
  hgtr::Snap a = column(6);
  hgtr::Flux f = t2flux(1.0, 0.0, 0.0, 0.0, 0.0, 0.2);
  for (int l = 0; l < 6; l++) { f.decomp[l] = 5.0; f.ts[l] = 2.0; }
  hgtr::Snap b = a; for (int l = 0; l < 6; l++) { b.soil[l] = a.soil[l] - 5.0; }
  auto run = [&](const hgtr::Params& q, double zpf) { hgtr::State s; hgtr::initialize(q, s, a); s.zpf = zpf; s.h_init = s.total(); hgtr::step(q, s, a, b, f); return s; };
  hgtr::State s0 = run(p, 0.25);
  hgtr::Params q = p; q.t2_thawrel = 1;
  hgtr::State s1 = run(q, -1.0), s2 = run(q, 5.0);
  for (int l = 0; l < 6; l++) { CHECK(s1.soil[l] == s0.soil[l]); CHECK(s1.diss[l] == s0.diss[l]); CHECK(s2.soil[l] == s0.soil[l]); }
  CHECK(s1.m.evasion == s0.m.evasion); CHECK(s1.m.export_ == s0.m.export_); CHECK(s1.m.thawev == 0.0);
  hgtr::State s3 = run(q, 0.25);   // layers 4, 5 (tops .25, .35)
  const double r4 = 5.0 * p.r_min, r5 = 5.0 * p.r_min;   // rel = decomp x Hg:C (mineral initial ratio)
  CHECK(s3.m.thawev == doctest::Approx(0.16 * (r4 + r5)).epsilon(1e-12)); CHECK(s3.m.thawex == doctest::Approx(0.015 * (r4 + r5)).epsilon(1e-12));
  CHECK(s3.m.export_ - s0.m.export_ == doctest::Approx(s3.m.thawex).epsilon(1e-9));
  CHECK(s3.m.evasion > s0.m.evasion); CHECK(s3.m.soilev - s0.m.soilev == doctest::Approx(s3.m.evasion - s0.m.evasion).epsilon(1e-12));
  for (int l = 0; l < 4; l++) { CHECK(s3.soil[l] == s0.soil[l]); CHECK(s3.diss[l] == s0.diss[l]); }
  CHECK(s3.diss[4] < s0.diss[4]);
  CHECK(std::fabs(s3.balance_error()) < 1e-13 * s3.total());
  // frozen (T <= 0) layers below zpf: frozen-soil decomposition re-sorbs in place as with the switch off (no thaw release)
  for (int l = 4; l < 6; l++) { f.ts[l] = (l == 4) ? 0.0 : -1.5; }
  hgtr::State s4 = run(q, 0.25), s5 = run(p, 0.25);
  CHECK(s4.m.thawev == 0.0); CHECK(s4.m.thawex == 0.0); CHECK(s4.m.export_ == s5.m.export_); CHECK(s4.m.evasion == s5.m.evasion);
  for (int l = 0; l < 6; l++) { CHECK(s4.soil[l] == s5.soil[l]); CHECK(s4.diss[l] == s5.diss[l]); }
  f.ts[4] = 0.5;   // only layer 4 thawed
  hgtr::State s6 = run(q, 0.25);
  CHECK(s6.m.thawev == doctest::Approx(0.16 * r4).epsilon(1e-12)); CHECK(s6.m.thawex == doctest::Approx(0.015 * r4).epsilon(1e-12));
  CHECK(std::fabs(s6.balance_error()) < 1e-13 * s6.total());
  std::remove(T2DEP);
}

TEST_CASE("hg_t2 phase 6b all on (deep spread + thaw release), random steps with geometry changes and thaw tracking: balance to roundoff, pools >= 0") {
  std::mt19937 rng(67); std::uniform_real_distribution<double> U(0.0, 1.0);
  hgtr::Params p = t2dparams(); p.t2_evasion = 1; p.t2_foliar = 2; p.t2_meltroute = 1; p.fire_ash = 0.2; p.t2_rts = 1; p.wash = 0.5;
  p.t2_kdsplit = 1; p.t2_dlog_kd_hum = 0.4; p.t2_xpath = 1; p.f_conn = 0.3; p.t2_snowsink = 1; p.v_snow = 4e-4;
  p.t2_deeproute = 1; p.f_deep = 0.8; p.z_deep = 0.5; p.t2_deepspread = 1; p.t2_thawrel = 1; p.thaw_ref_year = 2050;
  hgtr::State s; hgtr::Snap a = column(8, 3000.0); hgtr::initialize(p, s, a);
  double thaw = 0.0;
  for (int it = 0; it < 1500; it++) {
    const int mth = it % 12;
    hgtr::Snap b = a; hgtr::Flux f = t2flux(3 * U(rng), U(rng), (mth < 5 || mth > 9) ? 200 * U(rng) : 0.0, 30 * U(rng), 5 * U(rng), 0.4 * U(rng));
    for (int i = 0; i < NUM_PFT; i++) { f.gpp[i] = 50 * U(rng); }
    f.ta = 30 * U(rng) - 15; f.co2 = 300 + 200 * U(rng); f.rain = 60 * U(rng);
    for (int l = 0; l < a.nsl; l++) { f.ts[l] = 20 * U(rng) - 10; f.decomp[l] = 5 * U(rng); b.soil[l] = std::fmax(0.0, a.soil[l] - f.decomp[l] + 4 * U(rng)); }
    if (U(rng) < 0.3) {
      int n = 3 + int(U(rng) * 12); b.nsl = n; double c0 = 0; for (int l = 0; l < a.nsl; l++) c0 += a.soil[l];
      for (int l = 0; l < MAX_SOI_LAY; l++) { b.soil[l] = (l < n) ? c0 / n : 0; b.dz[l] = (l < n) ? 0.05 + 0.1 * U(rng) : 0; b.type[l] = (l == 0) ? 0 : ((l < 2) ? 1 : ((l < 3) ? 2 : 3)); }
    }
    hgtr::step(p, s, a, b, f);
    const int ff = (mth < 5 || mth > 9) ? int(U(rng) * 3) : 3 + int(U(rng) * ((it < 900) ? 3 : 8));   // deeper thaw later
    hgtr::Flux t = dflux(b.nsl, 30, 20 * U(rng), 10 * U(rng), 3 * U(rng), 5 + 40 * U(rng), ff);
    t.year = 2000 + it / 12;
    for (int d = 0; d < 30; d++) t.d_melt[d] = (U(rng) < 0.3) ? 5 * U(rng) : 0.0;
    hgtr::thaw_track(p, s, b, t);
    hgtr::step(p, s, b, b, t);
    a = b; thaw += s.m.thawev + s.m.thawex; s.m.clear();
  }
  CHECK(std::fabs(s.balance_error()) < 1e-8 * std::fmax(1.0, s.total()));
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s.soil[l] >= 0.0); CHECK(s.diss[l] >= 0.0); }
  CHECK(s.zmax >= s.zpf); CHECK(s.zpf > 0.0); CHECK(s.cum.export_ > 0.0); CHECK(thaw > 0.0);
  std::remove(T2DEP);
}

// ---- Phase M1 (MeHg literature items, docs_src/WORK_NOTES_hg_t2.md section 13): substrate switch, rate and export diagnostics
TEST_CASE("hg_t2_mehg_subst parse: default 0, range 0-2, needs hg_mehg, hg_t2 and hg_t2_diss; hg_q10_* usable with hg_t2") {
  hgtr::Params p = t2params(); CHECK(p.t2_mehg_subst == 0);   // writes T2DEP
  Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_t2"] = 1; m["hg_dep_file"] = T2DEP; m["hg_mehg"] = 1;
  hgtr::Params r; REQUIRE(hgtr::parse(m, r) == ""); CHECK(r.t2_mehg_subst == 0); CHECK(r.q10_meth == 1.0); CHECK(r.q10_demeth == 2.0);
  Json::Value a = m; a["hg_t2_mehg_subst"] = 2; a["hg_q10_meth"] = 3.0; a["hg_q10_demeth"] = 1.0;
  REQUIRE(hgtr::parse(a, r) == ""); CHECK(r.t2_mehg_subst == 2); CHECK(r.q10_meth == 3.0); CHECK(r.q10_demeth == 1.0);
  Json::Value b = m; b["hg_t2_mehg_subst"] = 3; CHECK(hgtr::parse(b, r) != "");
  Json::Value c = m; c["hg_t2_mehg_subst"] = 1; c["hg_mehg"] = 0; CHECK(hgtr::parse(c, r) != "");   // needs the index
  Json::Value d = m; d["hg_t2_mehg_subst"] = 1; d["hg_t2_diss"] = 0; CHECK(hgtr::parse(d, r) != "");   // needs the dissolved pool
  Json::Value e = ms0(); e["hg_tracer"] = 1; e["hg_mehg"] = 1; e["hg_t2_mehg_subst"] = 1; CHECK(hgtr::parse(e, r) != "");   // needs hg_t2
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_mehg_subst: substrate = layer Hg (0, unchanged), dissolved Hg (1), Kd-equilibrium share (2); rate diagnostics") {
  hgtr::Params p; p.tracer = 1; p.mehg = 1; p.k_meth_proxy = 0.6; p.k_meth_ch4 = 30.0;
  hgtr::Snap a = column(4); hgtr::State s0; hgtr::initialize(p, s0, a);
  for (int l = 0; l < 4; l++) { s0.diss[l] = 0.1 * (l + 1) * s0.soil[l]; s0.phieq[l] = 0.05 * (l + 1); }
  double A[MAX_SOI_LAY] = {}, T[MAX_SOI_LAY];
  for (int l = 0; l < MAX_SOI_LAY; l++) { A[l] = 1.0e-3 * (l + 1); T[l] = p.t_ref; }
  for (int k = 0; k < 3; k++) {
    hgtr::Params q = p; q.t2_mehg_subst = k; hgtr::State s = s0;
    for (int m = 0; m < 60; m++) { hgtr::mehg_month(q, s, A, T, 4, true); }
    for (int l = 0; l < 4; l++) {
      const double hl = s.soil[l] + s.diss[l], hs = (k == 0) ? hl : (k == 1) ? s.diss[l] : hl * s.phieq[l];
      CHECK(s.mehg[l] == doctest::Approx(q.k_meth_ch4 * A[l] * hs / q.k_demeth).epsilon(1e-6));
    }
    CHECK(s.total() == s0.total());   // the index never changes Hg
    // one month from zeroed counters: k_m,eff = sum P / sum Hg / 30.44, kd_eff = k_demeth / 30.44 at t_ref
    s.mehg_prod_m = 0; s.mehg_hsub_m = 0; s.mehg_kdp_m = 0; hgtr::mehg_month(q, s, A, T, 4, true);
    double P = 0, Hs = 0; for (int l = 0; l < 4; l++) { const double hl = s.soil[l] + s.diss[l]; const double hs = (k == 0) ? hl : (k == 1) ? s.diss[l] : hl * s.phieq[l]; P += q.k_meth_ch4 * A[l] * hs; Hs += hl; }
    CHECK(s.mehg_prod_m == doctest::Approx(P).epsilon(1e-12));
    CHECK(hgtr::mehg_kmeff_d(s.mehg_prod_m, s.mehg_hsub_m) == doctest::Approx(P / Hs / 30.44).epsilon(1e-12));
    CHECK(hgtr::mehg_kdeff_d(s.mehg_kdp_m, s.mehg_prod_m) == doctest::Approx(q.k_demeth / 30.44).epsilon(1e-12));
  }
  CHECK(hgtr::mehg_kmeff_d(1.0, 0.0) == 0.0); CHECK(hgtr::mehg_kdeff_d(1.0, 0.0) == 0.0);
  // q10 swap (Q10_meth 3, Q10_demeth 1) at t_ref - 10: production / 3, demethylation unchanged
  hgtr::Params q = p; q.q10_meth = 3.0; q.q10_demeth = 1.0; hgtr::State s = s0;
  for (int l = 0; l < MAX_SOI_LAY; l++) { T[l] = p.t_ref - 10.0; }
  for (int m = 0; m < 60; m++) { hgtr::mehg_month(q, s, A, T, 4, false); }
  CHECK(s.mehg[1] == doctest::Approx(q.k_meth_proxy * A[1] * (s.soil[1] + s.diss[1]) / 3.0 / q.k_demeth).epsilon(1e-6));
}

TEST_CASE("HG_MEHGEXP diagnostic: MeHg share x dissolved export, zero without flow or index; Hg pools and the index unchanged") {
  hgtr::Params p = t2dparams(); p.d_cryo = 0; p.d_bio = 0; p.mehg = 1;
  hgtr::Snap a = column(6); hgtr::Flux f = dflux(6, 10, 2.0, 1.0, 1.0, 20.0, 5);
  hgtr::State s; hgtr::initialize(p, s, a);
  for (int l = 0; l < 6; l++) { s.diss[l] = 0.01 * s.soil[l]; s.soil[l] -= s.diss[l]; s.mehg[l] = 0.01 * (s.soil[l] + s.diss[l]); }
  hgtr::State s1 = s, s2 = s; hgtr::Params p0 = p; p0.mehg = 0;
  hgtr::step(p, s1, a, a, f); hgtr::step(p0, s2, a, a, f);
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s1.soil[l] == s2.soil[l]); CHECK(s1.diss[l] == s2.diss[l]); }
  CHECK(s1.m.export_ == s2.m.export_); CHECK(s1.m.export_ > 0.0);
  CHECK(s2.mehg_exp_m == 0.0);   // index off: no diagnostic
  CHECK(s1.mehg_exp_m >= 0.01 * s1.m.export_ * (1.0 - 1e-9)); CHECK(s1.mehg_exp_m <= 0.0105 * s1.m.export_);   // share 1 % (grows as Hg leaves)
  CHECK(s1.mehg_exp_y == s1.mehg_exp_m);
  hgtr::State s3 = s; hgtr::Flux g = dflux(6, 10, 0.0, 0.0, 0.0, 20.0, 5); hgtr::step(p, s3, a, a, g);
  CHECK(s3.mehg_exp_m == 0.0);
  // subst 2: monthly Kd-equilibrium share in [0, 1], unset (0) with subst 0
  hgtr::State s4 = s; hgtr::Params p2 = p; p2.t2_mehg_subst = 2; hgtr::step(p2, s4, a, a, f);
  for (int l = 0; l < 6; l++) { CHECK(s4.phieq[l] > 0.0); CHECK(s4.phieq[l] <= 1.0); CHECK(s1.phieq[l] == 0.0); }
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s4.soil[l] == s1.soil[l]); CHECK(s4.diss[l] == s1.diss[l]); }
  std::remove(T2DEP);
}

// ---- Phase M2 item 1 (WORK_NOTES_hg_t2.md section 14): frozen-layer methylation freeze factor, MEHGPRODL
TEST_CASE("hg_t2_mehg_frz parse: default 0 (t0 0, dt 1), range 0-2, needs hg_mehg and hg_t2") {
  hgtr::Params p = t2params(); CHECK(p.t2_mehg_frz == 0); CHECK(p.t2_mehg_frz_t0 == 0.0); CHECK(p.t2_mehg_frz_dt == 1.0);
  Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_t2"] = 1; m["hg_dep_file"] = T2DEP; m["hg_mehg"] = 1;
  hgtr::Params r;
  Json::Value a = m; a["hg_t2_mehg_frz"] = 2; a["hg_t2_mehg_frz_t0"] = -1.0; a["hg_t2_mehg_frz_dt"] = 2.0;
  REQUIRE(hgtr::parse(a, r) == ""); CHECK(r.t2_mehg_frz == 2); CHECK(r.t2_mehg_frz_t0 == -1.0); CHECK(r.t2_mehg_frz_dt == 2.0);
  Json::Value b = m; b["hg_t2_mehg_frz"] = 3; CHECK(hgtr::parse(b, r) != "");
  Json::Value c = m; c["hg_t2_mehg_frz"] = 1; c["hg_mehg"] = 0; CHECK(hgtr::parse(c, r) != "");
  Json::Value d = ms0(); d["hg_tracer"] = 1; d["hg_mehg"] = 1; d["hg_t2_mehg_frz"] = 1; CHECK(hgtr::parse(d, r) != "");
  Json::Value e = m; e["hg_t2_mehg_frz_dt"] = 0.0; CHECK(hgtr::parse(e, r) != "");
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_mehg_frz: factor algebra; off bit-identical; cutoff stops production at T <= t0 (MeHg decays); MEHGPRODL sums to MEHGPROD") {
  hgtr::Params p; p.tracer = 1; p.mehg = 1; p.k_meth_proxy = 0.6;
  hgtr::Params q1 = p; q1.t2_mehg_frz = 1; hgtr::Params q2 = p; q2.t2_mehg_frz = 2; q2.t2_mehg_frz_t0 = -1.0; q2.t2_mehg_frz_dt = 2.0;
  CHECK(hgtr::mehg_frz_factor(p, -30.0) == 1.0);
  CHECK(hgtr::mehg_frz_factor(q1, 0.0) == 0.0); CHECK(hgtr::mehg_frz_factor(q1, 1e-9) == 1.0); CHECK(hgtr::mehg_frz_factor(q1, -5.0) == 0.0);
  CHECK(hgtr::mehg_frz_factor(q2, -1.0) == 0.0); CHECK(hgtr::mehg_frz_factor(q2, 0.0) == doctest::Approx(0.5)); CHECK(hgtr::mehg_frz_factor(q2, 1.0) == 1.0);
  CHECK(hgtr::mehg_frz_factor(q2, -3.0) == 0.0); CHECK(hgtr::mehg_frz_factor(q2, 5.0) == 1.0);
  hgtr::Snap a = column(4); hgtr::State s0; hgtr::initialize(p, s0, a);
  for (int l = 0; l < 4; l++) { s0.mehg[l] = 0.002 * s0.soil[l]; }
  double A[MAX_SOI_LAY] = {}, T[MAX_SOI_LAY] = {};
  for (int l = 0; l < MAX_SOI_LAY; l++) { A[l] = 1.0e-3 * (l + 1); }
  T[0] = -4.0; T[1] = 0.0; T[2] = 0.5; T[3] = 6.0;
  hgtr::State s = s0, s1 = s0, s2 = s0; hgtr::Params p0 = p; p0.t2_mehg_frz_t0 = 3.0; p0.t2_mehg_frz_dt = 7.0;   // mode 0 ignores t0 / dt
  hgtr::mehg_month(p, s, A, T, 4, false); hgtr::State sx = s0; hgtr::mehg_month(p0, sx, A, T, 4, false);
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(sx.mehg[l] == s.mehg[l]); CHECK(sx.mehg_prodl[l] == s.mehg_prodl[l]); }
  CHECK(sx.mehg_prod_m == s.mehg_prod_m);
  hgtr::mehg_month(q1, s1, A, T, 4, false); hgtr::mehg_month(q2, s2, A, T, 4, false);
  for (int l = 0; l < 4; l++) {
    const double kd = p.k_demeth * std::pow(p.q10_demeth, (T[l] - p.t_ref) / 10.0);
    const double f2 = std::fmin(1.0, std::fmax(0.0, (T[l] + 1.0) / 2.0));
    CHECK(s1.mehg_prodl[l] == doctest::Approx(((T[l] > 0.0) ? 1.0 : 0.0) * s.mehg_prodl[l]).epsilon(1e-14));
    CHECK(s2.mehg_prodl[l] == doctest::Approx(f2 * s.mehg_prodl[l]).epsilon(1e-14));
    if (T[l] <= 0.0) { CHECK(s1.mehg[l] == doctest::Approx(s0.mehg[l] * std::exp(-kd)).epsilon(1e-12)); }   // pure decay
    else { CHECK(s1.mehg[l] == s.mehg[l]); }
  }
  double sp = 0, s1p = 0; for (int l = 0; l < MAX_SOI_LAY; l++) { sp += s.mehg_prodl[l]; s1p += s1.mehg_prodl[l]; }
  CHECK(sp == doctest::Approx(s.mehg_prod_m).epsilon(1e-14)); CHECK(s1p == doctest::Approx(s1.mehg_prod_m).epsilon(1e-14));
  CHECK(s.mehg_prodl[0] > 0.0); CHECK(s1.total() == s0.total());   // index never changes Hg
}

// ---- Phase M2 items 2 and h (WORK_NOTES_hg_t2.md section 14): MeHg transport, depth-attenuated thaw evasion
TEST_CASE("hg_t2_mehg_transport / mobfac / hg_t2_thaw_evz parse: defaults off, dependencies") {
  hgtr::Params p = t2params(); CHECK(p.t2_mehg_transport == 0); CHECK(p.t2_mehg_mobfac == 1.0); CHECK(p.t2_thaw_evz == 0.0);
  Json::Value m = ms0(); m["hg_tracer"] = 1; m["hg_t2"] = 1; m["hg_dep_file"] = T2DEP; m["hg_mehg"] = 1;
  hgtr::Params r;
  Json::Value a = m; a["hg_t2_mehg_transport"] = 1; a["hg_t2_mehg_mobfac"] = 3.0; a["hg_t2_thaw_evz"] = 0.1;
  REQUIRE(hgtr::parse(a, r) == ""); CHECK(r.t2_mehg_transport == 1); CHECK(r.t2_mehg_mobfac == 3.0); CHECK(r.t2_thaw_evz == 0.1);   // preset: thawrel 1
  Json::Value b = m; b["hg_t2_mehg_transport"] = 2; CHECK(hgtr::parse(b, r) != "");
  Json::Value c = m; c["hg_t2_mehg_transport"] = 1; c["hg_mehg"] = 0; CHECK(hgtr::parse(c, r) != "");
  Json::Value d = m; d["hg_t2_mehg_transport"] = 1; d["hg_t2_diss"] = 0; CHECK(hgtr::parse(d, r) != "");
  Json::Value e = m; e["hg_t2_thaw_evz"] = 0.1; e["hg_t2_thawrel"] = 0; CHECK(hgtr::parse(e, r) != "");
  Json::Value g = m; g["hg_t2_thaw_evz"] = -1.0; CHECK(hgtr::parse(g, r) != "");
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_mehg_transport: Hg pools identical to off; MeHg conserved (start - end = removal + clamp); no flow == off; export <= M1 diagnostic") {
  hgtr::Params p = t2dparams(); p.d_cryo = 0; p.d_bio = 0; p.mehg = 1;
  hgtr::Snap a = column(6); hgtr::Flux f = dflux(6, 10, 2.0, 1.0, 1.0, 20.0, 5);
  hgtr::State s; hgtr::initialize(p, s, a);
  for (int l = 0; l < 6; l++) { s.diss[l] = 0.01 * s.soil[l]; s.soil[l] -= s.diss[l]; s.mehg[l] = 0.01 * (l + 1) * (s.soil[l] + s.diss[l]); }
  hgtr::Params q = p; q.t2_mehg_transport = 1; hgtr::Params q3 = q; q3.t2_mehg_mobfac = 3.0; hgtr::Params p0m = p; p0m.t2_mehg_mobfac = 7.0;
  hgtr::State s0 = s, s1 = s, s3 = s, sm = s;
  hgtr::step(p, s0, a, a, f); hgtr::step(q, s1, a, a, f); hgtr::step(q3, s3, a, a, f); hgtr::step(p0m, sm, a, a, f);
  for (int l = 0; l < MAX_SOI_LAY; l++) {
    CHECK(s1.soil[l] == s0.soil[l]); CHECK(s1.diss[l] == s0.diss[l]); CHECK(s3.diss[l] == s0.diss[l]);
    CHECK(sm.mehg[l] == s0.mehg[l]);   // mobfac ignored when off
  }
  CHECK(s1.m.export_ == s0.m.export_); CHECK(s1.balance_error() == s0.balance_error());
  CHECK(s0.mehg_tbal_m == 0.0); CHECK(sm.mehg_exp_m == s0.mehg_exp_m);
  CHECK(s1.mehg_exp_m > 0.0); CHECK(s1.mehg_exp_m <= s0.mehg_exp_m * (1.0 + 1e-12));
  CHECK(std::fabs(s1.mehg_tbal_m) < 1e-12 * s.mehgtot()); CHECK(std::fabs(s3.mehg_tbal_m) < 1e-12 * s.mehgtot());
  CHECK(s3.mehg_exp_m > s1.mehg_exp_m);   // more mobile MeHg -> more export
  CHECK(s1.mehgtot() < s0.mehgtot());
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(s1.mehg[l] >= 0.0); CHECK(s1.mehg[l] <= s1.soil[l] + s1.diss[l]); }
  hgtr::Flux g = dflux(6, 10, 0.0, 0.0, 0.0, 20.0, 5); hgtr::State t0 = s, t1 = s; hgtr::step(p, t0, a, a, g); hgtr::step(q, t1, a, a, g);
  for (int l = 0; l < MAX_SOI_LAY; l++) { CHECK(t1.mehg[l] == t0.mehg[l]); }
  CHECK(t1.mehg_exp_m == 0.0);
  std::remove(T2DEP);
}

TEST_CASE("hg_t2_thaw_evz: thaw evasion fraction f_thaw_evade x exp(-(z_mid - zpf) / evz); 0 == constant (bit-identical); release conserved") {
  hgtr::Params p = t2dparams(); p.t2_evasion = 1; p.t2_thawrel = 1;
  CHECK(hgtr::thaw_fev(p, 3.0, 0.2) == p.f_thaw_evade);
  hgtr::Params q = p; q.t2_thaw_evz = 0.1;
  CHECK(hgtr::thaw_fev(q, 0.3, 0.25) == doctest::Approx(0.16 * std::exp(-0.5)).epsilon(1e-14)); CHECK(hgtr::thaw_fev(q, 0.2, 0.25) == 0.16);
  hgtr::Snap a = column(6);
  hgtr::Flux f = t2flux(1.0, 0.0, 0.0, 0.0, 0.0, 0.2);
  for (int l = 0; l < 6; l++) { f.decomp[l] = 5.0; f.ts[l] = 2.0; }
  hgtr::Snap b = a; for (int l = 0; l < 6; l++) { b.soil[l] = a.soil[l] - 5.0; }
  auto run = [&](const hgtr::Params& x) { hgtr::State s; hgtr::initialize(x, s, a); s.zpf = 0.25; s.h_init = s.total(); hgtr::step(x, s, a, b, f); return s; };
  hgtr::State s0 = run(p), s1 = run(q);
  double z = 0.0, w = 0.0; for (int l = 0; l < 6; l++) { const double d = a.dz[l]; if (z >= 0.25) { w += 5.0 * p.r_min * std::exp(-std::fmax(0.0, z + 0.5 * d - 0.25) / 0.1); } z += d; }
  CHECK(s1.m.thawev == doctest::Approx(0.16 * w).epsilon(1e-12)); CHECK(s1.m.thawev < s0.m.thawev);
  CHECK(s1.m.thawex == doctest::Approx(s0.m.thawex).epsilon(1e-14));
  CHECK(std::fabs(s1.balance_error()) < 1e-13 * s1.total());
  CHECK(s1.total() > s0.total());   // less evaded -> more retained
  std::remove(T2DEP);
}
