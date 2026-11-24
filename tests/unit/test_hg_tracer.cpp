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
  CHECK(p.f_evade == doctest::Approx(0.02)); CHECK(p.f_export == doctest::Approx(0.015)); CHECK(p.f_leach == 0.0);
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
  // geometry mismatch: refused, state unchanged
  hgtr::State q = s; CHECK_FALSE(hgtr::spinup_jump(p, q, yr, c.nsl + 1)); CHECK(q.total() == s.total());
}

TEST_CASE("hg restart layout v2 round trip (spin, spun) and v1 still read") {
  hgtr::Params p; hgtr::State s; hgtr::Snap a = column(6); hgtr::initialize(p, s, a);
  s.cum.spin = 12.5; s.spun = 1; s.mehg[1] = 0.25;
  double v[hgtr::NSTATE]; hgtr::pack(s, v); CHECK(v[0] == 2.0);
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
