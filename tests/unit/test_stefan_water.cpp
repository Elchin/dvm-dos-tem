// Area 6: Stefan front handling (all-fronts-swept, rebuild when deque empty);
// Area 8: water bookkeeping (freeze-expansion clip -> expelled water, puddle kept);
// Area 7: Richards array sizing and fail-fast guards.
#include "doctest.h"
#include "test_helpers.h"

using ut::make_soil_layer;

namespace {
// 5 x 0.1 m mineral layers; per-layer (liq fraction of pores, ice fraction, tem, frozen)
struct Spec { double fl, fi, tem; int frozen; double ff; };
std::vector<SoilLayer*> build(const std::vector<Spec>& s, double poro = 0.4, double dz = 0.1) {
  std::vector<SoilLayer*> v;
  for (size_t i = 0; i < s.size(); ++i) {
    v.push_back(make_soil_layer(i * dz, dz, poro, s[i].fl * poro * DENLIQ * dz,
                                s[i].fi * poro * DENICE * dz, s[i].tem, s[i].frozen, s[i].ff));
  }
  return v;
}
}

TEST_SUITE("stefan") {

TEST_CASE("rebuildFrontsFromLayerStates: fronts at every frozen/unfrozen boundary") {
  // unfrozen, unfrozen, frozen, frozen, unfrozen(talik below?), frozen
  ut::Column c(build({{0.6,0,2,-1,0},{0.6,0,1,-1,0},{0,0.6,-1,1,1},{0,0.6,-2,1,1},{0.6,0,0.5,-1,0},{0,0.6,-1,1,1}}));
  Stefan st; st.setGround(&c.g);
  st.rebuildFrontsFromLayerStates();
  REQUIRE(c.g.frontsz.size() == 3);
  CHECK(c.g.frontsz[0] == doctest::Approx(0.2)); CHECK(c.g.frontstype[0] == -1);   // thawing front: unfrozen above
  CHECK(c.g.frontsz[1] == doctest::Approx(0.4)); CHECK(c.g.frontstype[1] == 1);    // freezing front: frozen above
  CHECK(c.g.frontsz[2] == doctest::Approx(0.5)); CHECK(c.g.frontstype[2] == -1);
  // fronts alternate in type
  for (size_t i = 1; i < c.g.frontstype.size(); ++i) CHECK(c.g.frontstype[i] != c.g.frontstype[i-1]);
}

TEST_CASE("rebuildFrontsFromLayerStates: partially frozen layers classified by frozenfrac; uniform column -> no fronts") {
  ut::Column a(build({{0.6,0,1,-1,0},{0.3,0.3,0,0,0.7},{0,0.6,-1,1,1}}));
  Stefan s1; s1.setGround(&a.g); s1.rebuildFrontsFromLayerStates();
  REQUIRE(a.g.frontsz.size() == 1);
  CHECK(a.g.frontsz[0] == doctest::Approx(0.1)); CHECK(a.g.frontstype[0] == -1);
  ut::Column b(build({{0,0.6,-1,1,1},{0,0.6,-1,1,1},{0,0.6,-3,1,1}}));
  Stefan s2; s2.setGround(&b.g); s2.rebuildFrontsFromLayerStates();
  CHECK(b.g.frontsz.empty());
}

TEST_CASE("rebuilt fronts reproduce the layer states (no phase change)") {
  ut::Column c(build({{0.6,0,2,-1,0},{0.6,0,1,-1,0},{0,0.6,-1,1,1},{0,0.6,-2,1,1},{0,0.6,-2,1,1}}));
  Stefan st; st.setGround(&c.g);
  st.rebuildFrontsFromLayerStates();
  double w0 = c.total_water();
  st.all_fronts_swept = false;
  st.updateLayerFrozenState(c.g.toplayer, 1);
  // Frozen fractions (hence ice/liq) are reproduced exactly. The layer whose BOTTOM
  // boundary carries the rebuilt front gets frozen = 0 (a front at z+dz counts as
  // "within" the layer in updateLayerFrozenState: 0 < dz <= layer dz) with frozenfrac 0;
  // this is the same convention as for any front lying exactly on a layer boundary.
  int exp_frozen[5] = {-1, 0, 1, 1, 1};
  double exp_ff[5] = {0, 0, 1, 1, 1};
  for (int i = 0; i < 5; ++i) { CAPTURE(i); CHECK(c.L[i]->frozen == exp_frozen[i]); CHECK(c.L[i]->frozenfrac == doctest::Approx(exp_ff[i])); }
  double ice0[5]; for (int i = 0; i < 5; ++i) ice0[i] = c.L[i]->ice;
  st.updateWaterAfterFront(c.g.toplayer);
  for (int i = 0; i < 5; ++i) { CAPTURE(i); CHECK(c.L[i]->ice == doctest::Approx(ice0[i])); }
  CHECK(c.total_water() == doctest::Approx(w0));
}

TEST_CASE("updateFronts with an empty deque does not thaw deep permafrost (fronts rebuilt)") {
  // active layer thawed, permafrost below; deque empty (e.g. after a sweep). Mild warm day.
  ut::Column c(build({{0.5,0,3,-1,0},{0.5,0,2,-1,0},{0,0.9,-1,1,1},{0,0.9,-2,1,1},{0,0.9,-3,1,1}}));
  Stefan st; st.setGround(&c.g); st.initpce();
  REQUIRE(c.g.frontsz.empty());
  double w0 = c.total_water();
  st.updateFronts(2.0, 86400.);
  // a mild day thaws at most a few cm into the top permafrost layer; the layers below
  // stay fully frozen (without the rebuild, the single new front would flip them all)
  CHECK(c.L[2]->frozenfrac > 0.9); CHECK(c.L[2]->ice > 0.0);
  for (int i = 3; i < 5; ++i) { CAPTURE(i); CHECK(c.L[i]->frozen == 1); CHECK(c.L[i]->frozenfrac == 1.0); }
  REQUIRE(!c.g.frontsz.empty());
  CHECK(c.g.frontstype.back() == -1);                // thawing front on top of permafrost
  CHECK(c.g.frontsz.back() > 0.2); CHECK(c.g.frontsz.back() < 0.21);
  CHECK(c.total_water() == doctest::Approx(w0));
}

TEST_CASE("all fronts swept by a freezing front: every soil layer frozen, water conserved") {
  // talik (0-0.2 m) above permafrost; one thawing front at 0.2 m. The unfrozen layers
  // are still warm (tem > 0) and one permafrost layer is at +0.01 C (slightly above 0).
  ut::Column c(build({{0.5,0,0.5,-1,0},{0.5,0,0.3,-1,0},{0,0.6,-0.5,1,1},{0,0.6,0.01,1,1},{0,0.6,-1,1,1}}));
  c.g.frontsz.push_back(0.2); c.g.frontstype.push_back(-1);
  Stefan st; st.setGround(&c.g); st.initpce();
  double w0 = c.total_water();
  c.g.expelled_water_mm = 0.0;
  st.updateFronts(-100.0, 86400.);                   // very strong forcing: freezes the whole talik in one step
  CHECK(st.all_fronts_swept);
  CHECK(c.g.frontsz.empty());
  for (int i = 0; i < 5; ++i) {
    CAPTURE(i);
    CHECK(c.L[i]->frozen == 1);                      // without the fix layers 0,1,3 would be -1
    CHECK(c.L[i]->frozenfrac == doctest::Approx(1.0));
    CHECK(c.L[i]->liq == doctest::Approx(0.0));
  }
  CHECK(c.total_water() + c.g.expelled_water_mm == doctest::Approx(w0).epsilon(1e-12));
  // phase-change energy recorded for the freezing layers only
  CHECK(c.L[0]->pce_f > 0.0); CHECK(c.L[1]->pce_f > 0.0); CHECK(c.L[2]->pce_f == 0.0);
}

TEST_CASE("no fronts at all and flag unset: layer states from temperature sign (legacy path)") {
  ut::Column c(build({{0.5,0,0.5,-1,0},{0,0.6,-0.5,1,1}}));
  Stefan st; st.setGround(&c.g);
  st.all_fronts_swept = false;
  c.L[1]->tem = 0.2;                                 // frozen layer slightly warm, no fronts
  st.updateLayerFrozenState(c.g.toplayer, 1);
  CHECK(c.L[0]->frozen == -1);
  CHECK(c.L[1]->frozen == -1);                       // documents the legacy flip that the flag prevents
}

TEST_CASE("combineExtraFronts: one surplus pair (the daily case, MAX_NUM_FNT + 1 fronts) is merged") {
  ut::Column c(build({{0.5,0,1,-1,0},{0.5,0,1,-1,0},{0.5,0,1,-1,0},{0.5,0,1,-1,0},{0.5,0,1,-1,0}}));
  for (int i = 0; i < MAX_NUM_FNT + 1; ++i) { c.g.frontsz.push_back(0.01 + 0.03 * i + (i == 4 ? 0.02 : 0.0)); c.g.frontstype.push_back(i % 2 ? 1 : -1); }
  Stefan st; st.setGround(&c.g); st.combineExtraFronts();
  REQUIRE((int)c.g.frontsz.size() == MAX_NUM_FNT - 1);
  // the closest pair (fronts 4 and 5, gap 0.01 m) was removed
  CHECK(c.g.frontsz[4] == doctest::Approx(0.01 + 0.03 * 6));
  for (size_t i = 1; i < c.g.frontstype.size(); ++i) { CHECK(c.g.frontstype[i] != c.g.frontstype[i-1]); CHECK(c.g.frontsz[i] > c.g.frontsz[i-1]); }
}

// KNOWN BUG (reported, not fixed): Stefan::combineExtraFronts initialises 'mindz' once,
// outside its while loop. After the first pair is removed, a later pass finds no gap
// strictly smaller than the previous minimum, removefntid stays -1, and the loop never
// ends. Unreachable from the daily front update (at most one new front per call ->
// at most MAX_NUM_FNT+1 fronts), but reachable via rebuildFrontsFromLayerStates() if a
// column has more than MAX_NUM_FNT+2 frozen/unfrozen boundaries. Run in a forked child
// with a 3 s alarm so the hang cannot block the suite.
TEST_CASE("combineExtraFronts terminates with several surplus pairs [KNOWN BUG: xfail]"
          * doctest::should_fail(true)) {
  int st_code = ut::exit_status_of([]() {
    alarm(3);
    ut::Column c(build({{0.5,0,1,-1,0},{0.5,0,1,-1,0},{0.5,0,1,-1,0},{0.5,0,1,-1,0},{0.5,0,1,-1,0}}));
    for (int i = 0; i < MAX_NUM_FNT + 4; ++i) { c.g.frontsz.push_back(0.01 + 0.03 * i + (i % 3) * 0.001); c.g.frontstype.push_back(i % 2 ? 1 : -1); }
    Stefan st; st.setGround(&c.g); st.combineExtraFronts();
    _exit((int)c.g.frontsz.size() <= MAX_NUM_FNT ? 0 : 1);
  });
  CHECK(st_code == 0);   // currently 128+SIGALRM (hang)
}

} // TEST_SUITE stefan

TEST_SUITE("water") {

TEST_CASE("freeze-expansion clip: excess over pore space is booked as expelled water (mass conserved)") {
  // saturated layer freezes completely: ice would need 1000/917 of the pore volume
  ut::Column c(build({{1.0,0,-1,0,1.0}}));
  Stefan st; st.setGround(&c.g);
  c.L[0]->frozenfrac = 1.0;
  double w0 = c.total_water();
  c.g.expelled_water_mm = 0.0;
  st.updateWaterAfterFront(c.g.toplayer);
  const double poro = 0.4, dz = 0.1;
  CHECK(c.L[0]->ice == doctest::Approx(poro * DENICE * dz));                     // capped at maxice
  CHECK(c.g.expelled_water_mm == doctest::Approx(poro * dz * (DENLIQ - DENICE)));  // 3.32 mm
  CHECK(c.total_water() + c.g.expelled_water_mm == doctest::Approx(w0).epsilon(1e-12));
  CHECK(c.L[0]->pce_f == doctest::Approx(w0 * LHFUS));
}

TEST_CASE("freeze-expansion clip: partially frozen layer, liquid fraction keeps its pore share") {
  ut::Column c(build({{0.9,0.1,-0.1,0,0.5}}));
  Stefan st; st.setGround(&c.g);
  double w0 = c.total_water();
  c.g.expelled_water_mm = 0.0;
  st.updateWaterAfterFront(c.g.toplayer);
  CHECK(c.total_water() + c.g.expelled_water_mm == doctest::Approx(w0).epsilon(1e-12));
  CHECK(c.L[0]->ice <= c.L[0]->maxice + 1e-12);
  CHECK(c.g.expelled_water_mm >= 0.0);
}

TEST_CASE("no clip and nothing expelled when the pores have room") {
  ut::Column c(build({{0.6,0,-1,0,1.0}}));
  Stefan st; st.setGround(&c.g);
  double w0 = c.total_water();
  c.g.expelled_water_mm = 0.0;
  st.updateWaterAfterFront(c.g.toplayer);
  CHECK(c.g.expelled_water_mm == 0.0);
  CHECK(c.L[0]->ice == doctest::Approx(w0));
}

TEST_CASE("surface puddle is kept across years (grnd_beginOfYear no longer zeroes it)") {
  EnvData ed;
  ed.d_soi2l.magic_puddle = 3.25;
  ed.grnd_beginOfYear();
  CHECK(ed.d_soi2l.magic_puddle == doctest::Approx(3.25));
}

} // TEST_SUITE water

TEST_SUITE("richards") {

TEST_CASE("work arrays are sized MAX_SNW_LAY + MAX_SOI_LAY + 2") {
  CHECK(MAX_RICH_IND == MAX_SNW_LAY + MAX_SOI_LAY + 2);
  Richards r;
  CHECK(sizeof(r.qtrans) / sizeof(double) == (size_t)MAX_RICH_IND);
  CHECK(sizeof(r.layer_drain) / sizeof(double) == (size_t)MAX_RICH_IND);
  CHECK(sizeof(r.percolation) / sizeof(double) == (size_t)MAX_RICH_IND);
  CHECK(sizeof(r.nodemm) / sizeof(double) == (size_t)MAX_RICH_IND);
  CHECK(sizeof(r.tridiag_error) / sizeof(double) == (size_t)MAX_RICH_IND);
  // largest index used: drainind + 1 with all snow layers above a full soil column
  CHECK(MAX_SNW_LAY + MAX_SOI_LAY + 1 < MAX_RICH_IND);
}

TEST_CASE("clearRichardsArrays clears every element, including the last one (off-by-one fix)") {
  Richards r;
  for (int i = 0; i < MAX_RICH_IND; ++i) { r.qtrans[i] = 1.; r.qin[i] = 1.; r.qout[i] = 1.; r.amx[i] = 1.; r.rmx[i] = 1.; r.tridiag_error[i] = 1.; r.percolation[i] = 1.; }
  r.clearRichardsArrays();
  for (int i = 0; i < MAX_RICH_IND; ++i) {
    CAPTURE(i);
    CHECK(r.qtrans[i] == 0.0); CHECK(r.qin[i] == 0.0); CHECK(r.qout[i] == 0.0);
    CHECK(r.amx[i] == 0.0); CHECK(r.rmx[i] == 0.0); CHECK(r.tridiag_error[i] == 0.0);
  }
}

static std::vector<SoilLayer*> wet_column(int n) {
  std::vector<SoilLayer*> v;
  for (int i = 0; i < n; ++i) v.push_back(make_soil_layer(0.1 * i, 0.1, 0.4, 0.5 * 0.4 * DENLIQ * 0.1, 0.0, 2.0, -1, 0.0));
  return v;
}

TEST_CASE("fail-fast guard: top index < 1 throws instead of writing out of bounds") {
  ut::Column c(wet_column(4));
  c.L[0]->indl = 0;                                   // invalid
  Richards r; double trans[MAX_SOI_LAY] = {0};
  CHECK_THROWS_AS(r.update(c.L[0], c.L[3], 0.35, 0.0, 1.0, trans, 0.0, 0.0, 0.0), std::runtime_error);
}

TEST_CASE("fail-fast guard: drainage index + 1 >= MAX_RICH_IND throws") {
  ut::Column c(wet_column(4));
  c.L[3]->indl = MAX_RICH_IND - 1;
  Richards r; double trans[MAX_SOI_LAY] = {0};
  CHECK_THROWS_AS(r.update(c.L[0], c.L[3], 0.35, 0.0, 1.0, trans, 0.0, 0.0, 0.0), std::runtime_error);
  // largest legal drainage index does not trip the guard
  c.L[3]->indl = MAX_RICH_IND - 2;
  bool msg_ok = true;
  try { r.update(c.L[0], c.L[3], 0.35, 0.0, 1.0, trans, 0.0, 0.0, 0.0); }
  catch (const std::runtime_error& e) { msg_ok = std::string(e.what()).find("layer index out of range") == std::string::npos; }
  CHECK(msg_ok);
}

TEST_CASE("fail-fast guard: non-finite or zero substep throws instead of looping forever") {
  for (double dt : {0.0, std::nan(""), -86400.0}) {
    CAPTURE(dt);
    ut::Column c(wet_column(4));
    Richards r; r.delta_t = dt; double trans[MAX_SOI_LAY] = {0};
    CHECK_THROWS_WITH_AS(r.update(c.L[0], c.L[3], 0.35, 0.0, 1.0, trans, 0.0, 0.0, 0.0),
                         doctest::Contains("substep loop cannot finish"), std::runtime_error);
  }
}

// KNOWN BUG (reported, not fixed; pre-dates daily mode): when the drainage depth is
// at or above the top of the Richards domain (e.g. a freezing front at the top of
// the first shallow layer in autumn, Ground::setDrainL sets draindepth = fstfntl->z),
// Richards::update returns early and adds the day's infiltration to the puddle/qover
// of its OWN EnvData member 'Richards::ed', which is never connected to the cohort's
// EnvData (grep: no code reads richards.ed). Soil_Env has already reported the water
// as infiltration (qinfl), so it leaves the system without being booked anywhere.
TEST_CASE("early-return branch: infiltrated water is not lost [KNOWN BUG: xfail]"
          * doctest::should_fail(true)) {
  ut::Column c(wet_column(3));
  Richards r;
  double trans[MAX_SOI_LAY] = {0};
  double w0 = c.total_water();
  r.update(c.L[0], c.L[2], 0.0, 0.0, 1.0, trans, 0.0, 5.0 / SEC_IN_DAY, 0.0);
  // caller-visible water: soil layers + Richards drainage/percolation outputs
  double out = r.qdrain * SEC_IN_DAY;
  CHECK(c.total_water() + out == doctest::Approx(w0 + 5.0));
}

TEST_CASE("early-return branch: documents where the water currently goes (Richards' private EnvData)") {
  ut::Column c(wet_column(3));
  Richards r; r.ed.d_soi2l.magic_puddle = 0.0; r.ed.d_soi2l.qover = 0.0;
  double trans[MAX_SOI_LAY] = {0};
  double w0 = c.total_water();
  r.update(c.L[0], c.L[2], 0.0, 0.0, 1.0, trans, 0.0, 5.0 / SEC_IN_DAY, 0.0);
  CHECK(r.ed.d_soi2l.magic_puddle == doctest::Approx(5.0));
  CHECK(c.total_water() == doctest::Approx(w0));
  CHECK(r.qdrain == 0.0);
}

} // TEST_SUITE richards
