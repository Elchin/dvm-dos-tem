// Lightweight unit tests for RTS initiation helpers (no external framework).
// Build/run: make rts-unit-test
#include <cmath>
#include <cstdio>
#include <vector>
#include "Thermokarst.h"
#include "Ground.h"

static int g_fail = 0, g_pass = 0;
#define CHECK(cond) do { if (cond) { g_pass++; } else { g_fail++; \
  std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)
#define CHECK_NEAR(a, b, tol) CHECK(std::fabs((a) - (b)) <= (tol))

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

static void test_partition_mass_balance() {
  const double fr[] = {0.0, 0.1, 0.5, 1.0, -0.3, 1.7};
  for (double f : fr) {
    double rc, rn, lc, ln;
    Thermokarst::partitionLostCN(5000.0, 175.0, f, rc, rn, lc, ln);
    CHECK_NEAR(rc + lc, 5000.0, 1e-9);   // removed = respired + lateral
    CHECK_NEAR(rn + ln, 175.0, 1e-9);
    CHECK(rc >= 0 && lc >= 0 && rn >= 0 && ln >= 0);
  }
  double rc, rn, lc, ln;
  Thermokarst::partitionLostCN(1000.0, 10.0, 0.25, rc, rn, lc, ln);
  CHECK_NEAR(rc, 250.0, 1e-9); CHECK_NEAR(lc, 750.0, 1e-9);
  CHECK_NEAR(rn, 2.5, 1e-9);   CHECK_NEAR(ln, 7.5, 1e-9);
}

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
  CHECK(thawed[0] && thawed[1] && thawed[2] && !thawed[3]);
  CHECK_NEAR(ice[0] + ice[1] + ice[2], 0.0, 1e-12);
  CHECK_NEAR(ice[3], 80.0, 1e-12);                 // below thaw depth: unchanged
  CHECK_NEAR(liq[1], 100.0, 1e-12);                // capped at maxliq
  CHECK_NEAR(drained, 20.0, 1e-9);                 // 120 ice -> 100 kept + 20 drained
  double total_after = drained; for (size_t i = 0; i < 4; i++) total_after += ice[i] + liq[i];
  CHECK_NEAR(total_after, total_before, 1e-9);     // water conserved incl. drainage
  // partial layer is not thawed (whole layers only)
  std::vector<double> ice2 = {10, 10, 10, 10}, liq2 = {0, 0, 0, 0};
  Ground::thawColumn(ztop, dz, ice2, liq2, maxliq, 0.25, drained, thawed);
  CHECK(thawed[0] && !thawed[1]);
  CHECK_NEAR(ice2[1], 10.0, 1e-12);
  // zero depth: no-op
  std::vector<double> ice3 = {10, 10, 10, 10}, liq3 = {1, 1, 1, 1};
  double m3 = Ground::thawColumn(ztop, dz, ice3, liq3, maxliq, 0.0, drained, thawed);
  CHECK_NEAR(m3, 0.0, 1e-12); CHECK_NEAR(drained, 0.0, 1e-12); CHECK(!thawed[0]);
}

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
  CHECK(dz == dz0);                                  // matrix thickness unchanged
  CHECK_NEAR(col0 - col1, sub, 1e-12);               // column thickness drops by subsidence
  CHECK_NEAR(w1, w0, 1e-9);                          // water mass conserved
  // no excess ice -> no-op
  std::vector<double> e0 = {0, 0, 0}, l0 = {1, 2, 3};
  sub = Ground::meltExcessColumn(e0, l0, maxliq, meltfrac, melted, drained);
  CHECK_NEAR(sub, 0.0, 1e-15); CHECK_NEAR(melted, 0.0, 1e-15); CHECK_NEAR(drained, 0.0, 1e-15);
  CHECK(l0[0] == 1 && l0[1] == 2 && l0[2] == 3);
  // frozen layers (meltfrac 0) -> no-op
  std::vector<double> e1 = {100, 100, 100}, l1 = {1, 1, 1}, f0 = {0, 0, 0};
  sub = Ground::meltExcessColumn(e1, l1, maxliq, f0, melted, drained);
  CHECK_NEAR(sub, 0.0, 1e-15); CHECK(e1[0] == 100 && e1[2] == 100);
}

int main() {
  test_excess_ice_init();
  test_excess_melt_subsidence();
  test_removal_thickness();
  test_partition_mass_balance();
  test_thaw_column();
  std::printf("RTS unit tests: %d passed, %d failed\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}
