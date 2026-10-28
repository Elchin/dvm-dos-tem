// v5q-negpool: switchable fix NP (default OFF = legacy).
//   model_settings.soil_negpool_fix -> Integrator::soil_negpool_fix
// Legacy y2cstate_soi() clamps a negative layer AVLN with ORGN += AVLN - 0.001,
// AVLN = 0.001, which writes ORGN = -0.001 into an empty moss layer; the next
// soil month then cannot take any RKF45 step (checkPools) and gives up at t0 = 0
// for the rest of the run. F1 (clampLayerN) borrows only N the layer has; F2
// (healNegative) zeroes negative pools at month start; both return the removed N
// so it is booked as debt (L1 n_debt) and total N is conserved.
// Refs (NP): Shampine et al. 2005, doi:10.1016/j.amc.2004.12.011.
#include "doctest.h"
#include "test_helpers.h"
#include <cmath>

TEST_CASE("NP default is off (legacy)") {
  CHECK(Integrator::soil_negpool_fix == 0);
}

TEST_CASE("NP F1: non-negative AVLN untouched") {
  float a = 0.25f, o = 3.0f;
  CHECK(Integrator::clampLayerN(a, o) == 0.0);
  CHECK(a == 0.25f); CHECK(o == 3.0f);
}

TEST_CASE("NP F1: legacy arithmetic when the layer can afford it") {
  float a = -0.0004f, o = 2.0f;
  float la = a, lo = o;                       // legacy reference
  lo += la - 0.001; la = 0.001;
  CHECK(Integrator::clampLayerN(a, o) == 0.0);
  CHECK(a == la); CHECK(o == lo);
}

TEST_CASE("NP F1: empty moss layer (CAP125 trigger) gives no negative ORGN") {
  float a = -1.0e-9f, o = 0.0f;
  float lo = o + a - 0.001;                   // legacy would write -0.001
  CHECK(lo < 0.0f);
  double debt = Integrator::clampLayerN(a, o);
  CHECK(a == 0.0f); CHECK(o == 0.0f);
  CHECK(debt == doctest::Approx(1.0e-9).epsilon(1e-6));
}

TEST_CASE("NP F1: partial ORGN, conserved without debt") {
  float a = -0.0005f, o = 0.0008f;
  double tot = (double)a + (double)o;
  double debt = Integrator::clampLayerN(a, o);
  CHECK(debt == 0.0);
  CHECK(a == 0.0f);
  CHECK((double)o == doctest::Approx(tot).epsilon(1e-6));
  CHECK(o >= 0.0f);
}

TEST_CASE("NP F2: healNegative zeroes negatives and returns the deficit") {
  float p[4] = {-0.001f, 0.5f, 0.0f, -0.002f};
  double d = Integrator::healNegative(p, 4);
  CHECK(d == doctest::Approx(0.003).epsilon(1e-6));
  CHECK(p[0] == 0.0f); CHECK(p[1] == 0.5f); CHECK(p[2] == 0.0f); CHECK(p[3] == 0.0f);
  float q[2] = {0.1f, 0.2f};
  CHECK(Integrator::healNegative(q, 2) == 0.0);
}

TEST_CASE("NP F1+F2: N mass conserved incl. debt over random layers") {
  unsigned s = 12345u; double before = 0, after = 0, debt = 0;
  for (int k = 0; k < 1000; ++k) {
    s = s * 1103515245u + 12345u; float a = ((int)(s >> 16) % 2001 - 1000) * 1e-6f;
    s = s * 1103515245u + 12345u; float o = ((int)(s >> 16) % 1001) * 1e-6f;
    before += (double)a + (double)o;
    debt += Integrator::clampLayerN(a, o);
    CHECK(a >= 0.0f); CHECK(o >= 0.0f);
    after += (double)a + (double)o;
  }
  CHECK(after - debt == doctest::Approx(before).epsilon(1e-5));
}
