// n_safetynet_conserve 2 (default 0 = legacy): L1 plus exact accumulation in
// the soil adapt(). The soil state y[] is float; an accepted step adds dt*f to
// each pool in float, so a large pool (SOMCR ~1e3 g/m2) drops every increment
// below half its float spacing while the small RH / N-flux accumulators keep
// it. When the step size collapses (months that end in a give-up) that is a
// systematic C and N gain (RTS ice case: +0.39 gC/m2 in one month). Level 2
// sums the same float increments in double (Integrator::acceptStep).
// Refs: Higham 2002, doi:10.1137/1.9780898718027 (ch. 4); Kahan 1965,
//   doi:10.1145/363707.363723.
#include "doctest.h"
#include "test_helpers.h"
#include <json/json.h>
#include <cmath>
#define private public            // tk_set_month_flags (private), as test_access.h does
#include "../../include/Cohort.h"
#undef private

namespace {
struct LvGuard { ~LvGuard() { Integrator::n_safetynet_conserve = 0; } };
// pool (index 0) loses f per unit time to an accumulator (index 1); n steps of dt
void run(bool dacc, int n, float dt, float f, float &pool, float &acc) {
  float st[2] = {1500.0f, 0.0f}, d4[2], yp[2] = {-f, f};
  double yd[2] = {st[0], st[1]};
  for (int k = 0; k < n; ++k) {
    for (int i = 0; i < 2; ++i) d4[i] = st[i] + dt * yp[i];   // what rkf45 returns
    Integrator::acceptStep(st, d4, yd, yp, dt, true, dacc, 2);
  }
  pool = st[0]; acc = st[1];
}
}

TEST_SUITE("snc2") {

TEST_CASE("acceptStep legacy (level 0/1) copies the float RKF result unchanged (bit-identical)") {
  float st[3] = {1.f, 2.f, 3.f}, d4[3] = {1.5f, 2.25f, -0.125f}, yp[3] = {9.f, 9.f, 9.f};
  double yd[3] = {7., 7., 7.};
  Integrator::acceptStep(st, d4, yd, yp, 0.5f, true, false, 3);
  for (int i = 0; i < 3; ++i) { CHECK(st[i] == d4[i]); CHECK(yd[i] == 7.); }
}

TEST_CASE("float absorption: legacy pool keeps increments below half its spacing, accumulator books them") {
  float pool, acc;
  const int n = 100000; const float dt = std::ldexp(1.0f, -17), f = 0.01f;   // dt*f ~ 7.6e-8 << 6e-5
  run(false, n, dt, f, pool, acc);
  CHECK(pool == 1500.0f);                                    // pool never changes
  CHECK(acc == doctest::Approx(n * dt * f).epsilon(1e-3));   // but the flux is booked
  CHECK((pool + acc) - 1500.0 > 0.007);                      // -> mass gain ~ 0.0076
}

TEST_CASE("level 2: same increments summed in double; pool + accumulator conserved") {
  float pool, acc;
  const int n = 100000; const float dt = std::ldexp(1.0f, -17), f = 0.01f;
  run(true, n, dt, f, pool, acc);
  CHECK(acc == doctest::Approx(n * dt * f).epsilon(1e-3));
  CHECK(std::fabs((double)pool + acc - 1500.0) < 1.3e-4);    // one float rounding of the pool
  CHECK((double)pool == doctest::Approx(1500.0 - n * dt * f).epsilon(1e-7));
}

TEST_CASE("level 2: a step forced after a negative pool (stepok false) adds nothing") {
  float st[2] = {1.f, 2.f}, d4[2] = {1.f, 2.f}, yp[2] = {100.f, -100.f};   // yprime is partial garbage
  double yd[2] = {1.0, 2.0};
  Integrator::acceptStep(st, d4, yd, yp, 0.25f, false, true, 2);
  CHECK(st[0] == 1.f); CHECK(st[1] == 2.f); CHECK(yd[0] == 1.0); CHECK(yd[1] == 2.0);
}

TEST_CASE("tk_set_month_flags: after an RTS initiation level 2 (unless the key is false); before: global level") {
  LvGuard g; ModelData md; Cohort c; c.md = &md;
  struct Row { int glob; bool tk; int msi; int exp; };
  Row rows[] = {{0, true, -1, 0}, {0, true, 0, 2}, {0, true, 5, 2}, {0, false, 5, 0},
                {1, true, -1, 1}, {1, true, 5, 2}, {1, false, 5, 1}, {2, true, -1, 2}, {2, false, 5, 2}};
  for (const Row &r : rows) {
    md.n_safetynet_conserve = r.glob; md.tk_safetynet_conserve = r.tk; c.ground.tk_months_since_init = r.msi;
    c.tk_set_month_flags();
    CHECK(Integrator::n_safetynet_conserve == r.exp);
  }
}

} // TEST_SUITE snc2
