// Lightweight unit tests for RTS initiation helpers (no external framework).
// Build/run: make rts-unit-test
#include <cmath>
#include <cstdio>
#include <vector>
#include "Thermokarst.h"
#include "RTSMethane.h"
#include "RTSRegrow.h"
#include "Ground.h"
#include "Climate.h"
#include "TEMUtilityFunctions.h"
#include <algorithm>
// white-box access to Soil_Bgc internals for the burial conservation test
#define private public
#include "Soil_Bgc.h"
#include "Integrator.h"
#undef private
#include "CohortData.h"
#include "EnvData.h"
#include "BgcData.h"

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

static void test_m2d_mean_preserving() {
  struct V { double mean, amp, noise; bool nonneg, clip; } vars[] = {
    {-8.0, 22.0, 6.0, false, false},   // tair (Gydan-like)
    {0.4, 0.5, 0.2, true, true},       // vapo (kPa) with zeros
    {80.0, 110.0, 30.0, true, true},   // nirr W/m2, polar night zeros
    {40.0, 55.0, 15.0, true, true}};   // par
  int nyr = 20;
  for (const V& v : vars) {
    std::vector<float> ser = synth_series(nyr, v.mean, v.amp, v.noise, 7u, v.clip);
    int bad_mean = 0, bad_neg = 0, bad_len = 0;
    for (int y = 0; y < nyr; ++y) {
      std::vector<float> d = Climate::monthly2daily_mp(Climate::mp_window(ser, y), v.nonneg);
      if (d.size() != 365) { bad_len++; continue; }
      check_means(d, std::vector<float>(ser.begin() + 12 * y, ser.begin() + 12 * y + 12), bad_mean);
      if (v.nonneg) for (float x : d) if (x < 0) bad_neg++;
    }
    CHECK(bad_len == 0); CHECK(bad_mean == 0); CHECK(bad_neg == 0);
  }
  // climatology path (eq stage)
  std::vector<float> c12 = synth_series(1, 50.0, 70.0, 0.0, 3u, true);
  std::vector<float> d = Climate::monthly2daily_mp(Climate::mp_window_clim(c12), true);
  int bad = 0; check_means(d, c12, bad); CHECK(bad == 0);
  CHECK(*std::min_element(d.begin(), d.end()) >= 0.0f);
}

static void test_m2d_continuity() {
  // Dec 16..31 and Jan 1..15 lie on the same linear segment, so the jump
  // Dec31->Jan1 across years must equal the Dec30->Dec31 step.
  std::vector<float> ser = synth_series(15, -8.0, 22.0, 6.0, 11u, false);
  double worst = 0;
  for (int y = 0; y + 1 < 15; ++y) {
    std::vector<float> a = Climate::monthly2daily_mp(Climate::mp_window(ser, y), false);
    std::vector<float> b = Climate::monthly2daily_mp(Climate::mp_window(ser, y + 1), false);
    worst = std::max(worst, (double)std::fabs((b[0] - a[364]) - (a[364] - a[363])));
  }
  CHECK(worst < 1e-5);
  // periodic climatology: Dec31 -> Jan1 of same curve continuous
  std::vector<float> c12 = synth_series(1, -8.0, 22.0, 0.0, 3u, false);
  std::vector<float> d = Climate::monthly2daily_mp(Climate::mp_window_clim(c12), false);
  CHECK(std::fabs((d[0] - d[364]) - (d[364] - d[363])) < 1e-5);
  // within-year: no step larger than the max monthly-mean jump
  double maxstep = 0; for (int i = 1; i < 365; ++i) maxstep = std::max(maxstep, (double)std::fabs(d[i] - d[i-1]));
  CHECK(maxstep < 3.0);
}

static void test_m2d_constant_and_tdd() {
  // constant input -> constant output exactly; and the convex July peak is
  // no longer flattened (legacy loses TDD)
  std::vector<float> k(36, 4.25f);
  std::vector<float> d = Climate::monthly2daily_mp(k, false);
  double mx = 0; for (float x : d) mx = std::max(mx, (double)std::fabs(x - 4.25f));
  CHECK(mx < 1e-6);
  std::vector<float> m12 = {-25, -26, -22, -15, -5, 5, 12, 9, 3, -7, -18, -23};
  Climate clm;
  std::vector<float> leg = clm.monthly2daily(std::vector<float>{-23, -25, -26, -22, -15, -5, 5, 12, 9, 3, -7, -18, -23, -25});
  std::vector<float> mp = Climate::monthly2daily_mp(Climate::mp_window_clim(m12), false);
  double tdd_true = 0, tdd_leg = 0, tdd_mp = 0;
  for (int m = 0; m < 12; ++m) tdd_true += std::max(0.0f, m12[m]) * DIM12[m];
  for (int i = 0; i < 365; ++i) { tdd_leg += std::max(0.0f, leg[i]); tdd_mp += std::max(0.0f, mp[i]); }
  CHECK(leg.size() == 366);                 // legacy container length unchanged
  CHECK(tdd_leg < tdd_mp);                  // legacy under-estimates thaw
  CHECK(std::fabs(tdd_mp - tdd_true) < 0.06 * tdd_true);
}

static void test_m2d_legacy_bit_identical() {
  // The legacy routine must match a verbatim re-implementation of the old code.
  static const float arr[] = { -15.5, 15.5, 45.0, 74.5, 105.0, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 380.5 };
  std::vector<float> v14 = {-23.1f, -25.2f, -26.3f, -22.4f, -15.5f, -5.6f, 5.7f, 12.8f, 9.9f, 3.1f, -7.2f, -18.3f, -23.4f, -25.5f};
  std::vector<float> ref;
  for (int i = 1; i < 14; ++i) {
    std::vector<float> psd = temutil::resample(std::make_pair(arr[i-1], v14[i-1]), std::make_pair(arr[i], v14[i]),
                                               int(floor(arr[i-1])), int(floor(arr[i])), 1);
    ref.insert(ref.end(), psd.begin(), psd.end());
  }
  ref = std::vector<float>(ref.begin() + 16, ref.end() - 14);
  Climate clm;
  std::vector<float> got = clm.monthly2daily(v14);
  CHECK(got.size() == ref.size());
  bool same = got.size() == ref.size();
  for (size_t i = 0; same && i < got.size(); ++i) same = (got[i] == ref[i]);
  CHECK(same);
}

// deep_carbon_v2 burial (Soil_Bgc::deltastate part II) must only move C and N
// between layers. Column: fibric | humic | mineral 2.6 m | mineral 2.6 m, so the
// second mineral layer exhausts xtopmlthick (5.2 m) exactly.
static void test_burial_conserves_c_and_n() {
  CohortData cd; EnvData ed; BgcData bd; Soil_Bgc sb;
  sb.setCohortData(&cd); sb.setEnvData(&ed); sb.setBgcData(&bd);
  const int n = 4; cd.m_soil.numsl = n;
  int type[n] = {1, 2, 3, 3}; double dz[n] = {0.05, 0.10, 2.6, 2.6};
  double z = 0.0;
  for (int i = 0; i < n; ++i) { cd.m_soil.type[i] = type[i]; cd.m_soil.dz[i] = dz[i]; cd.m_soil.z[i] = z; z += dz[i]; }
  cd.m_soil.deepthick = 0.10;
  sb.bgcpar.som2co2 = 0.5; sb.bgcpar.fsoma = 0.3; sb.bgcpar.fsompr = 0.5; sb.bgcpar.fsomcr = 0.2;
  sb.set_nfeed(1); sb.set_avlnflg(0); sb.set_baseline(0);
  sb.d2wdebrisc = 0.; sb.d2wdebrisn = 0.; sb.totdzliq = 0.; sb.del_soi2l.avlnlost = 0.;
  sb.tmp_sois.wdebrisc = 0.; sb.tmp_sois.wdebrisn = 0.; sb.del_soi2a.rhwdeb = 0.;
  bd.m_a2soi.avlninput = 0.; bd.m_a2soi.orgninput = 0.;
  double rhr[n] = {1.0, 0.2, 0.0, 0.0}, rha[n] = {0.0, 1.0, 0.0, 0.0}, lc[n] = {5.0, 0.5, 0.0, 0.0};
  double expC = 0.0;
  for (int i = 0; i < n; ++i) {
    sb.del_soi2a.rhrawc[i] = rhr[i]; sb.del_soi2a.rhsoma[i] = rha[i];
    sb.del_soi2a.rhsompr[i] = 0.1; sb.del_soi2a.rhsomcr[i] = 0.05;
    sb.ltrflc[i] = lc[i]; sb.ltrfln[i] = 0.0; sb.rtnextract[i] = 0.0;
    sb.del_soi2soi.netnmin[i] = 0.0; bd.m_v2soi.rtlfalfrac[i] = 0.0;
    sb.tmp_sois.rawc[i] = 100.0 * (i + 1); sb.tmp_sois.soma[i] = 50.0; sb.tmp_sois.sompr[i] = 80.0;
    sb.tmp_sois.somcr[i] = 120.0; sb.tmp_sois.orgn[i] = 10.0 * (i + 1); sb.tmp_sois.avln[i] = 0.1;
    double rhsum = rhr[i] + rha[i] + 0.15;
    expC += lc[i] - rhsum * (1.0 + 0.5) + rhsum * 0.5 * (0.3 + 0.5 + 0.2);   // part (I) only
  }
  sb.deltastate();
  double dC = 0.0, dN = 0.0;
  for (int i = 0; i < n; ++i) {
    dC += sb.del_sois.rawc[i] + sb.del_sois.soma[i] + sb.del_sois.sompr[i] + sb.del_sois.somcr[i];
    dN += sb.del_sois.orgn[i];
  }
  CHECK_NEAR(dC, expC, 1e-10);      // burial moves C only (no loss at the 5.2 m boundary)
  CHECK_NEAR(dN, 0.0, 1e-10);       // burial moves N only (no overwrite in humic layers)
  double c3 = sb.del_sois.rawc[3] + sb.del_sois.soma[3] + sb.del_sois.sompr[3] + sb.del_sois.somcr[3];
  CHECK(c3 > 0.0);                  // layer exhausting xtopmlthick receives its share
  CHECK(sb.del_sois.orgn[1] > -sb.tmp_sois.orgn[1]);  // humic N not over-drawn
}

// Negative layer ORGN (AVLN-deficit borrow in Integrator::y2cstate_soi) is covered
// from other layers, root zone first; column N is conserved exactly.
static void test_cover_negative_orgn() {
  double o[5] = {2.0, -0.5, 1.0, 10.0, -0.25}; double rf[5] = {0.5, 0.3, 0.2, 0.0, 0.0};
  double before = 0; for (int i = 0; i < 5; i++) before += o[i];
  double moved = Integrator::cover_negative_orgn(o, rf, 5);
  double after = 0; bool nonneg = true; for (int i = 0; i < 5; i++) { after += o[i]; nonneg = nonneg && o[i] >= 0; }
  CHECK_NEAR(moved, 0.75, 1e-12); CHECK_NEAR(after, before, 1e-12); CHECK(nonneg);
  CHECK_NEAR(o[3], 10.0, 1e-12);                 // root-zone layers covered it
  double o2[3] = {0.1, -1.0, 5.0}; double rf2[3] = {1.0, 1.0, 0.0};
  Integrator::cover_negative_orgn(o2, rf2, 3);   // root zone short -> rest from below
  CHECK_NEAR(o2[0] + o2[1] + o2[2], 4.1, 1e-12); CHECK(o2[0] >= 0 && o2[1] >= 0 && o2[2] >= 0);
  double o3[2] = {0.2, -1.0}; double rf3[2] = {1.0, 1.0};
  Integrator::cover_negative_orgn(o3, rf3, 2);   // column short: no N created
  CHECK_NEAR(o3[0] + o3[1], -0.8, 1e-12);
  double o4[3] = {1.0, 2.0, 3.0}; double rf4[3] = {1, 1, 1};
  CHECK(Integrator::cover_negative_orgn(o4, rf4, 3) == 0.0 && o4[1] == 2.0);  // no-op
}

// ---- RTS CH4 scheme (a): reduced-form partition of RH ----
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
    CHECK(f >= p.oxid_min - 1e-15 && f <= p.oxid_max + 1e-15 && f >= 0.0 && f <= 1.0);
    CHECK(f >= prev - 1e-15); prev = f;
  }
  RTSMethaneParams q = p; q.oxid_min = -0.5; q.oxid_max = 1.7;     // bad parameters are clamped
  CHECK(RTSMethane::oxidFraction(q, 0.0) >= 0.0 && RTSMethane::oxidFraction(q, 5.0) <= 1.0);
  // ponding: 0 before initiation, pond_sat0 at initiation, 1/e after pond_tau years
  CHECK(RTSMethane::pondSaturation(p, -1) == 0.0);
  CHECK_NEAR(RTSMethane::pondSaturation(p, 0), p.pond_sat0, 1e-15);
  CHECK_NEAR(RTSMethane::pondSaturation(p, (int)(12 * p.pond_tau)), p.pond_sat0 / std::exp(1.0), 1e-12);
  // methanogen recovery: 1 before initiation, mg0 at initiation, 1-(1-mg0)/e after mg_tau years, -> 1
  CHECK(RTSMethane::methanogenFactor(p, -1) == 1.0);
  CHECK_NEAR(RTSMethane::methanogenFactor(p, 0), p.mg0, 1e-15);
  CHECK_NEAR(RTSMethane::methanogenFactor(p, (int)(12 * p.mg_tau)), 1.0 - (1.0 - p.mg0) / std::exp(1.0), 1e-12);
  CHECK(RTSMethane::methanogenFactor(p, 12000) > 0.999999);
  { RTSMethaneParams q = p; q.mg0 = 1.0; CHECK(RTSMethane::methanogenFactor(q, 5) == 1.0); }
  // production <= RH (also with an extreme scale factor) and emis + oxid = prod
  std::vector<RTSMethaneLayer> col(6);
  double rh = 0.0;
  for (int i = 0; i < 6; i++) { col[i].ztop = 0.1 * i; col[i].dz = 0.1; col[i].tem = 25.0 - 4 * i; col[i].rh = 1.0 + i; rh += col[i].rh; }
  for (double frac : {0.0, 0.01, 0.5, 5.0, 100.0}) {
    for (double wtd : {0.0, 0.25, 2.0}) {
      for (double pond : {0.0, 0.5, 1.0}) {
        RTSMethaneParams r = p; r.ch4_frac = frac; r.rel_q10 = 3.0;
        RTSMethaneFluxes f = RTSMethane::compute(col, r, wtd, pond, 2.0);
        CHECK(f.prod >= 0.0 && f.prod <= rh + 1e-12);
        CHECK(f.oxid >= 0.0 && f.oxid <= f.prod + 1e-12 && f.emis >= 0.0);
        CHECK_NEAR(f.emis + f.oxid, f.prod, 1e-12);
        CHECK(f.plant <= f.emis + 1e-12);
      }
    }
  }
  // fully saturated, frac=100: prod == RH of thawed layers exactly (cap at 1)
  { RTSMethaneParams r = p; r.ch4_frac = 100.0;
    double rh_thawed = 0.0; for (auto & l : col) if (l.tem > 0) rh_thawed += l.rh;
    CHECK_NEAR(RTSMethane::compute(col, r, 0.0, 0.0, 0.0).prod, rh_thawed, 1e-12); }
  // frozen or no RH -> no CH4; pond=0 + deep WT -> no CH4 (A_RTS=0 / no anoxia baseline)
  { std::vector<RTSMethaneLayer> fz(col); for (auto & l : fz) l.tem = -1.0;
    CHECK(RTSMethane::compute(fz, p, 0.0, 1.0, 1.0).prod == 0.0);
    CHECK(RTSMethane::compute(col, p, 10.0, 0.0, 1.0).prod == 0.0); }
  // ponding raises production with a deep water table
  CHECK(RTSMethane::compute(col, p, 0.3, 0.5, 1.0).prod > RTSMethane::compute(col, p, 0.3, 0.0, 1.0).prod);
}

static void test_rts_floor() {
  RTSFloorParams p;
  // v0.2 defaults: all switches on, calibrated values; all off -> any() false (bitwise v0.1 path)
  CHECK(p.rh_cap && p.exp_on && p.anox && p.any());
  CHECK(p.cap0 == 0.045 && p.cap_inf == 0.002 && p.cap_tau == 1.5 && p.lat0 == 0.02 && p.lat_tau == 5.0 && p.anox_f == 0.3);
  { RTSFloorParams q = p; q.rh_cap = q.exp_on = q.anox = false; CHECK(q.any() == q.rmin_on); q.rmin_on = false; CHECK(!q.any()); }
  { RTSMethaneParams mp; CHECK(mp.mg0 == 1.0); CHECK(mp.ch4_frac == 0.135); }
  CHECK(p.rmin_on && p.rmin0 == 0.3 && p.rmin_tau == 10.0 && p.rmin_years == 5);
  { RTSFloorParams q = p; q.rh_cap = q.exp_on = q.anox = q.rmin_on = false; CHECK(!q.any()); }
  // rmin factor: 1 before initiation, rmin0 at initiation, 1-(1-rmin0)/e after rmin_tau years, -> 1
  CHECK(RTSFloor::rminFactor(p, -1) == 1.0);
  CHECK_NEAR(RTSFloor::rminFactor(p, 0), p.rmin0, 1e-15);
  CHECK_NEAR(RTSFloor::rminFactor(p, (int)(12 * p.rmin_tau)), 1.0 - (1.0 - p.rmin0) / std::exp(1.0), 1e-12);
  CHECK(RTSFloor::rminFactor(p, 120000) > 0.999999);
  { RTSFloorParams q = p; q.rmin0 = -2.0; CHECK(RTSFloor::rminFactor(q, 0) == 0.0); q.rmin0 = 3.0; CHECK(RTSFloor::rminFactor(q, 0) == 1.0); }
  // minimum respiration: raises RH to the target from the pools, conserves C, never above rmin_maxfrac of C
  { const int n = 3; double w[n] = {1.0, 0.5, 0.0};
    double cr[n] = {10, 0, 5}, ca[n] = {100, 200, 50}, cp[n] = {1000, 500, 0}, cc[n] = {5000, 3000, 100};
    double rr[n] = {0.01, 0, 0}, ra[n] = {0.1, 0.2, 0}, rp[n] = {0.05, 0.02, 0}, rc[n] = {0.01, 0.01, 0};
    double r0 = 0; for (int i = 0; i < n; i++) r0 += rr[i] + ra[i] + rp[i] + rc[i];
    const double add = RTSFloor::applyMinimum(10.0, n, w, cr, ca, cp, cc, rr, ra, rp, rc);
    double r1 = 0; for (int i = 0; i < n; i++) r1 += rr[i] + ra[i] + rp[i] + rc[i];
    CHECK_NEAR(r1, 10.0, 1e-9); CHECK_NEAR(r1 - r0, add, 1e-12);
    CHECK(rr[2] == 0 && ra[2] == 0 && rc[2] == 0);                 // zero weight (e.g. no activity): untouched
    for (int i = 0; i < n; i++) CHECK(rr[i] + ra[i] + rp[i] + rc[i] <= RTSFloor::rmin_maxfrac * (cr[i] + ca[i] + cp[i] + cc[i]) + 1e-12);
    // above target: unchanged; huge target: limited by rmin_maxfrac (no C creation beyond pools)
    double before = rr[0]; CHECK(RTSFloor::applyMinimum(1.0, n, w, cr, ca, cp, cc, rr, ra, rp, rc) == 0.0); CHECK(rr[0] == before);
    RTSFloor::applyMinimum(1e9, n, w, cr, ca, cp, cc, rr, ra, rp, rc);
    for (int i = 0; i < n; i++) CHECK(rr[i] + ra[i] + rp[i] + rc[i] <= RTSFloor::rmin_maxfrac * (cr[i] + ca[i] + cp[i] + cc[i]) + 1e-9);
    // empty pools: nothing added
    double z[n] = {0, 0, 0}, a2[n] = {0, 0, 0}, b2[n] = {0, 0, 0}, c2[n] = {0, 0, 0}, d2[n] = {0, 0, 0};
    CHECK(RTSFloor::applyMinimum(5.0, n, w, z, z, z, z, a2, b2, c2, d2) == 0.0); }
  // cap rate: cap0 at initiation, cap_inf + (cap0-cap_inf)/e after cap_tau years
  CHECK_NEAR(RTSFloor::capRate(p, 0), p.cap0, 1e-15);
  CHECK_NEAR(RTSFloor::capRate(p, (int)(12 * p.cap_tau + 0.5)) , p.cap_inf + (p.cap0 - p.cap_inf) * std::exp(-std::floor(12 * p.cap_tau + 0.5) / (12 * p.cap_tau)), 1e-12);
  CHECK(RTSFloor::capRate(p, 100000) - p.cap_inf < 1e-12);
  // layer cap: no change before initiation or below the cap, limited above it
  CHECK(RTSFloor::capLayerRH(p, -1, 50.0, 100.0) == 50.0);
  CHECK(RTSFloor::capLayerRH(p, 0, 1e-6, 1000.0) == 1e-6);
  CHECK_NEAR(RTSFloor::capLayerRH(p, 0, 1e6, 1200.0), p.cap0 / 12.0 * 1200.0, 1e-12);
  CHECK(RTSFloor::capLayerRH(p, 3, 0.0, 1000.0) == 0.0);
  // export fraction: lat0/12 at initiation, e-folding lat_tau, 0 before
  CHECK(RTSFloor::exportFraction(p, -1) == 0.0);
  CHECK_NEAR(RTSFloor::exportFraction(p, 0), p.lat0 / 12.0, 1e-15);
  CHECK_NEAR(RTSFloor::exportFraction(p, (int)(12 * p.lat_tau)), p.lat0 / 12.0 / std::exp(1.0), 1e-12);
  // anoxia: 1 for unsaturated, anox_f for saturated, linear between, clipped
  CHECK(RTSFloor::anoxFactor(p, 0.0) == 1.0);
  CHECK_NEAR(RTSFloor::anoxFactor(p, 1.0), p.anox_f, 1e-15);
  CHECK_NEAR(RTSFloor::anoxFactor(p, 0.5), 0.5 * (1.0 + p.anox_f), 1e-15);
  CHECK(RTSFloor::anoxFactor(p, 2.0) == RTSFloor::anoxFactor(p, 1.0));
  // export of a pool set conserves mass
  { double c[4] = {10.0, 200.0, 3000.0, 40000.0}, out = 0.0, tot0 = 0.0, tot1 = 0.0;
    const double f = RTSFloor::exportFraction(p, 7);
    for (double & x : c) { tot0 += x; const double d = f * x; x -= d; out += d; tot1 += x; }
    CHECK_NEAR(tot0, tot1 + out, 1e-9); }
}

static void test_rts_regrow() {
  RTSRegrowParams p;
  CHECK(!p.on);                                            // default OFF
  // Xia et al. 2026: tau = 1.35 GPP^-1.68 (GPP kgC/m2/yr): 0.3 -> ~10.2 yr, s ~1.46
  CHECK_NEAR(RTSRegrow::tau(p, 0.3), 1.35 * std::pow(0.3, -1.68), 1e-12);
  CHECK(RTSRegrow::tau(p, 0.3) > 10.0 && RTSRegrow::tau(p, 0.3) < 10.5);
  CHECK_NEAR(RTSRegrow::tau(p, 0.3) / p.tau_ref, 10.2 / 7.0, 0.02);  // site scale s (diagnostic in v2)
  CHECK(RTSRegrow::tau(p, 0.0) == p.tau_max);              // no GPP -> slowest
  CHECK(RTSRegrow::tau(p, 100.0) == p.tau_min);            // clamped
  // PFT classes
  CHECK(RTSRegrow::pftType(0, 0, 3) == RG_LICHEN);
  CHECK(RTSRegrow::pftType(0, 0, 2) == RG_MOSS && RTSRegrow::pftType(0, 0, 1) == RG_MOSS);
  CHECK(RTSRegrow::pftType(1, 1, 0) == RG_SHRUB && RTSRegrow::pftType(1, 0, 0) == RG_EVERGREEN);
  CHECK(RTSRegrow::pftType(0, 0, 0) == RG_GRAM);
  // time axis: stabilization 12 months after initiation by default
  CHECK(RTSRegrow::stabMonth(p) == 12);
  CHECK_NEAR(RTSRegrow::tYears(p, 36), 2.0, 1e-12);
  // v2: no logistic growth ceiling exists (targetFrac/limitMult removed by design; flagged in commit).
  // shrub boost: B0 at t=0, decays with Tb, only for deciduous shrubs, none while active
  CHECK_NEAR(RTSRegrow::boost(p, RG_SHRUB, 0.0), p.boost0, 1e-12);
  CHECK_NEAR(RTSRegrow::boost(p, RG_SHRUB, p.boost_T), 1.0 + (p.boost0 - 1.0) * std::exp(-1.0), 1e-12);
  CHECK(RTSRegrow::boost(p, RG_GRAM, 0.0) == 1.0 && RTSRegrow::boost(p, RG_SHRUB, -0.1) == 1.0);
  // OL rebuild and N supply rates (Pizano et al. 2014)
  CHECK(RTSRegrow::olRate(p, -0.1) == 0.0 && RTSRegrow::olRate(p, 10.0) == 32.0 && RTSRegrow::olRate(p, 61.0) == 6.5);
  CHECK(RTSRegrow::nRate(p, -0.1) == 0.0 && RTSRegrow::nRate(p, 59.9) == 2.2 && RTSRegrow::nRate(p, 60.1) == 0.4);
  // avln part + nominal OL N = total N supply (annual), never negative
  CHECK_NEAR(12.0 * RTSRegrow::avlnSupplyMonth(p, 10.0) + 32.0 / 30.0, 2.2, 1e-12);
  CHECK_NEAR(12.0 * RTSRegrow::avlnSupplyMonth(p, 70.0) + 6.5 / 30.0, 0.4, 1e-12);
  RTSRegrowParams q; q.ol_cn = 5.0; CHECK(RTSRegrow::avlnSupplyMonth(q, 10.0) == 0.0);
  // v2 GPP-recovery taper: full supply when GPP collapsed, ~zero when recovered, monotone in between
  CHECK(p.n_taper && p.n_taper_pow == 4.0);
  CHECK(RTSRegrow::taperFactor(p, 0.0) == 1.0 && RTSRegrow::taperFactor(p, -1.0) == 1.0);
  CHECK_NEAR(RTSRegrow::taperFactor(p, 0.5), 1.0 - 0.0625, 1e-12);   // plateau: near-full at half recovery
  CHECK_NEAR(RTSRegrow::taperFactor(p, 0.9), 1.0 - 0.6561, 1e-12);
  CHECK(RTSRegrow::taperFactor(p, 1.0) == 0.0 && RTSRegrow::taperFactor(p, 2.0) == 0.0);
  { double prev = 2.0; bool mono = true;
    for (int i = 0; i <= 100; i++) { const double v = RTSRegrow::taperFactor(p, i / 100.0); mono &= v <= prev; prev = v; }
    CHECK(mono); }
  { RTSRegrowParams q2; q2.n_taper_pow = 2.0; CHECK_NEAR(RTSRegrow::taperFactor(q2, 0.5), 0.75, 1e-12); }
  { RTSRegrowParams q3; q3.n_taper = false; CHECK(RTSRegrow::taperFactor(q3, 0.99) == 1.0); }  // untapered v1 supply
  // tapered avln supply scales linearly with the factor and is never negative
  CHECK_NEAR(RTSRegrow::avlnSupplyMonth(p, 10.0, 0.5), 0.5 * RTSRegrow::avlnSupplyMonth(p, 10.0, 1.0), 1e-15);
  CHECK(RTSRegrow::avlnSupplyMonth(p, 10.0, 0.0) == 0.0 && RTSRegrow::avlnSupplyMonth(p, 10.0, -1.0) == 0.0);
}

int main() {
  test_rts_regrow();
  test_rts_floor();
  test_rts_ch4();
  test_excess_ice_init();
  test_excess_melt_subsidence();
  test_removal_thickness();
  test_partition_mass_balance();
  test_thaw_column();
  test_m2d_mean_preserving();
  test_m2d_continuity();
  test_m2d_constant_and_tdd();
  test_m2d_legacy_bit_identical();
  test_burial_conserves_c_and_n();
  test_cover_negative_orgn();
  std::printf("RTS unit tests: %d passed, %d failed\n", g_pass, g_fail);
  return g_fail == 0 ? 0 : 1;
}
