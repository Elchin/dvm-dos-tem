// Stage B compute-cost benchmark: replay daily column states dumped from the legacy model
// (bench_patch.py) through the standalone enthalpy solver (include/EnthalpySolver.h).
// usage: bench_enthalpy [-r REPEATS] dump1.csv [dump2.csv ...]
// For each day: column = legacy start-of-day layers (T, water, solid capacity, texture, organic
// fraction, unfrozen/frozen conductivity), top BC = legacy surface temperature with the legacy
// Robin conductance (cn = 100 W m-2 K-1), zero bottom flux, dt = 1 day. Modes: NCZ sfcc 0,
// NCZ sfcc 1, Katzenelson sfcc 0. Legacy timing/solves per day are read from the dump (measured
// in-model around soilenv.updateDailyGroundT with std::chrono::steady_clock).
#include "../../include/EnthalpySolver.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace enthalpy;

struct Day { std::string file; int doy; double ts, legacy_ns; long legacy_solves; int snow; std::vector<Cell> c0, c1;
             std::vector<double> T0, Tleg; };

static std::vector<double> split(const std::string& s) {
  std::vector<double> v; std::stringstream ss(s); std::string tok;
  std::getline(ss, tok, ',');
  while (std::getline(ss, tok, ',')) v.push_back(atof(tok.c_str()));
  return v;
}

static void stats(std::vector<double> v, double& mean, double& p95, double& mx) {
  if (v.empty()) { mean = p95 = mx = 0; return; }
  double s = 0; for (double x : v) s += x; mean = s / v.size();
  std::sort(v.begin(), v.end()); p95 = v[(size_t)std::min(v.size() - 1, (size_t)ceil(0.95 * v.size()) - 1)]; mx = v.back();
}

int main(int argc, char** argv) {
  int rep = 3; std::vector<std::string> files;
  for (int i = 1; i < argc; ++i) { if (!strcmp(argv[i], "-r")) rep = atoi(argv[++i]); else files.push_back(argv[i]); }
  std::vector<Day> days;
  for (const std::string& f : files) {
    std::ifstream in(f.c_str()); std::string line;
    while (std::getline(in, line)) {
      if (line.compare(0, 2, "D,")) continue;
      std::vector<double> d = split(line);
      Day D; D.file = f; D.doy = (int)d[1]; D.ts = d[2]; D.legacy_ns = d[3]; D.legacy_solves = (long)d[4]; D.snow = 0;
      const int nl = (int)d[5];
      bool ok = true;
      for (int k = 0; k < nl; ++k) {
        if (!std::getline(in, line) || line.compare(0, 2, "L,")) { ok = false; break; }
        std::vector<double> L = split(line);   // kind dz T liq ice Cs W liqmin clay silt sand forg lu lf Tafter
        Cell c; c.kind = (int)L[0]; c.dz = L[1]; c.Cs = L[5]; c.W = std::max(0.0, L[6]); c.liqmin = L[7];
        c.clay = L[8]; c.silt = L[9]; c.sand = L[10]; c.forg = L[11]; c.lam_u = L[12]; c.lam_f = L[13]; c.eps = 0.01;
        if (c.kind == SNOW) { c.liqmin = 0.0; c.lam_f = c.lam_u; D.snow = 1; }
        if (c.kind == ROCK) { c.W = 0.0; c.lam_f = c.lam_u; }
        c.sfcc = SFCC_FREE; D.c0.push_back(c);
        if (c.kind == SOIL) c.sfcc = SFCC_SJ2016;
        D.c1.push_back(c); D.T0.push_back(L[2]); D.Tleg.push_back(L[14]);
      }
      if (ok && nl > 0) days.push_back(D);
    }
  }
  const size_t N = days.size();
  printf("days=%zu files=%zu repeats=%d\n", N, files.size(), rep);
  if (!N) return 1;
  const char* names[3] = {"NCZ_sfcc0", "NCZ_sfcc1", "Katz_sfcc0"};
  // phase-change classification (NCZ sfcc 0): any latent cell's liquid water changes by > 1e-6 kg m-2
  std::vector<int> pc(N, 0);
  {
    ColumnSolver s; Stats st; BC bc;
    for (size_t i = 0; i < N; ++i) {
      std::vector<double> T = days[i].T0; bc.Ts = days[i].ts; s.stepNCZ(days[i].c0, T, 86400.0, bc, 0, st);
      for (size_t k = 0; k < T.size(); ++k)
        if (hasLatent(days[i].c0[k]) && fabs(liq(days[i].c0[k], T[k]) - liq(days[i].c0[k], days[i].T0[k])) > 1e-6) { pc[i] = 1; break; }
    }
  }
  for (int cls = 0; cls < 3; ++cls) {   // 0 all, 1 phase-change days, 2 other days
    std::vector<double> lu, ls;
    for (size_t i = 0; i < N; ++i) if (cls == 0 || (cls == 1) == (pc[i] == 1)) { lu.push_back(days[i].legacy_ns * 1e-3); ls.push_back(days[i].legacy_solves); }
    double a, b, c, d, e, f; stats(lu, a, b, c); stats(ls, d, e, f);
    printf("legacy class=%s n=%zu us/day mean=%.2f p95=%.2f max=%.2f | tridiag solves/day mean=%.1f p95=%.0f max=%.0f\n",
           cls == 0 ? "all" : cls == 1 ? "phase" : "nophase", lu.size(), a, b, c, d, e, f);
  }
  for (int m = 0; m < 3; ++m) {
    std::vector<double> us(N, 1e30), sol(N), out(N), dtl(N, 0.0); double worst = 0; int nconv = 0;
    for (int r = 0; r < rep; ++r) {
      ColumnSolver s; Stats st; BC bc;
      for (size_t i = 0; i < N; ++i) {
        std::vector<double> T = days[i].T0; bc.Ts = days[i].ts;
        const std::vector<Cell>& c = (m == 1) ? days[i].c1 : days[i].c0;
        const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
        const bool ok = (m == 2) ? s.stepKatz(c, T, 86400.0, bc, 0, st) : s.stepNCZ(c, T, 86400.0, bc, 0, st);
        const double t = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
        us[i] = std::min(us[i], t);
        if (r == 0) { sol[i] = st.solves; out[i] = st.outer; worst = std::max(worst, fabs(st.closure) / 86400.0); nconv += ok;
                      for (size_t k = 0; k < T.size(); ++k) if (c[k].kind != ROCK) dtl[i] = std::max(dtl[i], fabs(T[k] - days[i].Tleg[k])); }
      }
    }
    double dm, dp, dx; stats(dtl, dm, dp, dx);
    printf("%s converged=%d/%zu max|closure|=%.3g W/m2 | per-day max|T-T_legacy| over snow/soil cells: mean=%.3f p95=%.3f max=%.2f C\n",
           names[m], nconv, N, worst, dm, dp, dx);
    for (int cls = 0; cls < 3; ++cls) {
      std::vector<double> u, so, ou;
      for (size_t i = 0; i < N; ++i) if (cls == 0 || (cls == 1) == (pc[i] == 1)) { u.push_back(us[i]); so.push_back(sol[i]); ou.push_back(out[i]); }
      double a, b, c, d, e, f, g, h, k; stats(u, a, b, c); stats(so, d, e, f); stats(ou, g, h, k);
      printf("  class=%s n=%zu us/day mean=%.2f p95=%.2f max=%.2f | solves/day mean=%.2f p95=%.0f max=%.0f | outer mean=%.2f p95=%.0f max=%.0f\n",
             cls == 0 ? "all" : cls == 1 ? "phase" : "nophase", u.size(), a, b, c, d, e, f, g, h, k);
    }
  }
  return 0;
}
