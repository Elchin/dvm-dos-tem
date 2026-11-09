#!/usr/bin/env python3
"""Stage B compute-cost evaluation: instrument a THROW-AWAY COPY of the model (never the branch)
to time the legacy thermal step and dump the daily column state.
usage: bench_patch.py SRC_TREE   (edits SRC_TREE/src/Cohort.cpp and SRC_TREE/src/CrankNicholson.cpp in place)
At run time: DVM_THERMAL_BENCH=/path/prefix  -> writes prefix.<pid>.csv (first DVM_THERMAL_BENCH_DAYS days,
default 3650) and prints 'THERMAL_BENCH ...' totals to stderr at exit. Unset: only one getenv per day.
Timed region: the soilenv.updateDailyGroundT(tdrv, daylength) call only (the whole legacy thermal step:
surface flux, Stefan fronts, Crank-Nicolson iterations, layer state update, front/ALD retrieval);
the dump is built before the timer starts and written after it stops."""
import sys, re
root = sys.argv[1]
cn = f"{root}/src/CrankNicholson.cpp"
s = open(cn).read()
assert "g_cn_solves" not in s
s = s.replace('#include "../include/CrankNicholson.h"\n',
              '#include "../include/CrankNicholson.h"\nlong g_cn_solves = 0; // BENCH\n', 1)
s = s.replace("  double denm_inv;\n\n  for (int il = endind-1 ; il>=startind+1; il--) {",
              "  double denm_inv;\n  g_cn_solves++; // BENCH\n\n  for (int il = endind-1 ; il>=startind+1; il--) {", 1)
assert s.count("g_cn_solves") == 2
open(cn, "w").write(s)

co = f"{root}/src/Cohort.cpp"
s = open(co).read()
old = "    soilenv.updateDailyGroundT(tdrv, daylength);\n"
assert s.count(old) == 1
new = r'''    { // BENCH (throw-away instrumentation, see tools/enthalpy_bench/bench_patch.py)
      static const char* bpre = getenv("DVM_THERMAL_BENCH");
      static BenchTot btot;
      static std::ofstream* bout = 0; static long bday = 0; static long bmax = 3650;
      std::ostringstream pre;
      const bool dump = bpre && bday < bmax;
      if (bpre && !bout) {
        if (getenv("DVM_THERMAL_BENCH_DAYS")) bmax = atol(getenv("DVM_THERMAL_BENCH_DAYS"));
        std::ostringstream fn; fn << bpre << "." << getpid() << ".csv"; bout = new std::ofstream(fn.str().c_str());
        bout->precision(17);
      }
      int nl = 0;
      if (dump) {
        pre.precision(17);
        for (Layer* l = ground.toplayer; l != NULL; l = l->nextl) {
          int kind = l->isSnow ? 1 : (l->isRock ? 2 : 0);
          double Cs = l->isSoil ? l->vhcsolid * (1.0 - l->poro) * l->dz : (l->isRock ? l->vhcsolid * l->dz : 0.0);
          double cl = 0, si = 0, sa = 0, forg = 1.0;
          MineralLayer* ml = dynamic_cast<MineralLayer*>(l);
          if (ml) { cl = ml->pctclay; si = ml->pctsilt; sa = ml->pctsand;
                    forg = Soil_Bgc::sjOrgFraction(l->rawc + l->soma + l->sompr + l->somcr, l->dz); }
          pre << "L," << kind << "," << l->dz << "," << l->tem << "," << l->liq << "," << l->ice << "," << Cs << ","
              << (kind == 2 ? 0.0 : l->liq + l->ice) << "," << (l->isSoil ? l->minliq : 0.0) << "," << cl << "," << si << ","
              << sa << "," << forg << "," << l->getUnfThermCond() << "," << l->getFrzThermCond() << ",";
          nl++;
        }
      }
      const long s0 = g_cn_solves;
      const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
      soilenv.updateDailyGroundT(tdrv, daylength);
      const double ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
      const long ds = g_cn_solves - s0;
      if (bpre) { btot.ns += ns; btot.days++; btot.solves += ds; }
      if (dump) {
        const double ts = ground.seb_on ? ground.seb_ts : tdrv * edall->d_soid.nfactor;
        *bout << "D," << bday << "," << doy << "," << ts << "," << ns << "," << ds << "," << nl << "," << (edall->d_snw2soi.melt > 0. ? 1 : 0) << "\n";
        std::string p = pre.str(); size_t a = 0; Layer* l = ground.toplayer;
        for (int k = 0; k < nl; ++k) { size_t e = p.find("L,", a + 1); if (e == std::string::npos) e = p.size();
          *bout << p.substr(a, e - a) << (l ? l->tem : -999.0) << "\n"; a = e; if (l) l = l->nextl; }
        bday++;
      }
    }
'''
s = s.replace(old, new, 1)
s = s.replace('#include <limits>\n', '''#include <limits>
#include <chrono>
#include <sstream>
#include <unistd.h>
#include "../include/MineralLayer.h"
extern long g_cn_solves; // BENCH
struct BenchTot { double ns; long days, solves; BenchTot(): ns(0), days(0), solves(0) {}
  ~BenchTot() { if (days) fprintf(stderr, "THERMAL_BENCH days=%ld thermal_s=%.3f us_per_day=%.3f cn_solves=%ld\\n",
                days, ns * 1e-9, ns * 1e-3 / days, solves); } }; // BENCH
''', 1)
open(co, "w").write(s)
print("patched", root)
