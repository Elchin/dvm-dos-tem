#ifndef BNF_H_
#define BNF_H_
// v5q: biological N fixation (BNF), model_settings.bnf_scheme (default 0 = off).
// DVM-DOS-TEM has no BNF; with nfeed = 1 the only N input to the soil
// (m_a2soi.avlninput, output NINPUT) is 0. With bnf_scheme > 0 an annual BNF
// rate is computed from the PREVIOUS model year and given to the soil
// inorganic N pool (AVLN of the top soil layer) month by month:
//   1 PRO   B = c * (1 - exp(d * NPP))  (gN/m2/yr; NPP in gC/m2/yr, < 0 -> 0)
//           c = 1.8 gN/m2/yr, d = -0.003 m2 yr/gC (Meyerholt et al. 2016,
//           Eq. A2 / Table A1; from Cleveland et al. 1999 as used in CLM-CN,
//           Thornton et al. 2007)
//   2 AET   B = max(0, a * ET + b)  (ET in mm/yr: transpiration + canopy
//           evaporation, i.e. actual ET without soil (and moss) evaporation
//           and sublimation) a = 0.00234 gN/mm/m2, b = -0.0172 gN/m2/yr
//           (Meyerholt et al. 2016, Eq. A1 / Table A1; Cleveland et al. 1999)
//   3 CONST B = bnf_const (gN/m2/yr; per CMT through the case configuration)
// plus an optional constant free-living (asymbiotic, e.g. moss-associated
// cyanobacteria) term bnf_asym (gN/m2/yr, default 0; the soil temperature,
// light and moisture modulation of Meyerholt et al. 2016 Eqs. A20-A23 is NOT
// implemented) and an optional taper on soil inorganic N bnf_avln_cap
// (default 0 = no taper; else factor max(0, 1 - AVLN/cap), AVLN = column
// total at the start of the month; Meyerholt's FOR scheme only states that
// BNF converges to 0 above 2 gN/m2 inorganic N, the linear form is ours).
// Month m receives B * w_m, w_m = max(NPP_m,0)/sum(max(NPP_k,0)) of the
// previous year (uniform 1/12 if that sum is 0). The first model year of a
// process (no previous year in memory; not stored in the restart file) gets 0.
// Phase 2, model_settings.bnf_to_plant (0/1, default 0; 1 needs bnf_scheme > 0):
// 1 = the scheme part of B (the symbiotic-type rate: PRO, AET or CONST) goes to
// the plant labile N pool (LABN) instead of the soil, as symbiotic BNF in
// Meyerholt et al. 2016 (the plant receives the fixed N), shared among the PFTs
// with vegcov > 0 by their previous-year positive annual NPP (equal shares if
// all are <= 0); bnf_asym (free-living) stays in the soil. Output BNF_VEG.
// The NPP-proportional split is ours (no published source; consistent with the
// NPP-driven PRO rate).
// Pure helpers (no model state) so they can be unit tested
// (tests/unit/test_bnf_v5q.cpp).
// Refs: Meyerholt, Zaehle & Smith 2016, doi:10.5194/bg-13-1491-2016;
//   Cleveland et al. 1999, doi:10.1029/1999GB900014; Thornton et al. 2007,
//   doi:10.1029/2006GB002868; Davies-Barnard & Friedlingstein 2020,
//   doi:10.1029/2019GB006387; moss-associated fixation in the Arctic:
//   DeLuca et al. 2002, doi:10.1038/nature01051; Stewart et al. 2011,
//   doi:10.1657/1938-4246-43.2.267; Rousk et al. 2013, doi:10.3389/fmicb.2013.00150.
//   See docs/daily-mode-REFERENCES.md.
#include <cmath>
#include <string>
#include <json/json.h>

namespace bnf {

struct Params {
  int scheme = 0;          // 0 off, 1 PRO, 2 AET, 3 CONST
  double pro_c = 1.8;      // gN/m2/yr
  double pro_d = -0.003;   // m2 yr/gC
  double aet_a = 0.00234;  // gN/mm/m2
  double aet_b = -0.0172;  // gN/m2/yr
  double cnst  = 0.0;      // gN/m2/yr (scheme 3)
  double asym  = 0.0;      // gN/m2/yr, free-living term (scheme > 0 only)
  double avln_cap = 0.0;   // gN/m2, 0 = no taper
  int to_plant = 0;        // Phase 2: 1 = scheme part to plant labile N
};

const int MAX_PFT = 10;    // = NUM_PFT (cohortconst.h); checked in Soil_Bgc.cpp

// An integer-valued number in [lo, hi] (booleans, strings, fractions rejected).
inline bool valid_int(const Json::Value& v, int lo, int hi) {
  if (v.isBool() || !v.isNumeric()) { return false; }
  double d = v.asDouble();
  return std::isfinite(d) && d == std::floor(d) && d >= lo && d <= hi;
}
// A finite number (booleans and strings rejected).
inline bool valid_num(const Json::Value& v) {
  return !v.isBool() && v.isNumeric() && std::isfinite(v.asDouble());
}

// Parse model_settings; returns "" or the FATAL message.
inline std::string parse(const Json::Value& ms, int n_rootless_fix, Params& p) {
  p = Params();
  Json::Value s = ms.get("bnf_scheme", 0);
  if (!valid_int(s, 0, 3)) {
    return "FATAL: model_settings.bnf_scheme must be 0 (off), 1 (PRO), 2 (AET) or 3 (CONST)";
  }
  p.scheme = s.asInt();
  const char* keys[] = {"bnf_pro_c", "bnf_pro_d", "bnf_aet_a", "bnf_aet_b",
                        "bnf_const", "bnf_asym", "bnf_avln_cap"};
  double* dst[] = {&p.pro_c, &p.pro_d, &p.aet_a, &p.aet_b, &p.cnst, &p.asym, &p.avln_cap};
  for (int k = 0; k < 7; ++k) {
    if (ms.isMember(keys[k])) {
      if (!valid_num(ms[keys[k]])) {
        return std::string("FATAL: model_settings.") + keys[k] + " must be a finite number";
      }
      *dst[k] = ms[keys[k]].asDouble();
    }
  }
  Json::Value tp = ms.get("bnf_to_plant", 0);
  if (!valid_int(tp, 0, 1)) {
    return "FATAL: model_settings.bnf_to_plant must be 0 (soil) or 1 (plant labile N)";
  }
  p.to_plant = tp.asInt();
  if (p.to_plant == 1 && p.scheme == 0) {
    return "FATAL: model_settings.bnf_to_plant = 1 needs bnf_scheme > 0";
  }
  if (p.pro_c < 0 || p.pro_d > 0 || p.aet_a < 0 || p.cnst < 0 || p.asym < 0 || p.avln_cap < 0) {
    return "FATAL: model_settings BNF parameters: need bnf_pro_c >= 0, bnf_pro_d <= 0, "
           "bnf_aet_a >= 0, bnf_const >= 0, bnf_asym >= 0, bnf_avln_cap >= 0";
  }
  if (p.scheme == 0 && (p.asym != 0.0 || p.cnst != 0.0 || p.avln_cap != 0.0)) {
    return "FATAL: model_settings bnf_asym / bnf_const / bnf_avln_cap need bnf_scheme > 0";
  }
  if (p.scheme > 0 && n_rootless_fix != 2) {
    return "FATAL: model_settings.bnf_scheme > 0 needs n_rootless_fix = 2 (BNF enters AVLN of "
           "the top soil layer, also when it has no roots)";
  }
  return "";
}

// PRO, Meyerholt et al. 2016 Eq. A2 (NPP < 0 treated as 0).
inline double pro(double npp, double c, double d) {
  return c * (1.0 - std::exp(d * std::fmax(npp, 0.0)));
}
// AET, Meyerholt et al. 2016 Eq. A1, clamped at 0.
inline double aet(double et, double a, double b) {
  return std::fmax(0.0, a * et + b);
}
// Annual symbiotic-type rate of the scheme (without the free-living term).
inline double scheme_annual(const Params& p, double npp_yr, double et_yr) {
  switch (p.scheme) {
    case 1: return pro(npp_yr, p.pro_c, p.pro_d);
    case 2: return aet(et_yr, p.aet_a, p.aet_b);
    case 3: return p.cnst;
    default: return 0.0;
  }
}
// Taper on soil inorganic N (cap <= 0: none).
inline double cap_factor(double avln, double cap) {
  if (!(cap > 0.0)) { return 1.0; }
  return std::fmax(0.0, 1.0 - std::fmax(avln, 0.0) / cap);
}
// Month weight from the previous year's monthly NPP (positive part).
inline double month_weight(const double npp12[12], int m) {
  double s = 0.0;
  for (int k = 0; k < 12; ++k) { s += std::fmax(npp12[k], 0.0); }
  if (!(s > 0.0)) { return 1.0 / 12.0; }
  return std::fmax(npp12[m], 0.0) / s;
}

// Previous-year memory (12 monthly values of the current year, rolled after month 11).
struct YearMemory {
  double cur_npp[12] = {0}, cur_et[12] = {0};
  double prev_npp[12] = {0}, prev_et[12] = {0};
  double cur_pft_npp[MAX_PFT] = {0}, prev_pft_npp[MAX_PFT] = {0};  // Phase 2 (annual sums)
  bool have_prev = false;
  // Phase 2: add PFT ip's NPP of this month (call before record() of the month).
  void record_pft(int ip, double npp) {
    if (ip < 0 || ip >= MAX_PFT) { return; }
    cur_pft_npp[ip] += npp;
  }
  void record(int m, double npp, double et) {
    if (m < 0 || m > 11) { return; }
    cur_npp[m] = npp; cur_et[m] = et;
    if (m == 11) {
      for (int k = 0; k < 12; ++k) {
        prev_npp[k] = cur_npp[k]; prev_et[k] = cur_et[k];
        cur_npp[k] = 0.0; cur_et[k] = 0.0;
      }
      for (int k = 0; k < MAX_PFT; ++k) { prev_pft_npp[k] = cur_pft_npp[k]; cur_pft_npp[k] = 0.0; }
      have_prev = true;
    }
  }
  double prev_npp_yr() const { double s = 0; for (int k = 0; k < 12; ++k) { s += prev_npp[k]; } return s; }
  double prev_et_yr()  const { double s = 0; for (int k = 0; k < 12; ++k) { s += prev_et[k]; } return s; }
};

// BNF for month m (gN/m2/month): 0 without a previous year.
inline double month_bnf(const Params& p, const YearMemory& y, int m, double avln_tot) {
  if (p.scheme <= 0 || !y.have_prev) { return 0.0; }
  double ann = scheme_annual(p, y.prev_npp_yr(), y.prev_et_yr()) + p.asym;
  return ann * month_weight(y.prev_npp, m) * cap_factor(avln_tot, p.avln_cap);
}

// Phase 2: the two parts of month m (gN/m2/month): sym = scheme part, asym = free-living.
inline void month_parts(const Params& p, const YearMemory& y, int m, double avln_tot,
                        double& sym, double& asym) {
  sym = 0.0; asym = 0.0;
  if (p.scheme <= 0 || !y.have_prev) { return; }
  double f = month_weight(y.prev_npp, m) * cap_factor(avln_tot, p.avln_cap);
  sym = scheme_annual(p, y.prev_npp_yr(), y.prev_et_yr()) * f;
  asym = p.asym * f;
}
// Phase 2: PFT shares of the plant BNF: previous-year positive annual NPP among the
// active PFTs (equal if all <= 0). Returns the number of active PFTs (0: all shares 0).
inline int pft_shares(const double prev_pft_npp[], const bool active[], int n, double out[]) {
  int na = 0; double s = 0.0;
  for (int i = 0; i < n; ++i) {
    out[i] = 0.0;
    if (active[i]) { ++na; s += std::fmax(prev_pft_npp[i], 0.0); }
  }
  if (na == 0) { return 0; }
  for (int i = 0; i < n; ++i) {
    if (!active[i]) { continue; }
    out[i] = (s > 0.0) ? std::fmax(prev_pft_npp[i], 0.0) / s : 1.0 / na;
  }
  return na;
}

} // namespace bnf
#endif
