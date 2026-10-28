#ifndef DAILYCSPLIT_H_
#define DAILYCSPLIT_H_
// v5p: statistical daily split of the monthly carbon fluxes (output only).
// model_settings.daily_c_split (0 = off, default; 1 = on). The model still
// computes every C flux once per month; with the switch on, each monthly
// total M is divided over the n days of the calendar month as
//     X_d = M * w_d / sum_k w_k
// where w_d is the model's own environmental scalar for the flux evaluated
// with that day's drivers (sum_d X_d == M up to floating-point round-off).
// Nothing computed here is fed back into the model state, the restart files
// or the monthly/yearly outputs. Pure helpers (no model state) so they can be
// unit tested; see tests/unit/test_daily_c_split.cpp.
// Refs: no published source for this exact method. The idea of spreading a
//   coarse total over finer steps with a driver-based weight follows
//   Olsen & Randerson 2004, doi:10.1029/2003JD003968 (3-hourly GPP/respiration
//   from monthly CASA fluxes with radiation/temperature scalars) and
//   Fisher et al. 2016, doi:10.5194/bg-13-4271-2016 (3-hourly downscaling of
//   monthly MsTMIP model NEE). Exact preservation of the coarse
//   total (temporal disaggregation; principle only): Denton 1971,
//   doi:10.1080/01621459.1971.10482227; Chow & Lin 1971, doi:10.2307/1928739.
//   See docs/daily-mode-REFERENCES.md.
#include <cmath>
#include <string>
#include <json/json.h>

namespace dcsplit {

// Daily-only output variables (ecosystem totals, g/m2/day), in output order.
enum Var { GPP = 0, NPP, RM, RG, RHSOM, RHDWD, RECO, NEE, INGPP, INNPP, NVAR };
static const char* const VAR_NAMES[NVAR] = {
  "GPP_DS", "NPP_DS", "RM_DS", "RG_DS", "RHSOM_DS", "RHDWD_DS",
  "RECO_DS", "NEE_DS", "INGPP_DS", "INNPP_DS"
};
// Monthly diagnostic: bit flags of the fluxes that used a fallback weight.
static const char* const FALLBACK_NAME = "DS_FALLBACK";
const int FLAG_GPP    = 1;    // GPP/INGPP weights all zero: uniform M/n
const int FLAG_RM     = 2;    // RM weights all zero: uniform
const int FLAG_RG_GPP = 4;    // RG weight max(GPP-RM,0) all zero: GPP weight used
const int FLAG_RG_UNI = 8;    // RG (or RG*) weight and GPP weight all zero: uniform
const int FLAG_RHSOM  = 16;   // a soil layer with RH != 0 had all-zero weights: uniform
const int FLAG_RHDWD  = 32;   // wood-debris weights all zero: uniform
const int FLAG_LAYER_REMAP = 64; // not a fallback: Env-loop and Bgc soil layers differed,
                                 // Bgc layers mapped to Env layers by midpoint depth

// Validation of model_settings.daily_c_split: only the numbers 0 and 1 are
// accepted (booleans, strings, 0.5, 2, -1, ... are rejected -> FATAL).
inline bool valid_switch(const Json::Value& v) {
  if (v.isBool() || !v.isNumeric()) { return false; }
  double d = v.asDouble();
  return d == 0.0 || d == 1.0;
}

// True for the variables that only exist with daily_c_split = 1.
inline bool is_ds_name(const std::string& name) {
  if (name == FALLBACK_NAME) { return true; }
  for (int k = 0; k < NVAR; ++k) { if (name == VAR_NAMES[k]) { return true; } }
  return false;
}

// Empty if the output-spec request is allowed, else the FATAL message.
// requested: any timestep ticked; daily: the Daily column ticked.
inline std::string request_error(const std::string& name, bool requested,
                                 bool daily, int daily_c_split) {
  if (!requested || !is_ds_name(name)) { return ""; }
  if (daily_c_split != 1) {
    return "FATAL: output " + name + " requested in output_spec.csv but model_settings.daily_c_split is 0 (set it to 1)";
  }
  if (name != FALLBACK_NAME && !daily) {
    return "FATAL: output " + name + " is daily only (tick Daily in output_spec.csv)";
  }
  if (name == FALLBACK_NAME && daily) {
    return "FATAL: output " + name + " is monthly only";
  }
  return "";
}

// A weight that is non-finite or negative counts as 0.
inline double clean(double w) {
  return (std::isfinite(w) && w > 0.0) ? w : 0.0;
}

// Split the monthly total M over n days with weights w (n >= 1).
// M == 0 -> all zero; sum of the cleaned weights == 0 (or not finite) with
// M != 0 -> uniform M/n and returns true (fallback used); otherwise
// out[d] = M*w_d/sum(w) and returns false.
inline bool split(double M, const double* w, int n, double* out) {
  if (M == 0.0) {
    for (int d = 0; d < n; ++d) { out[d] = 0.0; }
    return false;
  }
  double s = 0.0;
  for (int d = 0; d < n; ++d) { s += clean(w[d]); }
  if (!(s > 0.0) || !std::isfinite(s)) {
    for (int d = 0; d < n; ++d) { out[d] = M / n; }
    return true;
  }
  for (int d = 0; d < n; ++d) { out[d] = M * (clean(w[d]) / s); }
  return false;
}

// Split with a primary weight and a secondary one used only if the primary
// weights are all zero. Returns 0 (primary), 1 (secondary) or 2 (uniform).
inline int split2(double M, const double* w1, const double* w2, int n, double* out) {
  if (M == 0.0) { split(M, w1, n, out); return 0; }
  double s1 = 0.0;
  for (int d = 0; d < n; ++d) { s1 += clean(w1[d]); }
  if (s1 > 0.0 && std::isfinite(s1)) { split(M, w1, n, out); return 0; }
  return split(M, w2, n, out) ? 2 : 1;
}

// GPP weight of one PFT for one day: the day-varying part of
// Vegetation_Bgc::getGPP (cmax and foliage are constant within the month):
//   ci/(kc+ci) * par/(ki+par) * ftemp * thaw,  ci = co2*gv.
inline double gpp_weight(double co2, double gv, double kc, double par,
                         double ki, double ftemp, double thaw) {
  double ci = co2 * gv;
  double w = ci / (kc + ci) * (par / (ki + par)) * ftemp * thaw;
  return clean(w);
}

// PET floor used by Vegetation_Bgc::prepareIntegration before getGV.
inline double pet_floor(double pet) { return (pet <= 0.0) ? 0.0001 : pet; }

// RG weight of one PFT for one day: max(GPP_d - RM_d, 0) of the split fluxes.
inline double rg_weight(double gpp_d, double rm_d) {
  return clean(gpp_d - rm_d);
}

// Index of the layer (0..n-1, top-down, z = top depth, dz = thickness) that
// contains depth zm; above the column -> 0, below the column -> n-1.
inline int layer_containing(double zm, const double* z, const double* dz, int n) {
  if (n <= 0) { return -1; }
  for (int j = 0; j < n; ++j) {
    if (zm < z[j] + dz[j]) { return j; }
  }
  return n - 1;
}

// True if two layer structures are identical (count, top depths, thicknesses).
inline bool same_layers(int n1, const double* z1, const double* dz1,
                        int n2, const double* z2, const double* dz2) {
  if (n1 != n2) { return false; }
  for (int j = 0; j < n1; ++j) {
    if (z1[j] != z2[j] || dz1[j] != dz2[j]) { return false; }
  }
  return true;
}

// Monthly totals of one PFT (bd[ip], g/m2/month).
struct PftMonth { double gpp, rm, rg, ingpp, innpp; };

// Split one PFT's monthly fluxes with its daily GPP weights wg and RM weights
// wr and add the daily values to acc[var][day] (ecosystem totals). Returns
// the fallback flags.
//   GPP, INGPP: GPP weight; RM: RM weight;
//   RG: max(GPP_d - RM_d, 0), then the GPP weight, then uniform;
//   RG* = INGPP - RM - INNPP: the same RG weight;
//   NPP_d = GPP_d - RM_d - RG_d; INNPP_d = INGPP_d - RM_d - RG*_d.
inline int add_pft(const PftMonth& m, const double* wg, const double* wr,
                   int n, double acc[][31]) {
  int flags = 0;
  double g[31], r[31], rg[31], ig[31], rgs[31], wrg[31];
  if (split(m.gpp, wg, n, g))   { flags |= FLAG_GPP; }
  if (split(m.rm, wr, n, r))    { flags |= FLAG_RM; }
  for (int d = 0; d < n; ++d) { wrg[d] = rg_weight(g[d], r[d]); }
  int lv = split2(m.rg, wrg, wg, n, rg);
  if (lv == 1) { flags |= FLAG_RG_GPP; } else if (lv == 2) { flags |= FLAG_RG_UNI; }
  if (split(m.ingpp, wg, n, ig)) { flags |= FLAG_GPP; }
  double rgstar = m.ingpp - m.rm - m.innpp;
  lv = split2(rgstar, wrg, wg, n, rgs);
  if (lv == 1) { flags |= FLAG_RG_GPP; } else if (lv == 2) { flags |= FLAG_RG_UNI; }
  for (int d = 0; d < n; ++d) {
    acc[GPP][d]   += g[d];
    acc[RM][d]    += r[d];
    acc[RG][d]    += rg[d];
    acc[NPP][d]   += g[d] - r[d] - rg[d];
    acc[INGPP][d] += ig[d];
    acc[INNPP][d] += ig[d] - r[d] - rgs[d];
  }
  return flags;
}

// Split one soil layer's monthly RH (or the wood-debris RH) with the daily
// weights w (rhmoist*rhq10 of that day) and add it to acc_d. Returns true if
// the uniform fallback was used.
inline bool add_rh(double M, const double* w, int n, double* acc_d) {
  double x[31];
  bool fb = split(M, w, n, x);
  for (int d = 0; d < n; ++d) { acc_d[d] += x[d]; }
  return fb;
}

// RECO_d = RHSOM_d + RHDWD_d + RM_d + RG_d (same terms as the monthly RECO
// output); NEE_d = RECO_d - GPP_d (decision 3 of the design: no NEE in the
// model, sign convention positive = source).
inline void finish(double acc[][31], int n) {
  for (int d = 0; d < n; ++d) {
    acc[RECO][d] = acc[RHSOM][d] + acc[RHDWD][d] + acc[RM][d] + acc[RG][d];
    acc[NEE][d]  = acc[RECO][d] - acc[GPP][d];
  }
}

// Monthly counters (Cohort::ds_fallback_months, ds_remap_months): a month
// counts as a fallback month if any fallback flag is set (FLAG_LAYER_REMAP is
// not a fallback and is counted separately).
inline void count_month(int flags, long& fallback_months, long& remap_months) {
  if (flags & ~FLAG_LAYER_REMAP) { fallback_months++; }
  if (flags & FLAG_LAYER_REMAP) { remap_months++; }
}

// Zero all daily values (all 31 slots).
inline void clear(double acc[][31]) {
  for (int k = 0; k < NVAR; ++k) { for (int d = 0; d < 31; ++d) { acc[k][d] = 0.0; } }
}

} // namespace dcsplit
#endif /*DAILYCSPLIT_H_*/
