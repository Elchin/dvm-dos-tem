#ifndef MOSSDYN_H_
#define MOSSDYN_H_
// model_settings.moss_dyn: dynamic live-moss layer thickness (default 0 = legacy static 0.01 m,
// bit-identical). See docs/moss_dyn_worknotes.md.
//   0  legacy: one moss layer, 0.01 m hard-coded (Ground.cpp initSnowSoilLayers / redivideMossLayers)
//   1  prescribed regrowth from years since disturbance (no carbon):
//        d = max(dmin, dmax * t / (t + t_half))      (DOS-TEM, Yi et al. 2009, doi:10.1029/2008JG000841;
//                                                     Jafarov et al. 2013 Eq. 1, doi:10.1088/1748-9326/8/3/035030)
//        t = cd.yrsdist (reset to 0 by fire and RTS, in restart). t < 0 (never disturbed) or t_half = 0 -> dmax.
//   2  carbon-coupled: d = clamp(C_moss / (bulkden * cfrac * 1000), dmin, dmax)
//        C_moss = live C (leaf+stem+root, gC/m2 ground) of the nonvascular PFTs with vegcov > 0
//        (cmt_dimvegetation nonvasculartype > 0). Growth = the moss PFTs' own NPP, litter = their
//        litterfall into the first non-moss (fibric) layer, removal = the existing fire / RTS vegetation
//        burn/removal bookkeeping. No new carbon pool, so C balance is unchanged by construction.
//        Cohorts whose CMT has no nonvascular PFT keep the legacy layer (warning on stderr).
// Thickness changes are applied in January (Cohort::updateMonthly_DIMgrd) and right after a fire or
// RTS removal. Never zero (dmin > 0) because ground->fstmossl is dereferenced unconditionally.
// Moss C-density refs: O'Donnell et al. 2009 doi:10.1097/SS.0b013e3181c4a7f8 (live moss bulk density
// ~10-50 kg/m3); Druel et al. 2017 doi:10.5194/gmd-10-4693-2017 (h = biomass / density).
#include <cmath>
#include <string>
#include <json/json.h>

namespace mossdyn {

struct Params {
  int    mode      = 0;      // moss_dyn
  double dmax      = -1.0;   // moss_dmax (m); < 0 = cmt_dimground maxdmossthick
  double dmin      = 0.001;  // moss_min_thick (m), > 0
  double t_half    = 5.0;    // moss_t_half (yr), >= 0 (0 = always dmax)
  double bulkden   = 25.0;   // moss_bulkden (kg dry mass / m3), > 0
  double cfrac     = 0.45;   // moss_cfrac (g C / g dry mass), (0, 1]
};

inline bool num_(const Json::Value& v) { return !v.isBool() && v.isNumeric() && std::isfinite(v.asDouble()); }

// Parse model_settings; returns "" or the FATAL message. Absent keys = defaults.
inline std::string parse(const Json::Value& ms, Params& p) {
  p = Params();
  Json::Value m = ms.get("moss_dyn", 0);
  if (!num_(m) || m.asDouble() != std::floor(m.asDouble()) || m.asInt() < 0 || m.asInt() > 2) {
    return "model_settings.moss_dyn must be 0 (legacy static 0.01 m, default), 1 (regrowth curve) or 2 (moss PFT carbon)";
  }
  p.mode = m.asInt();
  struct K { const char* k; double* d; };
  K keys[] = {{"moss_dmax", &p.dmax}, {"moss_min_thick", &p.dmin}, {"moss_t_half", &p.t_half},
              {"moss_bulkden", &p.bulkden}, {"moss_cfrac", &p.cfrac}};
  for (auto& k : keys) {
    if (ms.isMember(k.k)) {
      if (!num_(ms[k.k])) { return std::string("model_settings.") + k.k + " must be a finite number"; }
      *k.d = ms[k.k].asDouble();
    }
  }
  if (!(p.dmin > 0.0 && p.dmin <= 0.05)) { return "model_settings.moss_min_thick must be in (0, 0.05] m"; }
  if (!(p.dmax < 0.0 || (p.dmax >= p.dmin && p.dmax <= 1.0))) {
    return "model_settings.moss_dmax must be < 0 (use cmt_dimground maxdmossthick) or in [moss_min_thick, 1] m";
  }
  if (!(p.t_half >= 0.0 && p.t_half <= 1000.0)) { return "model_settings.moss_t_half must be in [0, 1000] yr"; }
  if (!(p.bulkden > 0.0 && p.bulkden <= 500.0)) { return "model_settings.moss_bulkden must be in (0, 500] kg/m3"; }
  if (!(p.cfrac > 0.0 && p.cfrac <= 1.0)) { return "model_settings.moss_cfrac must be in (0, 1]"; }
  return "";
}

// Mode 1. ysf = years since disturbance (< 0 = never disturbed).
inline double regrow_thick(double dmax, double dmin, double ysf, double t_half) {
  double d = dmax;
  if (ysf >= 0.0 && t_half > 0.0) { d = dmax * ysf / (ysf + t_half); }
  return std::fmin(dmax, std::fmax(dmin, d));
}

// Mode 2. cmoss in gC/m2 ground; bulkden kg/m3; cfrac gC/g.
inline double carbon_thick(double cmoss, double bulkden, double cfrac, double dmin, double dmax) {
  double d = std::fmax(cmoss, 0.0) / (bulkden * 1000.0 * cfrac);
  return std::fmin(dmax, std::fmax(dmin, d));
}

} // namespace mossdyn
#endif
