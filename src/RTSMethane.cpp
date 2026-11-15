#include "../include/RTSMethane.h"
#include <algorithm>
#include <cmath>
#include <json/value.h>

void RTSMethaneParams::load(const Json::Value & ms) {
  struct { const char* k; double* v; } m[] = {
    {"thermokarst_ch4_frac", &ch4_frac}, {"thermokarst_ch4_rel_q10", &rel_q10},
    {"thermokarst_ch4_tref", &tref}, {"thermokarst_ch4_plant_frac", &plant_frac},
    {"thermokarst_ch4_lai_ref", &lai_ref},
    {"thermokarst_ch4_oxid_min", &oxid_min}, {"thermokarst_ch4_oxid_max", &oxid_max},
    {"thermokarst_ch4_oxid_zscale", &oxid_zscale},
    {"thermokarst_ch4_pond_sat0", &pond_sat0}, {"thermokarst_ch4_pond_tau", &pond_tau},
    {"thermokarst_ch4_lat_frac", &lat_ch4_frac},
    {"thermokarst_ch4_mg0", &mg0}, {"thermokarst_ch4_mg_tau", &mg_tau}, {"thermokarst_ch4_bg_sat", &bg_sat} };
  for (auto & e : m) { if (ms.isMember(e.k)) { *e.v = ms[e.k].asDouble(); } }
}

double RTSMethane::saturatedFraction(const double ztop, const double dz, const double wtd) {
  if (dz <= 0.0) { return 0.0; }
  return std::max(0.0, std::min(1.0, ((ztop + dz) - wtd) / dz));
}

double RTSMethane::oxidFraction(const RTSMethaneParams & p, const double wtd) {
  const double z = std::max(0.0, wtd);
  const double f = p.oxid_min + (p.oxid_max - p.oxid_min) * (1.0 - std::exp(-z / std::max(1e-9, p.oxid_zscale)));
  return std::max(0.0, std::min(1.0, f));
}

double RTSMethane::pondSaturation(const RTSMethaneParams & p, const int months_since_init) {
  if (months_since_init < 0 || p.pond_sat0 <= 0.0) { return 0.0; }
  const double tau = std::max(1e-9, 12.0 * p.pond_tau);
  return std::max(0.0, std::min(1.0, p.pond_sat0 * std::exp(-months_since_init / tau)));
}

double RTSMethane::methanogenFactor(const RTSMethaneParams & p, const int months_since_init) {
  if (months_since_init < 0) { return 1.0; }
  const double m0 = std::max(0.0, std::min(1.0, p.mg0));
  const double tau = std::max(1e-9, 12.0 * p.mg_tau);
  return 1.0 - (1.0 - m0) * std::exp(-months_since_init / tau);
}

RTSMethaneFluxes RTSMethane::compute(const std::vector<RTSMethaneLayer> & col,
                                     const RTSMethaneParams & p, const double wtd,
                                     const double pond, const double lai) {
  RTSMethaneFluxes f;
  for (size_t i = 0; i < col.size(); i++) {
    const RTSMethaneLayer & l = col[i];
    if (l.rh <= 0.0 || l.tem <= 0.0) { continue; }
    double sat = saturatedFraction(l.ztop, l.dz, wtd);
    sat += (1.0 - sat) * pond;
    const double tf = std::pow(p.rel_q10, (l.tem - p.tref) / 10.0);
    // CH4 can never exceed the anaerobic RH it is partitioned from
    f.prod += std::max(0.0, std::min(1.0, p.ch4_frac * tf)) * sat * l.rh;
  }
  const double fl = std::max(0.0, std::min(1.0, lai / std::max(1e-9, p.lai_ref)));
  f.plant = f.prod * std::max(0.0, std::min(1.0, p.plant_frac)) * fl;
  f.oxid = (f.prod - f.plant) * oxidFraction(p, std::max(0.0, wtd) * (1.0 - pond));
  f.emis = f.prod - f.oxid;
  return f;
}
