#include "../include/RTSFloor.h"
#include <algorithm>
#include <cmath>
#include <json/value.h>

void RTSFloorParams::load(const Json::Value & ms) {
  if (ms.isMember("thermokarst_rh_cap")) { rh_cap = ms["thermokarst_rh_cap"].asBool(); }
  if (ms.isMember("thermokarst_floor_export")) { exp_on = ms["thermokarst_floor_export"].asBool(); }
  if (ms.isMember("thermokarst_floor_anox")) { anox = ms["thermokarst_floor_anox"].asBool(); }
  struct { const char* k; double* v; } m[] = {
    {"thermokarst_rh_cap0", &cap0}, {"thermokarst_rh_cap_inf", &cap_inf}, {"thermokarst_rh_cap_tau", &cap_tau},
    {"thermokarst_floor_lat0", &lat0}, {"thermokarst_floor_lat_tau", &lat_tau},
    {"thermokarst_floor_anox_f", &anox_f} };
  for (auto & e : m) { if (ms.isMember(e.k)) { *e.v = ms[e.k].asDouble(); } }
}

double RTSFloor::capRate(const RTSFloorParams & p, const int m) {
  const double e = std::exp(-std::max(0, m) / std::max(1e-9, 12.0 * p.cap_tau));
  return std::max(0.0, p.cap_inf + (p.cap0 - p.cap_inf) * e);
}

double RTSFloor::capLayerRH(const RTSFloorParams & p, const int m, const double rh, const double c) {
  if (m < 0 || rh <= 0.0) { return rh; }
  const double cap = capRate(p, m) / 12.0 * std::max(0.0, c);
  return std::min(rh, cap);
}

double RTSFloor::exportFraction(const RTSFloorParams & p, const int m) {
  if (m < 0) { return 0.0; }
  const double f = std::max(0.0, p.lat0) / 12.0 * std::exp(-m / std::max(1e-9, 12.0 * p.lat_tau));
  return std::min(1.0, f);
}

double RTSFloor::anoxFactor(const RTSFloorParams & p, const double sat) {
  const double s = std::max(0.0, std::min(1.0, sat));
  const double f = std::max(0.0, std::min(1.0, p.anox_f));
  return 1.0 - (1.0 - f) * s;
}
