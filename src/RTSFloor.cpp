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
    {"thermokarst_floor_anox_f", &anox_f}, {"rts_floor_rmin0", &rmin0}, {"rts_floor_rmin_tau", &rmin_tau} };
  if (ms.isMember("thermokarst_floor_rmin")) { rmin_on = ms["thermokarst_floor_rmin"].asBool(); }
  if (ms.isMember("rts_floor_rmin_years")) { rmin_years = std::max(1, std::min(10, ms["rts_floor_rmin_years"].asInt())); }
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

double RTSFloor::rminFactor(const RTSFloorParams & p, const int m) {
  if (m < 0) { return 1.0; }
  const double r0 = std::max(0.0, std::min(1.0, p.rmin0));
  return 1.0 - (1.0 - r0) * std::exp(-m / std::max(1e-9, 12.0 * p.rmin_tau));
}

double RTSFloor::applyMinimum(const double target, const int n, const double * w,
                              const double * crawc, const double * csoma, const double * csompr, const double * csomcr,
                              double * rrawc, double * rsoma, double * rsompr, double * rsomcr) {
  double cur = 0.0, wsum = 0.0;
  for (int i = 0; i < n; i++) {
    cur += rrawc[i] + rsoma[i] + rsompr[i] + rsomcr[i];
    const double c = std::max(0.0, crawc[i]) + std::max(0.0, csoma[i]) + std::max(0.0, csompr[i]) + std::max(0.0, csomcr[i]);
    wsum += std::max(0.0, w[i]) * c;
  }
  double need = target - cur;
  if (!(need > 0.0) || !(wsum > 0.0)) { return 0.0; }
  double added = 0.0;
  for (int i = 0; i < n; i++) {
    double cp[4] = { std::max(0.0, crawc[i]), std::max(0.0, csoma[i]), std::max(0.0, csompr[i]), std::max(0.0, csomcr[i]) };
    double * rp[4] = { &rrawc[i], &rsoma[i], &rsompr[i], &rsomcr[i] };
    const double c = cp[0] + cp[1] + cp[2] + cp[3];
    if (c <= 0.0 || w[i] <= 0.0) { continue; }
    double x = need * std::max(0.0, w[i]) * c / wsum;
    // keep the layer's total RH within rmin_maxfrac of its C
    const double rl = *rp[0] + *rp[1] + *rp[2] + *rp[3];
    x = std::min(x, std::max(0.0, rmin_maxfrac * c - rl));
    if (x <= 0.0) { continue; }
    for (int k = 0; k < 4; k++) { *rp[k] += x * cp[k] / c; }
    added += x;
  }
  return added;
}
