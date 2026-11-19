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
    {"thermokarst_ch4_mg0", &mg0}, {"thermokarst_ch4_mg_tau", &mg_tau}, {"thermokarst_ch4_bg_sat", &bg_sat},
    // scheme (b)
    {"thermokarst_ch4_kfrac", &kfrac}, {"thermokarst_ch4_prodq10", &prodq10},
    {"thermokarst_ch4_avail_frac", &avail_frac},
    {"thermokarst_ch4_oxid_vmax", &oxid_vmax}, {"thermokarst_ch4_oxid_km", &oxid_km},
    {"thermokarst_ch4_oxidq10", &oxidq10}, {"thermokarst_ch4_plant_k", &plant_k},
    {"thermokarst_ch4_froot_ref", &froot_ref},
    {"thermokarst_ch4_rhizo_ox", &rhizo_ox}, {"thermokarst_ch4_tau_sat", &tau_sat},
    {"thermokarst_ch4_tau_unsat", &tau_unsat} };
  for (auto & e : m) { if (ms.isMember(e.k)) { *e.v = ms[e.k].asDouble(); } }
}

void RTSMethaneParams::applySchemeBDefaults(const Json::Value & ms) {
  if (!ms.isMember("thermokarst_ch4_pond_sat0")) { pond_sat0 = 0.5; }
  if (!ms.isMember("thermokarst_ch4_pond_tau"))  { pond_tau = 10.0; }
  // kfrac recalibrated on the merged base (moss/rts-int-m2d, docs/moss_dyn_worknotes.md s.12) with the
  // rts-ch4-b protocol: CMT06 wet sedge, Gydan cell, bg_sat = 1, baseline two-stage run (hist 76, sc 30),
  // mean annual TKCH4EMIS over stage 2b = 5.0 gC/m2/yr -> 0.114 (monthly2daily_mean_preserving off, the
  // default; 0.109 with it on). The prototype value 0.29 (struct default) gives 11.9 on this base.
  if (!ms.isMember("thermokarst_ch4_kfrac"))     { kfrac = 0.114; }
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

/** Scheme (b) (merged from rts-ch4-b 6d49ad57, function body unchanged). */
RTSMethaneFluxes RTSMethane::computeB(std::vector<RTSMethaneLayer> & col, const RTSMethaneParams & p,
                                      const double wtd, const double pond, const double lai,
                                      const double days) {
  RTSMethaneFluxes f;
  const double wtd_eff = std::max(0.0, wtd) * (1.0 - pond);
  const double fox_transit = oxidFraction(p, wtd_eff);
  const double flai = std::max(0.0, std::min(1.0, lai / std::max(1e-9, p.lai_ref)));
  for (size_t i = 0; i < col.size(); i++) {
    RTSMethaneLayer & l = col[i];
    for (int k = 0; k < 4; k++) { l.dc[k] = 0.0; }
    double P = std::max(0.0, l.pool);
    if (l.tem > 0.0) {
      double sat = saturatedFraction(l.ztop, l.dz, wtd);
      sat += (1.0 - sat) * pond;
      // production
      double pot = 0.0, ctot = 0.0;
      for (int k = 0; k < 4; k++) {
        if (l.c[k] > 0.0 && l.kd[k] > 0.0) { pot += l.kd[k] * l.c[k]; }
        if (l.c[k] > 0.0) { ctot += l.c[k]; }
      }
      double prod = std::max(0.0, p.kfrac) * sat * std::pow(p.prodq10, (l.tem - p.tref) / 10.0) * pot;
      // split over the SOM pools in proportion to kd_k C_k; each pool gives at most
      // avail_frac of its C per month (so prod <= avail_frac * layer C <= available C)
      const double af = std::max(0.0, std::min(1.0, p.avail_frac));
      double ptot = 0.0;
      if (pot > 0.0 && prod > 0.0) {
        for (int k = 0; k < 4; k++) {
          if (l.c[k] > 0.0 && l.kd[k] > 0.0) {
            l.dc[k] = std::min(prod * (l.kd[k] * l.c[k]) / pot, af * l.c[k]);
            ptot += l.dc[k];
          }
        }
      }
      prod = ptot;
      P += prod; f.prod += prod;
      // in-layer oxidation (oxic part of the layer)
      const double kmz = std::max(1e-12, p.oxid_km * l.dz);
      double ox = (1.0 - sat) * p.oxid_vmax * l.dz * days * P / (kmz + P)
                  * std::pow(p.oxidq10, (l.tem - p.tref) / 10.0);
      ox = std::max(0.0, std::min(P, ox));
      P -= ox; f.oxid += ox;
      // plant-mediated transport (aerenchyma), part oxidized in the rhizosphere
      const double fr = std::max(0.0, std::min(1.0, l.froot / std::max(1e-9, p.froot_ref)));
      const double pt = P * (1.0 - std::exp(-std::max(0.0, p.plant_k) * flai * fr * days));
      P -= pt;
      const double rz = std::max(0.0, std::min(1.0, p.rhizo_ox));
      f.oxid += rz * pt; f.emis += (1.0 - rz) * pt; f.plant += (1.0 - rz) * pt;
      // escape by diffusion/ebullition; from below the WT through the oxic zone
      const double tau = (sat >= 0.5) ? p.tau_sat : p.tau_unsat;
      const double esc = P * (1.0 - std::exp(-days / std::max(1e-9, tau)));
      P -= esc;
      const double oxt = esc * sat * fox_transit;
      f.oxid += oxt; f.emis += esc - oxt;
    }
    l.pool = P; f.pool += P;
  }
  return f;
}
