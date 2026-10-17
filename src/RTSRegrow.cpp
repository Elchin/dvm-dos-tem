#include "../include/RTSRegrow.h"
#include <algorithm>
#include <cmath>
#include <json/value.h>

void RTSRegrowParams::load(const Json::Value & ms) {
  if (ms.isMember("thermokarst_regrow")) { on = ms["thermokarst_regrow"].asBool(); }
  if (ms.isMember("rts_regrow_n_taper")) { n_taper = ms["rts_regrow_n_taper"].asBool(); }
  if (ms.isMember("rts_regrow_navln_exempt")) { navln_exempt = ms["rts_regrow_navln_exempt"].asBool(); }
  struct { const char* k; double* v; } m[] = {
    {"rts_regrow_t_active", &t_active}, {"rts_regrow_f_reseed", &f_reseed}, {"rts_regrow_tau_ref", &tau_ref},
    {"rts_regrow_B0", &boost0}, {"rts_regrow_Tb", &boost_T}, {"rts_regrow_ol_rate1", &ol_rate1}, {"rts_regrow_ol_rate2", &ol_rate2},
    {"rts_regrow_ol_T", &ol_T}, {"rts_regrow_ol_cn", &ol_cn}, {"rts_regrow_ol_maxmult", &ol_maxmult}, {"rts_regrow_ol_frac_raw", &ol_frac_raw}, {"rts_regrow_n_rate1", &n_rate1},
    {"rts_regrow_n_rate2", &n_rate2}, {"rts_regrow_n_T", &n_T}, {"rts_regrow_n_taper_pow", &n_taper_pow} };
  for (auto & e : m) { if (ms.isMember(e.k)) { *e.v = ms[e.k].asDouble(); } }
}

int RTSRegrow::pftType(int ifwoody, int ifdeciwoody, int nonvascular) {
  if (nonvascular == 3) { return RG_LICHEN; }
  if (nonvascular == 1 || nonvascular == 2) { return RG_MOSS; }
  if (ifwoody) { return ifdeciwoody ? RG_SHRUB : RG_EVERGREEN; }
  return RG_GRAM;
}

double RTSRegrow::tau(const RTSRegrowParams & p, double gpp_kg) {
  if (!(gpp_kg > 0.0)) { return p.tau_max; }
  return std::max(p.tau_min, std::min(p.tau_max, 1.35 * std::pow(gpp_kg, -1.68)));
}
int RTSRegrow::stabMonth(const RTSRegrowParams & p) { return std::max(0, (int)std::lround(12.0 * p.t_active)); }
double RTSRegrow::tYears(const RTSRegrowParams & p, int months) { return (months - stabMonth(p)) / 12.0; }

double RTSRegrow::boost(const RTSRegrowParams & p, int type, double t) {
  if (type != RG_SHRUB || t < 0.0) { return 1.0; }
  return 1.0 + (p.boost0 - 1.0) * std::exp(-t / std::max(1e-9, p.boost_T));
}
double RTSRegrow::olRate(const RTSRegrowParams & p, double t) {
  if (t < 0.0) { return 0.0; }
  return std::max(0.0, t <= p.ol_T ? p.ol_rate1 : p.ol_rate2);
}
double RTSRegrow::nRate(const RTSRegrowParams & p, double t) {
  if (t < 0.0) { return 0.0; }
  return std::max(0.0, t <= p.n_T ? p.n_rate1 : p.n_rate2);
}
double RTSRegrow::taperFactor(const RTSRegrowParams & p, double gpp_ratio) {
  if (!p.n_taper || !(gpp_ratio > 0.0)) { return 1.0; }
  if (gpp_ratio >= 1.0) { return 0.0; }
  return 1.0 - std::pow(gpp_ratio, std::max(0.0, p.n_taper_pow));
}
double RTSRegrow::avlnSupplyMonth(const RTSRegrowParams & p, double t, double taper) {
  const double oln = p.ol_cn > 0.0 ? olRate(p, t) / p.ol_cn : 0.0;
  const double f = std::max(0.0, std::min(1.0, taper));
  return f * std::max(0.0, nRate(p, t) - oln) / 12.0;
}
