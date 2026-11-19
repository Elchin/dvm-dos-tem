#ifndef MOSSDYN_H_
#define MOSSDYN_H_
// model_settings.moss_dyn: dynamic live-moss layer thickness (default 0 = legacy static 0.01 m,
// bit-identical). See docs/moss_dyn_worknotes.md.
//   0  legacy: one moss layer, 0.01 m hard-coded (Ground.cpp initSnowSoilLayers / redivideMossLayers)
//   dmax: model_settings.moss_dmax if > 0, else the DOS-TEM maximum by drainage class (default_dmax):
//        4 cm upland (drainage class 0, well drained), 6 cm lowland (1, poorly drained); Jafarov et al. 2013
//        (sec. 2, after Eq. 1) for the DOS-TEM regrowth curve of Yi et al. 2009b, "according to our observations"
//        at their Interior Alaska upland/lowland sites. Not the cmt_dimground maxdmossthick (old dead-moss value).
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
  double dmax      = -1.0;   // moss_dmax (m); < 0 = default_dmax(drainage class)
  double dmin      = 0.001;  // moss_min_thick (m), > 0
  double t_half    = 5.0;    // moss_t_half (yr), >= 0 (0 = always dmax)
  double bulkden   = 25.0;   // moss_bulkden (kg dry mass / m3), > 0
  double cfrac     = 0.45;   // moss_cfrac (g C / g dry mass), (0, 1]
  // model_settings.moss_tc_scheme (moss layer thermal conductivity, any moss_dyn mode):
  //   0 (default) legacy (SoilLayer Kersten form, or org_tc_scheme 1 if on)
  //   1 Chadburn et al. 2015 (GMD 8, 1493, sec. 2.2.1; doi:10.5194/gmd-8-1493-2015): linear in moss
  //     water content between dry 0.06 and saturated 0.5 W/m/K (Soudzilovskaia et al. 2013 values),
  //     k = k_dry + (k_sat - k_dry) * theta/theta_sat, used for the frozen and unfrozen moss layer.
  //     Overrides org_tc_scheme for the moss layer only (fibric/humic layers keep their scheme).
  int    tc_scheme = 0;      // moss_tc_scheme
  double tc_kdry   = 0.06;   // moss_tc_kdry (W/m/K), Chadburn et al. 2015
  double tc_ksat   = 0.5;    // moss_tc_ksat (W/m/K), Chadburn et al. 2015
  //   2 Johansen (1975) / Farouki (1981) form as in CLM (Lawrence & Slater 2008) and GIPL (Nicolsky et al.
  //     2007 Eq. 6): saturated thawed/frozen k = ks^(1-eta) * k_liq^eta / ks^(1-eta) * k_ice^eta (eta = moss
  //     porosity), unsaturated k = Ke*k_sat + (1-Ke)*k_dry with Ke = max(0, log10(Sr)+1) thawed, Ke = Sr frozen,
  //     Sr = (theta_liq+theta_ice)/eta. The enthalpy solver blends the two as lam_u^(1-f) lam_f^f (f = ice/W).
  //     k_dry = moss_tc_kdry. Like scheme 1 it replaces org_tc_scheme (S21H) for the moss layer only.
  double tc_ksolid = 0.25;   // moss_tc_ksolid (W/m/K): organic matter solids, Farouki 1981 (CLM5, Lawrence & Slater 2008)
  // model_settings.moss_water_scheme:
  //   0 (default) legacy: the moss liquid is overwritten daily by the first non-moss layer's liquid
  //     (Soil_Env::updateDailySM), the moss holds no ice of its own.
  //   1 the moss layer keeps its own liquid + ice. Daily (Soil_Env::mossWaterStep): it intercepts
  //     infiltration up to moss_thetasat, drains liquid above the Chadburn et al. 2015 Eq. 1 equilibrium
  //     theta_eq = theta_sat (psi / psi_sat)^(-1/b) (psi = suction of the top soil layer; Brooks-Corey,
  //     Beringer et al. 2001 moss values) into the soil infiltration, is refilled from the top soil layer up
  //     to theta_eq, and supplies soil evaporation first. Freezing by the enthalpy solver (thermal_solver 1).
  //     Scheme 1 also forces richards_infil_fix = 1 and richards_dt_fix = 1 (Cohort.cpp): the moss drainage
  //     reaches Richards in bursts and on partly frozen columns, where the legacy B2 / R1 bugs lose / create
  //     about 1-2 mm/yr (docs/moss_dyn_worknotes.md section 9).
  int    water_scheme = 0;   // moss_water_scheme
  double w_bsw      = 1.0;   // moss_bsw (-), Beringer et al. 2001 via Chadburn et al. 2015
  double w_psisat   = 0.12;  // moss_psisat (m suction), Beringer et al. 2001 via Chadburn et al. 2015
  // moss_thetasat (-): -1 (default) = the moss layer porosity (Elchin 2026-11-18: cap matches porosity);
  // a value in (0, 1] is used instead, capped at porosity (Beringer et al. 2001 via Chadburn et al. 2015: 0.9).
  // Frozen moss (frozen == 1): today's surface water (infiltration, then the surface runoff split off for frozen
  // ground) is taken and refrozen up to the pore space and the cold content (refreeze_take); the rest goes on to
  // the existing frozen-ground path (ponding / surface runoff).
  double w_thetasat = -1.0;
};

// Chadburn et al. 2015 Eq. 1 (Brooks-Corey): equilibrium moss water content for a soil suction psi (m, >= 0).
inline double theta_eq(double psi, double thetasat, double psisat, double b) {
  if (!(psi > psisat)) { return thetasat; }
  return thetasat * std::pow(psi / psisat, -1.0 / b);
}

// One day of moss_water_scheme 1 (all mm = kg/m2). Pure function so the mass balance is testable.
struct WaterIn  { double liq, ice, dz, poro, thetasat, theq, minliq, infil, evap, fib_avail; bool active; };
struct WaterOut { double intercept = 0.0, drain = 0.0, rise = 0.0, evap = 0.0, liq = 0.0; };
inline WaterOut water_step(const WaterIn& w) {
  WaterOut o; o.liq = w.liq;
  if (!w.active || !(w.dz > 0.0)) { return o; }
  const double icevol = std::fmax(w.ice, 0.0) * (1000.0 / 917.0);            // mm of pore volume taken by ice
  const double tsat = std::fmin(w.thetasat, w.poro);
  const double cap = std::fmax(0.0, tsat * w.dz * 1000.0 - icevol - o.liq);   // room up to theta_sat
  o.intercept = std::fmin(std::fmax(w.infil, 0.0), cap);
  o.liq += o.intercept;
  // liquid target; residual water (minliq) does not drain
  const double eq = std::fmax(std::fmax(0.0, w.minliq), std::fmin(w.theq, tsat) * w.dz * 1000.0 - icevol);
  if (o.liq > eq) {
    o.drain = o.liq - eq; o.liq = eq;
  } else if (w.fib_avail > 0.0) {
    o.rise = std::fmin(eq - o.liq, w.fib_avail); o.liq += o.rise;
  }
  o.evap = std::fmin(std::fmax(w.evap, 0.0), std::fmax(0.0, o.liq - w.minliq));
  o.liq -= o.evap;
  return o;
}

// moss_water_scheme 1, frozen moss: water (kg/m2) the moss takes from supply and refreezes today, limited by
// the pore space left (as ice, 917 kg/m3) and the cold content hcap (-T) / (Lf + h_in); hcap in J/m2/K, T in C,
// h_in = specific enthalpy of the incoming liquid (J/kg, >= 0).
inline double refreeze_take(double supply, double dz, double poro, double liq, double ice, double hcap, double T,
                            double h_in, double lf) {
  if (!(supply > 0.0) || !(T < 0.0) || !(dz > 0.0) || !(hcap > 0.0)) { return 0.0; }
  const double room_mm = std::fmax(0.0, poro * dz * 1000.0 - std::fmax(ice, 0.0) * (1000.0 / 917.0) - std::fmax(liq, 0.0));
  const double cold = hcap * (-T) / (lf + std::fmax(h_in, 0.0));
  return std::fmax(0.0, std::fmin(supply, std::fmin(room_mm * 0.917, cold)));
}

// moss_tc_scheme 2: Johansen / Farouki moss conductivity (see Params). sr = (theta_liq+theta_ice)/eta.
inline double johansen_k(double sr, double eta, double ks, double kliq, double kice, double kdry, bool frozen) {
  sr = std::fmin(std::fmax(sr, 0.0), 1.0);
  if (sr < 1.0e-7) { return kdry; }
  const double ksat = std::pow(ks, 1.0 - eta) * std::pow(frozen ? kice : kliq, eta);
  const double ke = frozen ? sr : std::fmax(0.0, std::log10(sr) + 1.0);
  return ke * ksat + (1.0 - ke) * kdry;
}

// DOS-TEM maximum live-moss thickness (m) by drainage class: 1 (poorly drained) = lowland 0.06,
// anything else (0 = well drained) = upland 0.04. Jafarov et al. 2013, doi:10.1088/1748-9326/8/3/035030.
const double DMAX_UPLAND = 0.04, DMAX_LOWLAND = 0.06;
inline double default_dmax(int drainage_type) { return drainage_type == 1 ? DMAX_LOWLAND : DMAX_UPLAND; }

// Chadburn et al. 2015 moss conductivity; sat = theta / theta_sat (clamped to [0, 1]).
inline double chadburn_k(double sat, double kdry, double ksat) {
  sat = std::fmin(std::fmax(sat, 0.0), 1.0);
  return kdry + (ksat - kdry) * sat;
}

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
    return "model_settings.moss_dmax must be < 0 (DOS-TEM default by drainage class) or in [moss_min_thick, 1] m";
  }
  if (!(p.t_half >= 0.0 && p.t_half <= 1000.0)) { return "model_settings.moss_t_half must be in [0, 1000] yr"; }
  if (!(p.bulkden > 0.0 && p.bulkden <= 500.0)) { return "model_settings.moss_bulkden must be in (0, 500] kg/m3"; }
  if (!(p.cfrac > 0.0 && p.cfrac <= 1.0)) { return "model_settings.moss_cfrac must be in (0, 1]"; }
  Json::Value t = ms.get("moss_tc_scheme", 0);
  if (!num_(t) || t.asDouble() != std::floor(t.asDouble()) || t.asInt() < 0 || t.asInt() > 2) {
    return "model_settings.moss_tc_scheme must be 0 (legacy, default), 1 (Chadburn et al. 2015 linear) or 2 (Johansen/Farouki geometric mean)";
  }
  p.tc_scheme = t.asInt();
  K tk[] = {{"moss_tc_kdry", &p.tc_kdry}, {"moss_tc_ksat", &p.tc_ksat}, {"moss_tc_ksolid", &p.tc_ksolid},
            {"moss_bsw", &p.w_bsw}, {"moss_psisat", &p.w_psisat}, {"moss_thetasat", &p.w_thetasat}};
  for (auto& k : tk) {
    if (ms.isMember(k.k)) {
      if (!num_(ms[k.k])) { return std::string("model_settings.") + k.k + " must be a finite number"; }
      *k.d = ms[k.k].asDouble();
    }
  }
  if (!(p.tc_kdry > 0.0 && p.tc_kdry <= p.tc_ksat && p.tc_ksat <= 5.0)) {
    return "model_settings.moss_tc_kdry / moss_tc_ksat must satisfy 0 < kdry <= ksat <= 5 W/m/K";
  }
  if (!(p.tc_ksolid > 0.0 && p.tc_ksolid <= 5.0)) { return "model_settings.moss_tc_ksolid must be in (0, 5] W/m/K"; }
  Json::Value w = ms.get("moss_water_scheme", 0);
  if (!num_(w) || w.asDouble() != std::floor(w.asDouble()) || w.asInt() < 0 || w.asInt() > 1) {
    return "model_settings.moss_water_scheme must be 0 (legacy copy of the fibric liquid, default) or 1 (own moss water and ice)";
  }
  p.water_scheme = w.asInt();
  if (!(p.w_bsw > 0.0 && p.w_bsw <= 20.0)) { return "model_settings.moss_bsw must be in (0, 20]"; }
  if (!(p.w_psisat > 0.0 && p.w_psisat <= 10.0)) { return "model_settings.moss_psisat must be in (0, 10] m"; }
  if (!(p.w_thetasat == -1.0 || (p.w_thetasat > 0.0 && p.w_thetasat <= 1.0))) {
    return "model_settings.moss_thetasat must be -1 (layer porosity, default) or in (0, 1]";
  }
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
