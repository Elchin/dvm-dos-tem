#ifndef HGTRACER_H_
#define HGTRACER_H_
// Mercury (Hg) tracer that follows carbon (model_settings.hg_tracer, default 0 = off, nothing runs).
// Diagnostic only: it reads the C pools and C fluxes and never writes them, so C/N/water/energy are
// bit-identical with the switch on or off. See docs_src/WORK_NOTES_hg_tracer.md.
//
// Refs:
//   Schaefer et al. 2020, Nat. Commun. 11, 4650, doi:10.1038/s41467-020-18398-5 (SiBCASA-Hg: Hg moves
//     with every C transfer; of the Hg released by decay 16% evades as Hg(0), 1.5% is exported, the rest
//     re-binds to the soil organic matter of the same layer).
//   Schuster et al. 2018, GRL 45, 1463, doi:10.1002/2017GL075571 (RHgC median 1.6 ug Hg / g C).
//   T0 inventory (docs_src/HG_T0_inventory.md, Azaroff et al. 2026 synthesis, doi:10.17043/azaroff-2026-
//     permafrost-mercury-1): organic-horizon RHgC 0.32-0.35, mineral 0.97-1.69 (synthesis median 1.47).
//   Obrist et al. 2017, Nature 547, 201, doi:10.1038/nature22997; Olson et al. 2018, GBC 32, 1058 (tundra veg.
//     THg ~20-60 ng/g dry mass -> ~0.05-0.1 ug Hg / g C at 0.45 g C / g; doi:10.1029/2017GB005840).
//   Thompson et al. 2025, ES&T 59, 19457, doi:10.1021/acs.est.5c04510 (soil %MeHg 0.1 plateau, 0.7 bog,
//     2.0 fen) for the MeHg index (hg_mehg).
//
// Units: C in g C/m2, Hg and MeHg in ug/m2, ratios in ug Hg / g C, fluxes per model month.
//
// Pools (parallel to C): live vegetation per PFT x {leaf, stem, root}; standing dead per PFT; woody
// debris; one Hg pool per soil layer (pool-weighted: the 4 SOM pools of a layer share its Hg:C ratio).
//
// One generic step (hgtr::step) is applied around every model sub-step that moves C (environment with
// daily CH4, ground re-division, BGC, fire, RTS) from C snapshots before (a) and after (b) plus the
// sub-step's external C fluxes (Flux):
//   1. inputs: NPP > 0 of each PFT part carries r_leaf/r_stem/r_root (atmospheric Hg(0) uptake);
//      optional direct soil deposition hg_dep_soil (ug/m2/yr) to the top soil layer in the BGC step.
//   2. decomposition (RH, CH4 production, woody debris RH) of a pool releases C x (Hg:C of the pool);
//      f_evade of it evades, f_export is exported, the rest stays in the pool (re-binding);
//      DOC loss exports its Hg.
//   3. internal transfers: dC_int = C(b) - C(a) - inputs + outputs for each pool. Pools that lost C
//      (donors: litterfall, standing dead fall, burned / slumped material, ...) put Hg at their Hg:C ratio
//      into one pot; aggregate sinks (fire emission, RTS lateral export, RTS respiration) take their
//      C share of the pot; the rest goes to the pools that gained C (receivers) in proportion to dC_int.
//   4. soil geometry change (layer split/combine, burn, slump): the soil column is re-mapped bottom-up
//      on cumulative C (uniform Hg:C within an old layer); a net soil C loss is taken from the top, a net
//      gain is added on top as a receiver.
// The step is exactly conservative: total Hg change = uptake + deposition - evasion - export - fire.
// C that leaves a donor with no receiver or sink leaves its Hg behind (concentration rises), C that
// appears in a receiver with no donor brings no Hg (dilution); both are reported as unexplained C.
//
// MeHg index (hg_mehg = 1, needs hg_tracer = 1): per soil layer, no feedback on Hg or C:
//   P = k_meth x Hg x A x q10_meth^((T - t_ref)/10),   kd = k_demeth x q10_demeth^((T - t_ref)/10)
//   A = anaerobic C turnover fraction of the month: CH4 production / layer C where the CH4 module is on
//       (cd.ch4_on), else saturated fraction (water table) x RH / layer C (proxy); separate k_meth for
//       the two drivers (k_meth_ch4, k_meth_proxy).
//   MeHg(t+1) = Meq + (MeHg(t) - Meq) exp(-kd), Meq = P / kd, capped at the layer Hg.
#include <cmath>
#include <string>
#include <algorithm>
#include <json/json.h>
#include "cohortconst.h"
#include "layerconst.h"

namespace hgtr {

struct Params {
  int    tracer     = 0;      // hg_tracer: 0 off (default), 1 on
  double r_org      = 0.35;   // hg_r_org: initial RHgC of moss/fibric/humic layers (T0 Toolik-region organic median)
  double r_min      = 1.47;   // hg_r_min: initial RHgC of mineral layers (T0 synthesis mineral median)
  double r_leaf     = 0.07;   // hg_r_leaf: RHgC of new leaf tissue (NPP), ~30 ng/g dry mass
  double r_stem     = 0.02;   // hg_r_stem: RHgC of new stem tissue (wood ~5-10 ng/g)
  double r_root     = 0.07;   // hg_r_root: RHgC of new root tissue
  double f_evade    = 0.16;   // hg_f_evade: fraction of decay-released Hg that evades (Schaefer 2020)
  double f_export   = 0.015;  // hg_f_export: fraction of decay-released Hg exported (Schaefer 2020)
  double dep_soil   = 0.0;    // hg_dep_soil: direct Hg(II) deposition to the top soil layer (ug/m2/yr)
  int    mehg       = 0;      // hg_mehg: 0 off (default), 1 MeHg index
  int    mehg_driver = 0;     // hg_mehg_driver: 0 auto (CH4 production if the CH4 module is on), 1 proxy always
  // k_meth: calibrated (T1, 2026-11-23) to Thompson et al. 2025 growing-season (Jun-Sep) active-layer %MeHg in
  // steady state: CH4 driver -> 2.0 % (fen) at Toolik CMT06 wet sedge with the CH4 module on (k 30 gave 0.0189 %);
  // proxy driver -> 0.1 % (permafrost plateau) at HPC S21H (k 0.6 gave 0.00078 %). %MeHg is linear in k_meth.
  double k_meth_ch4   = 3200.0; // hg_k_meth_ch4 (-): MeHg production per unit Hg x CH4-C turnover fraction
  double k_meth_proxy = 77.0;   // hg_k_meth_proxy (-): same, per saturated-fraction RH C turnover fraction
  double k_demeth   = 3.0;    // hg_k_demeth (1/month at t_ref): ~0.1 /day, peat incubations
  double q10_meth   = 1.0;    // hg_q10_meth: extra T response of methylation (the C turnover already has one)
  double q10_demeth = 2.0;    // hg_q10_demeth
  double t_ref      = 10.0;   // hg_t_ref (C)
};

inline bool num_(const Json::Value& v) { return !v.isBool() && v.isNumeric() && std::isfinite(v.asDouble()); }

// Parse model_settings; returns "" or the FATAL message. Absent keys = defaults; unknown hg_* keys are errors.
inline std::string parse(const Json::Value& ms, Params& p) {
  p = Params();
  if (!ms.isObject()) { return ""; }
  struct D { const char* k; double* v; double lo, hi; };
  D ds[] = {{"hg_r_org", &p.r_org, 0, 1e3}, {"hg_r_min", &p.r_min, 0, 1e3}, {"hg_r_leaf", &p.r_leaf, 0, 1e3},
            {"hg_r_stem", &p.r_stem, 0, 1e3}, {"hg_r_root", &p.r_root, 0, 1e3}, {"hg_f_evade", &p.f_evade, 0, 1},
            {"hg_f_export", &p.f_export, 0, 1}, {"hg_dep_soil", &p.dep_soil, 0, 1e4},
            {"hg_k_meth_ch4", &p.k_meth_ch4, 0, 1e9}, {"hg_k_meth_proxy", &p.k_meth_proxy, 0, 1e9},
            {"hg_k_demeth", &p.k_demeth, 1e-9, 1e6}, {"hg_q10_meth", &p.q10_meth, 1e-3, 100},
            {"hg_q10_demeth", &p.q10_demeth, 1e-3, 100}, {"hg_t_ref", &p.t_ref, -50, 50}};
  struct I { const char* k; int* v; int lo, hi; };
  I is[] = {{"hg_tracer", &p.tracer, 0, 1}, {"hg_mehg", &p.mehg, 0, 1}, {"hg_mehg_driver", &p.mehg_driver, 0, 1}};
  for (const std::string& key : ms.getMemberNames()) {
    if (key.rfind("hg_", 0) != 0) { continue; }
    bool known = false;
    for (D& d : ds) {
      if (key == d.k) {
        known = true; const Json::Value& v = ms[key];
        if (!num_(v) || v.asDouble() < d.lo || v.asDouble() > d.hi) {
          return "model_settings." + key + " must be a number in [" + std::to_string(d.lo) + ", " + std::to_string(d.hi) + "]";
        }
        *d.v = v.asDouble();
      }
    }
    for (I& d : is) {
      if (key == d.k) {
        known = true; const Json::Value& v = ms[key];
        if (!num_(v) || v.asDouble() != std::floor(v.asDouble()) || v.asInt() < d.lo || v.asInt() > d.hi) {
          return "model_settings." + key + " must be an integer in [" + std::to_string(d.lo) + ", " + std::to_string(d.hi) + "]";
        }
        *d.v = v.asInt();
      }
    }
    if (!known) { return "model_settings." + key + ": unknown hg_* setting"; }
  }
  if (p.f_evade + p.f_export > 1.0) { return "model_settings.hg_f_evade + hg_f_export must be <= 1"; }
  if (p.mehg == 1 && p.tracer != 1) { return "model_settings.hg_mehg = 1 needs hg_tracer = 1"; }
  return "";
}

const int NPOOL = NUM_PFT * NUM_PFT_PART + NUM_PFT + 1;   // veg parts, standing dead, woody debris

// C snapshot of the column (g C/m2).
struct Snap {
  double veg[NUM_PFT][NUM_PFT_PART] = {};
  double dead[NUM_PFT] = {};
  double wd = 0.0;
  int nsl = 0;
  double soil[MAX_SOI_LAY] = {};
  double dz[MAX_SOI_LAY] = {};
  int type[MAX_SOI_LAY] = {};
};

// External C fluxes of one sub-step (g C/m2 over the sub-step).
struct Flux {
  double npp[NUM_PFT][NUM_PFT_PART] = {};   // net uptake (may be < 0: C loss without Hg)
  double decomp[MAX_SOI_LAY] = {};          // soil decomposition C (RH + CH4 production), by layer
  double doc[MAX_SOI_LAY] = {};             // DOC C exported, by layer
  double decomp_wd = 0.0;                   // woody debris RH
  double dep_soil = 0.0;                    // direct Hg deposition to the top soil layer (ug/m2), not C
  double burn = 0.0;                        // fire C emission (aggregate, from the donor pot)
  double lateral = 0.0;                     // RTS lateral C export (aggregate, from the donor pot)
  double pot_decomp = 0.0;                  // RTS respiration (aggregate decomposition from the pot)
};

struct Acc {
  double uptake = 0, dep = 0, evasion = 0, export_ = 0, fire = 0, unexpl_closs = 0, unexpl_cgain = 0;
  void clear() { *this = Acc(); }
  void add(const Acc& o) {
    uptake += o.uptake; dep += o.dep; evasion += o.evasion; export_ += o.export_; fire += o.fire;
    unexpl_closs += o.unexpl_closs; unexpl_cgain += o.unexpl_cgain;
  }
  double net() const { return uptake + dep - evasion - export_ - fire; }
};

struct State {
  bool init = false;
  double veg[NUM_PFT][NUM_PFT_PART] = {};
  double dead[NUM_PFT] = {};
  double wd = 0.0;
  int nsl = 0;
  double soil[MAX_SOI_LAY] = {};
  double mehg[MAX_SOI_LAY] = {};
  Acc m, y;               // this month / this year (outputs)
  Acc cum;                // since initialization (conservation diagnostic)
  double h_init = 0.0;    // total Hg at initialization
  double mehg_prod_m = 0.0, mehg_prod_y = 0.0;   // MeHg gross production (ug/m2) month / year
  double vegtot() const { double s = 0; for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) s += veg[i][j]; return s; }
  double deadtot() const { double s = wd; for (int i = 0; i < NUM_PFT; i++) s += dead[i]; return s; }
  double soiltot() const { double s = 0; for (int i = 0; i < MAX_SOI_LAY; i++) s += soil[i]; return s; }
  double mehgtot() const { double s = 0; for (int i = 0; i < MAX_SOI_LAY; i++) s += mehg[i]; return s; }
  double total() const { return vegtot() + deadtot() + soiltot(); }
  // conservation residual since initialization (ug/m2); roundoff only
  double balance_error() const { return total() - h_init - cum.net(); }
};

inline bool is_mineral(int type) { return type >= 3; }   // CohortData m_soil.type: 0 moss, 1 shallow, 2 deep, 3 mineral

// Initialize Hg from C x RHgC (soil by layer type; vegetation, standing dead and woody debris at the
// tissue ratios: leaf ratio for standing dead and woody debris).
inline void initialize(const Params& p, State& s, const Snap& c) {
  const double rv[NUM_PFT_PART] = {p.r_leaf, p.r_stem, p.r_root};
  for (int i = 0; i < NUM_PFT; i++) {
    for (int j = 0; j < NUM_PFT_PART; j++) { s.veg[i][j] = std::fmax(0.0, c.veg[i][j]) * rv[j]; }
    s.dead[i] = std::fmax(0.0, c.dead[i]) * p.r_leaf;
  }
  s.wd = std::fmax(0.0, c.wd) * p.r_leaf;
  s.nsl = c.nsl;
  for (int l = 0; l < MAX_SOI_LAY; l++) {
    s.soil[l] = (l < c.nsl) ? std::fmax(0.0, c.soil[l]) * (is_mineral(c.type[l]) ? p.r_min : p.r_org) : 0.0;
    s.mehg[l] = 0.0;
  }
  s.init = true;
  s.cum.clear();
  s.h_init = s.total();
}

inline double ratio_(double h, double c) { return (c > 1.0e-12 && h > 0.0) ? h / c : 0.0; }

inline bool same_geometry(const Snap& a, const Snap& b) {
  if (a.nsl != b.nsl) { return false; }
  for (int l = 0; l < a.nsl; l++) { if (a.dz[l] != b.dz[l]) { return false; } }
  return true;
}

// One sub-step. a, b: C before / after; f: the sub-step's external C fluxes. Updates s (and s.m, s.cum).
inline void step(const Params& p, State& s, const Snap& a, const Snap& b, const Flux& f) {
  if (!s.init) { return; }
  Acc acc;
  const double fr = p.f_evade + p.f_export;
  const double rv[NUM_PFT_PART] = {p.r_leaf, p.r_stem, p.r_root};
  const int na = std::min(a.nsl, MAX_SOI_LAY), nb = std::min(b.nsl, MAX_SOI_LAY);
  const bool same = same_geometry(a, b);

  // ---- non-soil pools (flattened): veg parts, standing dead, woody debris
  double* hp[NPOOL]; double dint[NPOOL]; double r[NPOOL];
  int k = 0;
  for (int i = 0; i < NUM_PFT; i++) {
    for (int j = 0; j < NUM_PFT_PART; j++, k++) {
      hp[k] = &s.veg[i][j];
      const double in = f.npp[i][j];
      const double hin = std::fmax(0.0, in) * rv[j];
      *hp[k] += hin; acc.uptake += hin;
      r[k] = ratio_(*hp[k], a.veg[i][j] + std::fmax(0.0, in));
      dint[k] = b.veg[i][j] - a.veg[i][j] - in;
    }
  }
  for (int i = 0; i < NUM_PFT; i++, k++) { hp[k] = &s.dead[i]; r[k] = ratio_(s.dead[i], a.dead[i]); dint[k] = b.dead[i] - a.dead[i]; }
  {
    hp[k] = &s.wd; r[k] = ratio_(s.wd, a.wd);
    const double rel = std::fmin(s.wd, std::fmax(0.0, f.decomp_wd) * r[k]);
    acc.evasion += rel * p.f_evade; acc.export_ += rel * p.f_export; s.wd -= rel * fr;
    dint[k] = b.wd - a.wd + f.decomp_wd;
    k++;
  }

  // ---- soil
  double sd[MAX_SOI_LAY] = {};      // same geometry: dC_int by layer
  double rs[MAX_SOI_LAY] = {};      // Hg:C at the start of the step
  double pot_decomp = f.pot_decomp, lateral = f.lateral;
  double soilnet = 0.0;             // different geometry: net soil C change (after the external fluxes)
  double h0soil[MAX_SOI_LAY];
  for (int l = 0; l < MAX_SOI_LAY; l++) { h0soil[l] = s.soil[l]; rs[l] = (l < na) ? ratio_(s.soil[l], a.soil[l]) : 0.0; }
  if (same) {
    for (int l = 0; l < na; l++) {
      double rel = std::fmin(s.soil[l], std::fmax(0.0, f.decomp[l]) * rs[l]);
      acc.evasion += rel * p.f_evade; acc.export_ += rel * p.f_export; s.soil[l] -= rel * fr;
      const double ex = std::fmin(s.soil[l], std::fmax(0.0, f.doc[l]) * rs[l]);
      acc.export_ += ex; s.soil[l] -= ex;
      sd[l] = b.soil[l] - a.soil[l] + f.decomp[l] + f.doc[l];
    }
  } else {
    double c0 = 0.0, c1 = 0.0, dec = 0.0, doc = 0.0;
    for (int l = 0; l < na; l++) { c0 += a.soil[l]; dec += f.decomp[l]; doc += f.doc[l]; }
    for (int l = 0; l < nb; l++) { c1 += b.soil[l]; }
    pot_decomp += dec; lateral += doc;          // per-layer fluxes with a geometry change: from the pot
    soilnet = c1 - c0 + dec + doc;
  }
  if (f.dep_soil > 0.0 && nb > 0) { acc.dep += f.dep_soil; }   // added to the top layer after the remap

  // ---- donors / receivers
  double D = 0.0, R = 0.0;
  for (int i = 0; i < NPOOL; i++) { if (dint[i] < 0) D -= dint[i]; else R += dint[i]; }
  if (same) { for (int l = 0; l < na; l++) { if (sd[l] < 0) D -= sd[l]; else R += sd[l]; } }
  else { if (soilnet < 0) D -= soilnet; else R += soilnet; }
  const double S = std::fmax(0.0, f.burn) + std::fmax(0.0, lateral) + std::fmax(0.0, pot_decomp);
  const double g = (D > 0.0) ? std::fmin(1.0, (R + S) / D) : 0.0;
  const double cflow = D * g;
  acc.unexpl_closs = D - cflow;

  double P = 0.0; double give[NPOOL] = {}; double give_soil[MAX_SOI_LAY] = {}; double give_col = 0.0;
  for (int i = 0; i < NPOOL; i++) {
    if (dint[i] < 0) { give[i] = std::fmin(*hp[i], -dint[i] * r[i] * g); *hp[i] -= give[i]; P += give[i]; }
  }
  // MeHg follows Hg losses of a layer proportionally (same geometry); with a geometry change it is mapped.
  double mehg_new[MAX_SOI_LAY] = {};
  if (same) {
    for (int l = 0; l < na; l++) {
      if (sd[l] < 0) { give_soil[l] = std::fmin(s.soil[l], -sd[l] * rs[l] * g); s.soil[l] -= give_soil[l]; P += give_soil[l]; }
      mehg_new[l] = (h0soil[l] > 0.0) ? s.mehg[l] * (s.soil[l] / h0soil[l]) : 0.0;
    }
  } else {
    // old column (top -> bottom), take the net loss from the top
    double cc[MAX_SOI_LAY], hh[MAX_SOI_LAY], mm[MAX_SOI_LAY];
    for (int l = 0; l < na; l++) { cc[l] = std::fmax(0.0, a.soil[l]); hh[l] = s.soil[l]; mm[l] = s.mehg[l]; }
    if (soilnet < 0) {
      double t = -soilnet * g;
      for (int l = 0; l < na && t > 0.0; l++) {
        const double x = std::fmin(cc[l], t);
        const double fx = (cc[l] > 0.0) ? x / cc[l] : 0.0;
        const double h = hh[l] * fx, m = mm[l] * fx;
        cc[l] -= x; hh[l] -= h; mm[l] -= m; t -= x; give_col += h; P += h;
      }
      if (t > 0.0) { acc.unexpl_closs += t; }
    }
    // bottom-up map onto the new layers
    double newh[MAX_SOI_LAY] = {}, newm[MAX_SOI_LAY] = {}, excess[MAX_SOI_LAY] = {};
    int i = na - 1;
    for (int j = nb - 1; j >= 0; j--) {
      double need = std::fmax(0.0, b.soil[j]);
      while (need > 0.0 && i >= 0) {
        if (cc[i] <= 0.0) { newh[j] += hh[i]; newm[j] += mm[i]; hh[i] = 0; mm[i] = 0; cc[i] = 0; i--; continue; }
        const double x = std::fmin(cc[i], need);
        const double fx = x / cc[i];
        const double h = (fx >= 1.0) ? hh[i] : hh[i] * fx, m = (fx >= 1.0) ? mm[i] : mm[i] * fx;
        newh[j] += h; newm[j] += m; hh[i] -= h; mm[i] -= m; cc[i] -= x; need -= x;
        if (fx >= 1.0) { cc[i] = 0; i--; }
      }
      excess[j] = need;
    }
    // leftover old material (C loss with no sink): Hg stays, in the top new layer
    double lh = 0.0, lm = 0.0;
    for (int l = 0; l <= i && l < na; l++) { lh += hh[l]; lm += mm[l]; }
    if (nb > 0) { newh[0] += lh; newm[0] += lm; } else { P += lh; give_col += lh; }
    for (int l = 0; l < MAX_SOI_LAY; l++) { s.soil[l] = (l < nb) ? newh[l] : 0.0; mehg_new[l] = (l < nb) ? newm[l] : 0.0; }
    // receivers: the top excess
    R = 0.0;
    for (int ii = 0; ii < NPOOL; ii++) { if (dint[ii] > 0) R += dint[ii]; }
    for (int j = 0; j < nb; j++) { sd[j] = excess[j]; R += excess[j]; }
  }
  for (int l = 0; l < MAX_SOI_LAY; l++) { s.mehg[l] = mehg_new[l]; }

  // ---- sinks from the pot
  if (P > 0.0 && cflow > 0.0) {
    const double sc = (S > cflow) ? cflow / S : 1.0;
    const double hb = P * std::fmax(0.0, f.burn) * sc / cflow;
    const double hl = P * std::fmax(0.0, lateral) * sc / cflow;
    const double hd = P * std::fmax(0.0, pot_decomp) * sc / cflow;
    acc.fire += hb; acc.export_ += hl + hd * p.f_export; acc.evasion += hd * p.f_evade;
    P -= hb + hl + hd * fr;
    if (P < 0.0) { P = 0.0; }
  }
  // ---- receivers
  if (P > 0.0) {
    if (R > 0.0) {
      double left = P;
      int last = -1;
      for (int i = 0; i < NPOOL; i++) { if (dint[i] > 0) last = i; }
      int lastl = -1;
      const int ns = same ? na : nb;
      for (int l = 0; l < ns; l++) { if (sd[l] > 0) lastl = l; }
      for (int i = 0; i < NPOOL; i++) {
        if (dint[i] > 0) { const double h = P * dint[i] / R; *hp[i] += h; left -= h; }
      }
      for (int l = 0; l < ns; l++) {
        if (sd[l] > 0) { const double h = P * sd[l] / R; s.soil[l] += h; left -= h; }
      }
      // roundoff remainder to the last receiver (exact conservation)
      if (lastl >= 0) { s.soil[lastl] += left; } else if (last >= 0) { *hp[last] += left; }
    } else {
      // nothing received C: the remaining pot goes back to the donors
      double tg = 0.0; for (int i = 0; i < NPOOL; i++) tg += give[i];
      for (int l = 0; l < MAX_SOI_LAY; l++) tg += give_soil[l];
      tg += give_col;
      double left = P;
      if (tg > 0.0) {
        for (int i = 0; i < NPOOL; i++) { if (give[i] > 0) { const double h = P * give[i] / tg; *hp[i] += h; left -= h; } }
        for (int l = 0; l < MAX_SOI_LAY; l++) { if (give_soil[l] > 0) { const double h = P * give_soil[l] / tg; s.soil[l] += h; left -= h; } }
      }
      if (nb > 0) { s.soil[0] += left; } else { s.wd += left; }
    }
  }
  if (acc.dep > 0.0) { s.soil[0] += acc.dep; }
  acc.unexpl_cgain = std::fmax(0.0, R - (cflow - std::fmin(S, cflow)));
  s.nsl = nb;
  s.m.add(acc); s.y.add(acc); s.cum.add(acc);
}

// MeHg index for one month. a_turn[l]: anaerobic C turnover fraction of the month (-); ts[l]: layer T (C).
inline void mehg_month(const Params& p, State& s, const double* a_turn, const double* ts, int nsl, bool ch4_driver) {
  const double km = ch4_driver ? p.k_meth_ch4 : p.k_meth_proxy;
  double prod = 0.0;
  for (int l = 0; l < MAX_SOI_LAY; l++) {
    if (l >= nsl || !(s.soil[l] > 0.0)) { s.mehg[l] = 0.0; continue; }
    const double t = std::isfinite(ts[l]) ? ts[l] : p.t_ref;
    const double pr = km * s.soil[l] * std::fmax(0.0, a_turn[l]) * std::pow(p.q10_meth, (t - p.t_ref) / 10.0);
    const double kd = p.k_demeth * std::pow(p.q10_demeth, (t - p.t_ref) / 10.0);
    const double meq = pr / kd;
    s.mehg[l] = meq + (s.mehg[l] - meq) * std::exp(-kd);
    s.mehg[l] = std::fmin(std::fmax(s.mehg[l], 0.0), s.soil[l]);
    prod += pr;
  }
  s.mehg_prod_m += prod; s.mehg_prod_y += prod;
}

// Restart layout (version 1): [0] version, veg (NUM_PFT x 3), dead (NUM_PFT), wd, soil (MAX_SOI_LAY),
// mehg (MAX_SOI_LAY), h_init, cum uptake, dep, evasion, export, fire (diagnostic continuity).
const int NSTATE = 1 + NUM_PFT * NUM_PFT_PART + NUM_PFT + 1 + 2 * MAX_SOI_LAY + 6;
inline void pack(const State& s, double* v) {
  int k = 0; v[k++] = 1.0;
  for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) v[k++] = s.veg[i][j];
  for (int i = 0; i < NUM_PFT; i++) v[k++] = s.dead[i];
  v[k++] = s.wd;
  for (int l = 0; l < MAX_SOI_LAY; l++) v[k++] = s.soil[l];
  for (int l = 0; l < MAX_SOI_LAY; l++) v[k++] = s.mehg[l];
  v[k++] = s.h_init; v[k++] = s.cum.uptake; v[k++] = s.cum.dep; v[k++] = s.cum.evasion; v[k++] = s.cum.export_; v[k++] = s.cum.fire;
}
inline bool unpack(State& s, const double* v, int nsl) {
  if (!(v[0] == 1.0)) { return false; }
  int k = 1;
  for (int i = 0; i < NUM_PFT; i++) for (int j = 0; j < NUM_PFT_PART; j++) s.veg[i][j] = v[k++];
  for (int i = 0; i < NUM_PFT; i++) s.dead[i] = v[k++];
  s.wd = v[k++];
  for (int l = 0; l < MAX_SOI_LAY; l++) s.soil[l] = v[k++];
  for (int l = 0; l < MAX_SOI_LAY; l++) s.mehg[l] = v[k++];
  s.h_init = v[k++]; s.cum.clear(); s.cum.uptake = v[k++]; s.cum.dep = v[k++]; s.cum.evasion = v[k++]; s.cum.export_ = v[k++]; s.cum.fire = v[k++];
  s.nsl = nsl; s.init = true; s.m.clear(); s.y.clear();
  return true;
}

} // namespace hgtr
#endif /*HGTRACER_H_*/
