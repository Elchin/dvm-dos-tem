/*! \file Soil_Env_enthalpy.cpp
 * Option C (enthalpy/stageC): the daily snow/soil/rock thermal step with the implicit enthalpy
 * solver (EnthalpySolver), used when model_settings.thermal_solver = 1. Nothing here runs with
 * thermal_solver = 0 (legacy Stefan + Crank-Nicolson path, bit-identical).
 *
 * One cell per Layer (snow, soil, rock). Each day:
 *  1. book the enthalpy change since the last step (non-thermal modules, Cohort windows);
 *  2. on switch-on / after a restart: project onto the freezing curve at fixed T (booked: init);
 *  3. put every cell back on its curve keeping H (water moved by other modules; equil ~ 0);
 *  4. backward-Euler steps (thermal_substeps per day; dt halving on non-convergence) with the
 *     radiation/SEB melt energy as a heat source in the top snow cell, Robin top boundary
 *     (legacy conductance cn = 100 W m-2 K-1 to the surface temperature, no fixed 0 C rule while
 *     melting) and the geothermal flux at the base;
 *  5. scatter T/liq/ice; heat of fully melted snow cells goes to the cell below; pore-volume ice
 *     cap as legacy (expelled water, booked: clip);
 *  6. rebuild fronts, frozen, frozenfrac from the profile; then the legacy daily diagnostics
 *     (rooting-zone T, unfrozen column, fronts to ed, ALD/ALC).
 * Design: tools/enthalpy_bench/enthalpy_solver_design.md (sections 4-12).
 * Refs: Tubini et al. 2021, doi:10.5194/tc-15-2541-2021; Casulli & Zanolli 2010,
 *   doi:10.1137/100786320; Schaefer & Jafarov 2016, doi:10.5194/bg-13-1991-2016.
 */
#include "../include/Soil_Env.h"
#include "../include/Snow_Env.h"
#include "../include/MineralLayer.h"
#include "../include/SJCurve.h"
#include "../include/physicalconst.h"

#include <cmath>
#include <vector>
#include <iostream>

#include "../include/TEMLogger.h"
extern src::severity_logger< severity_level > glg;

double Soil_Env::enthLayerCs(Layer* l) {
  if (l == NULL || l->isSnow) { return 0.0; }
  // rock: ParentLayer shadows Layer::vhcsolid (Layer's copy stays MISSING), use the virtual getter
  if (l->isRock) { return l->getFrzVolHeatCapa() * l->dz; }
  return l->vhcsolid * (1.0 - l->poro) * l->dz;
}

double Soil_Env::enthLayerH(Layer* l) {
  if (l == NULL) { return 0.0; }
  const double T = l->tem;
  const double liq = fmax(l->liq, 0.0), ice = fmax(l->ice, 0.0);
  return enthLayerCs(l) * T + liq * SHCLIQ * T + ice * (SHCICE * T - LHFUS);
}

double Soil_Env::columnEnthalpyE() const {
  double h = 0.0;
  for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) { h += enthLayerH(l); }
  return h;
}

void Soil_Env::enthalpyCell(Layer* l, enthalpy::Cell& c) const {
  c = enthalpy::Cell();
  c.dz = fmax(l->dz, 1.0e-6);
  c.Cs = enthLayerCs(l);
  c.eps = 0.01;   // legacy zerodegc
  if (l->isRock) {
    c.kind = enthalpy::ROCK; c.W = 0.0; c.liqmin = 0.0;
    c.lam_u = c.lam_f = l->getUnfThermCond();
    return;
  }
  c.W = fmax(l->liq, 0.0) + fmax(l->ice, 0.0);
  if (l->isSnow) {
    c.kind = enthalpy::SNOW; c.liqmin = 0.0; c.sfcc = enthalpy::SFCC_FREE;
    c.lam_u = c.lam_f = l->getUnfThermCond();
    return;
  }
  c.kind = enthalpy::SOIL;
  c.liqmin = fmax(l->minliq, 0.0);
  c.sfcc = (thermal_sfcc == 1) ? enthalpy::SFCC_SJ2016 : enthalpy::SFCC_FREE;
  c.forg = 1.0; c.clay = c.silt = c.sand = 0.0;   // moss/peat: organic (as Soil_Bgc::rhFrozenEfr)
  if (l->isMineral) {
    MineralLayer* ml = dynamic_cast<MineralLayer*>(l);
    if (ml != NULL && ml->pctclay + ml->pctsilt + ml->pctsand > 0.0) {
      c.clay = ml->pctclay; c.silt = ml->pctsilt; c.sand = ml->pctsand;
      c.forg = sjcurve::orgFraction(l->rawc + l->soma + l->sompr + l->somcr, l->dz);
    }
  }
  c.lam_u = l->getUnfThermCond();
  c.lam_f = l->getFrzThermCond();
}

double Soil_Env::enthMeltSource(const double tsurface) const {
  if (ground->toplayer == NULL || !ground->toplayer->isSnow) { return 0.0; }
  double melt_mm = 0.0;   // kg m-2 day-1
  if (ground->seb_on) {
    melt_mm = ground->seb_melt;          // SEB surplus at Ts = 0 (Soil_Env::solveSEB)
  } else if (snowenvp != NULL) {
    melt_mm = snowenvp->radMeltPotential(tsurface);   // Snow_Env::getSublimation radiation melt
  }
  return fmax(melt_mm, 0.0) * LHFUS / SEC_IN_DAY;      // W m-2
}

namespace {
// one backward-Euler step with dt halving on non-convergence (depth <= 6)
bool enthStep(Soil_Env& se, const std::vector<enthalpy::Cell>& cells, std::vector<double>& T,
              double dt, const enthalpy::BC& bc, const std::vector<double>& S, int depth,
              enthalpy::Stats& acc, int& halvings) {
  std::vector<double> T0 = T;
  enthalpy::Stats st;
  const bool ok = (se.thermal_nl_solver == 1) ? se.enth_solver.stepKatz(cells, T, dt, bc, &S, st)
                                              : se.enth_solver.stepNCZ(cells, T, dt, bc, &S, st);
  acc.solves += st.solves; acc.outer += st.outer;
  if (ok || depth >= 6) {
    acc.gtop += st.gtop; acc.gbot += st.gbot; acc.src += st.src; acc.dH += st.dH;
    acc.res_max = std::max(acc.res_max, st.res_max);
    acc.converged = acc.converged && ok;
    return ok;
  }
  halvings++;
  T = T0;
  bool a = enthStep(se, cells, T, 0.5 * dt, bc, S, depth + 1, acc, halvings);
  bool b = enthStep(se, cells, T, 0.5 * dt, bc, S, depth + 1, acc, halvings);
  return a && b;
}
}

// Abort with a layer dump if the column state is not finite (never expected; guards the budget).
void Soil_Env::enthCheckFinite(const char* where) const {
  if (std::isfinite(columnEnthalpyE())) { return; }
  std::cerr << "FATAL: thermal_solver 1: non-finite column enthalpy (" << where << "); layers:" << std::endl;
  int k = 0;
  for (Layer* l = ground->toplayer; l != NULL; l = l->nextl, ++k) {
    std::cerr << "  " << k << (l->isSnow ? " snow" : (l->isRock ? " rock" : " soil")) << " z=" << l->z << " dz=" << l->dz
              << " tem=" << l->tem << " liq=" << l->liq << " ice=" << l->ice << " poro=" << l->poro
              << " vhcsolid=" << l->vhcsolid << " minliq=" << l->minliq << " Cs=" << enthLayerCs(l)
              << " H=" << enthLayerH(l) << std::endl;
  }
  exit(EXIT_FAILURE);
}

void Soil_Env::enthalpyThermalStep(const double tsurface, const double qmelt_Wm2, EnergyDiag& q) {
  enthCheckFinite("start of the thermal step");
  // 2. switch-on / restart: fixed-T projection onto the curve (booked as init)
  if (enth_init_pending) {
    const double hb = columnEnthalpyE();
    double dmax = 0.0;
    for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) {
      if (l->isRock) { l->liq = 0.0; l->ice = 0.0; continue; }
      enthalpy::Cell c; enthalpyCell(l, c);
      const double lq = enthalpy::liq(c, l->tem);
      dmax = fmax(dmax, fabs(lq - l->liq));
      l->liq = lq; l->ice = fmax(0.0, c.W - lq);
    }
    q.init = columnEnthalpyE() - hb;
    enth_init_maxdliq = dmax;
    enth_init_pending = false;
    BOOST_LOG_SEV(glg, info) << "thermal_solver 1: switch-on projection onto the freezing curve at fixed T,"
                             << " max |d liq| = " << dmax << " kg/m2, dH = " << q.init << " J/m2";
  }
  // 3. H-conserving equilibration (cells changed by other modules)
  const double he0 = columnEnthalpyE();
  for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) {
    if (l->isRock) {
      if (l->liq != 0.0 || l->ice != 0.0) { l->liq = 0.0; l->ice = 0.0; }   // rock holds no water
      continue;
    }
    enthalpy::Cell c; enthalpyCell(l, c);
    if (!enthalpy::onCurve(c, l->tem, fmax(l->liq, 0.0), 1.0e-12 * fmax(1.0, c.W))) {
      double lq, ic;
      const double h = enthalpy::Hstate(c, l->tem, fmax(l->liq, 0.0), fmax(l->ice, 0.0));
      l->tem = enthalpy::equilibrate(c, h, lq, ic);
      l->liq = lq; l->ice = ic;
    }
  }
  enthCheckFinite("after equilibration");
  const double h0 = columnEnthalpyE();
  q.equil = h0 - he0;
  q.h0 = h0;

  // 4. solve
  std::vector<Layer*> L;
  for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) { L.push_back(l); }
  const int n = (int)L.size();
  std::vector<enthalpy::Cell> cells(n);
  std::vector<double> T(n), S(n, 0.0);
  for (int i = 0; i < n; ++i) { enthalpyCell(L[i], cells[i]); T[i] = L[i]->tem; }
  if (n > 0 && L[0]->isSnow) { S[0] = qmelt_Wm2; }
  enthalpy::BC bc;
  bc.top_type = 1; bc.cn_top = 1.0e2; bc.Ts = tsurface;   // legacy virtual top layer conductance
  bc.q_bot = geothermal_flux_Wm2;
  enth_solver.opt.face_mode = thermal_face_mode;
  const int nsub = (thermal_substeps > 0) ? thermal_substeps : 1;
  const double dt = SEC_IN_DAY / nsub;
  enthalpy::Stats acc; acc.converged = true;
  int halvings = 0;
  bool ok = true;
  for (int k = 0; k < nsub; ++k) {
    // cells re-assembled each substep (lagged conductivity from the current T)
    ok = enthStep(*this, cells, T, dt, bc, S, 0, acc, halvings) && ok;
  }
  if (halvings > 0) { enth_fallbacks++; }
  if (!ok) {
    enth_failures++;
    BOOST_LOG_SEV(glg, warn) << "thermal_solver 1: enthalpy step not converged after dt halving (res_max "
                            << acc.res_max << " J/m2)";
  }
  q.solves = acc.solves; q.fallbacks = ok ? 0 : 1;   // 1 = step not converged even after halving

  // 5. scatter
  for (int i = 0; i < n; ++i) {
    Layer* l = L[i];
    l->tem = T[i];
    l->tcond = enthalpy::conductivity(cells[i], T[i]);
    l->pce_f = 0.0; l->pce_t = 0.0;
    if (l->isRock) { l->liq = 0.0; l->ice = 0.0; continue; }
    const double lq = enthalpy::liq(cells[i], T[i]);
    l->liq = lq; l->ice = fmax(0.0, cells[i].W - lq);
  }
  // heat of fully melted snow cells (no ice left, T > 0) goes to the cell below (conservative)
  for (int i = 0; i + 1 < n; ++i) {
    Layer* l = L[i];
    if (!l->isSnow || l->ice > 0.0 || !(l->tem > 0.0)) { continue; }
    const double hx = enthLayerH(l);   // = liq c_w T > 0
    l->tem = 0.0;
    Layer* m = L[i + 1];
    enthalpy::Cell cm; enthalpyCell(m, cm);
    const double hm = enthalpy::Hstate(cm, m->tem, fmax(m->liq, 0.0), fmax(m->ice, 0.0)) + hx;
    double lq = 0.0, ic = 0.0;
    m->tem = enthalpy::equilibrate(cm, hm, lq, ic);
    if (!m->isRock) { m->liq = lq; m->ice = ic; }
  }
  enthCheckFinite("after the solve");
  const double hs = columnEnthalpyE();
  q.gtop = acc.gtop; q.gbot = acc.gbot; q.src = acc.src;
  q.res = (hs - h0) - (q.gtop + q.gbot + q.src);
  // pore-volume ice cap (as Stefan::updateWaterAfterFront): excess is expelled water
  for (int i = 0; i < n; ++i) {
    Layer* l = L[i];
    if (!l->isSoil) { continue; }
    const double icebylwc = l->getVolLiq() * DENICE * l->dz;
    if (l->ice >= (l->maxice - icebylwc)) {
      const double icecap = fmax(0., l->maxice - icebylwc);
      if (l->ice > icecap) {
        ground->expelled_water_mm += l->ice - icecap;
        l->ice = icecap;
      }
    }
  }

  // 6. layer states and fronts from the profile
  std::vector<Layer*> soil;
  for (int i = 0; i < n; ++i) {
    Layer* l = L[i];
    if (l->isSnow) {
      l->frozen = (l->ice <= 0.0) ? -1 : ((l->liq > 0.0) ? 0 : 1);
      l->frozenfrac = (l->liq + l->ice > 0.0) ? l->ice / (l->liq + l->ice) : 0.0;
    } else if (l->isRock) {
      l->frozen = (l->tem > 0.0) ? -1 : 1;
      l->frozenfrac = (l->frozen == 1) ? 1.0 : 0.0;
    } else if (l->isSoil) {
      soil.push_back(l);
    }
  }
  enthalpy::FrontIn fi; enthalpy::FrontOut fo;
  fi.eps = 0.01; fi.maxfronts = MAX_NUM_FNT;
  fi.Tabove = (soil.size() > 0 && soil[0]->prevl != NULL) ? soil[0]->prevl->tem : tsurface;
  fi.Tbelow = (soil.size() > 0 && soil.back()->nextl != NULL) ? soil.back()->nextl->tem : 0.0;
  for (size_t k = 0; k < soil.size(); ++k) {
    Layer* l = soil[k];
    enthalpy::Cell c; enthalpyCell(l, c);
    fi.z.push_back(l->z); fi.dz.push_back(l->dz); fi.T.push_back(l->tem); fi.ice.push_back(l->ice);
    fi.W.push_back(c.W);
    fi.ice_eps.push_back(fmax(0.0, c.W - enthalpy::liq(c, -c.eps)));
    fi.latent.push_back(enthalpy::hasLatent(c) ? 1 : 0);
  }
  enthalpy::reconstructFronts(fi, fo);
  ground->frontsz.clear(); ground->frontstype.clear();
  for (size_t k = 0; k < fo.fz.size(); ++k) {
    ground->frontsz.push_back(fo.fz[k]); ground->frontstype.push_back(fo.ftype[k]);
  }
  for (size_t k = 0; k < soil.size(); ++k) {
    soil[k]->frozen = fo.frozen[k]; soil[k]->frozenfrac = fo.frozenfrac[k];
  }
  ground->setFstLstFrontLayers();
  ground->checkWaterValidity();
  const double h1 = columnEnthalpyE();
  q.clip = h1 - hs;
  q.h1 = h1;
  q.nfront = (int)ground->frontstype.size();
  int ns = 0; for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) { if (l->isSnow) ns++; }
  q.nsnow = ns;
}

void Soil_Env::updateDailyGroundT_enthalpy(const double tsurface) {
  EnergyDiag q;
  q.used = energy_diag; q.enth = 1; q.tsurf = tsurface;
  // 1. between-step bookings
  const double hraw = columnEnthalpyE();
  q.bsnow = enth_bk_snow; q.bwater = enth_bk_water; q.blayer = enth_bk_layer; q.bfire = enth_bk_fire;
  if (enth_have_last) {
    q.unattr = (hraw - enth_h_last) - (q.bsnow + q.bwater + q.blayer + q.bfire);
  }
  enth_bk_snow = enth_bk_water = enth_bk_layer = enth_bk_fire = 0.0;
  q.rsnow = enth_rm_snow; q.rwater = enth_rm_water; q.wlayer = enth_wk_layer; q.remapfb = enth_remap_fallbacks;
  enth_rm_snow = enth_rm_water = enth_wk_layer = 0.0;
  stefan.initpce();
  const double qmelt = enthMeltSource(tsurface);
  enthalpyThermalStep(tsurface, qmelt, q);
  q.resday = q.res + q.equil + q.unattr;
  enth_h_last = q.h1; enth_have_last = true;
  if (energy_diag == 1) { energydiag = q; }

  ground->updateWholeFrozenStatus(); //for the whole column
  // at end of each day, 'ed' is updated for thermal properties (legacy functions)
  updateDailySoilThermal4Growth(ground->fstsoill, tsurface);
  updateLayerStateAfterThermal(ground->fstsoill, ground->lstsoill, ground->botlayer);
  retrieveDailyFronts();
}

// ---------------------------------------------------------------------------------------------
// Stage C2: energy-conserving remap of the snow and soil-water windows (design section 10).
// The modules move mass at fixed T; here every moved kg carries its enthalpy, mass that enters or
// leaves the column carries a physically defined enthalpy, and the cell temperatures are set from
// the resulting H at the module's liq/ice (H = C T - ice L is linear in T at fixed liq/ice, so this
// is exact; the next thermal step equilibrates onto the curve conserving H).
static const double ENTH_DIRECT_DEPTH = 0.10;   // m, soil depth that supplies the latent heat of direct snow melt

void Soil_Env::enthSnap(EnthSnap& s) const {
  s = EnthSnap();
  for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) {
    s.L.push_back(l); s.snow.push_back(l->isSnow ? 1 : 0); s.T.push_back(l->tem);
    s.liq.push_back(fmax(l->liq, 0.0)); s.ice.push_back(fmax(l->ice, 0.0)); s.H.push_back(enthLayerH(l));
  }
  s.hcol = columnEnthalpyE();
  s.topsnow = (ground->toplayer != NULL && ground->toplayer->isSnow);
  s.ttop = (ground->toplayer != NULL) ? ground->toplayer->tem : 0.0;
}

void Soil_Env::enthSetT(Layer* l, double h) {
  const double liq = fmax(l->liq, 0.0), ice = fmax(l->ice, 0.0);
  const double C = enthLayerCs(l) + liq * SHCLIQ + ice * SHCICE;
  if (C > 0.0) { l->tem = (h + ice * LHFUS) / C; }
}

double Soil_Env::columnWaterE() const {
  double w = 0.0;
  for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) {
    if (!l->isRock) { w += fmax(l->liq, 0.0) + fmax(l->ice, 0.0); }
  }
  return w;
}

void Soil_Env::enthRemapWaterCore(const EnthSnap& s, double infil_kg, double h_in) {
  const double hfix = columnEnthalpyE();
  std::vector<Layer*> now;
  for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) { now.push_back(l); }
  if (now != s.L) {   // layer structure changed inside the window: keep the fixed-T result, book it
    enth_remap_fallbacks++; enth_bk_water += hfix - s.hcol; return;
  }
  const int n = (int)now.size();
  std::vector<double> lossE(n, 0.0), gain(n, 0.0);
  double I = fmax(infil_kg, 0.0), P = I, E = I * h_in, G = 0.0;
  for (int i = 0; i < n; ++i) {
    Layer* l = now[i];
    if (l->isRock) { continue; }
    const double dl = fmax(l->liq, 0.0) - s.liq[i], di = fmax(l->ice, 0.0) - s.ice[i];
    if (dl < 0.0) { lossE[i] += -dl * SHCLIQ * s.T[i]; P += -dl; }
    if (di < 0.0) { lossE[i] += -di * (SHCICE * s.T[i] - LHFUS); P += -di; }
    gain[i] = fmax(dl, 0.0) + fmax(di, 0.0);
    G += gain[i]; E += lossE[i];
  }
  if (G > P) { const double x = G - P; P += x; E += x * h_in; I += x; }  // unreported inflow
  const double hp = (P > 0.0) ? E / P : h_in;     // mixed specific enthalpy of the moved water
  const double out = P - G;                        // drainage + runoff + ET leaving the column
  for (int i = 0; i < n; ++i) {
    if (lossE[i] == 0.0 && gain[i] == 0.0) { continue; }
    enthSetT(now[i], s.H[i] - lossE[i] + gain[i] * hp);
  }
  const double phys = I * h_in - out * hp;
  const double hnew = columnEnthalpyE();
  enth_bk_water += hnew - s.hcol;            // = phys up to round-off
  enth_rm_water += hnew - hfix;
  (void)phys;
}

void Soil_Env::enthRemapWater(const EnthSnap& s) {
  const double rnth = fmax((ed->d_v2g.rthfl + ed->d_v2g.rdrip) + (1.0 - cd->m_vegd.fpcsum) * ed->d_a2l.rnfl, 0.0);
  const double melt = fmax(ed->d_snw2soi.melt, 0.0);
  const double frain = (rnth + melt > 0.0) ? rnth / (rnth + melt) : 0.0;
  // rain at max(Tair, 0 C); snow melt and ponded water at 0 C (liquid reference: h = 0)
  const double h_in = SHCLIQ * fmax(ed->d_atms.ta, 0.0) * frain;
  enthRemapWaterCore(s, ed->d_soi2l.qinfl, h_in);
}

void Soil_Env::enthRemapSnowCore(const EnthSnap& s, double tsurface, double sublim_kg, double melt_kg) {
  const double hfix = columnEnthalpyE();
  double Mb = 0.0, Ma = 0.0;
  for (size_t i = 0; i < s.L.size(); ++i) { if (s.snow[i]) { Mb += s.liq[i] + s.ice[i]; } }
  for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) {
    if (l->isSnow) { Ma += fmax(l->liq, 0.0) + fmax(l->ice, 0.0); }
  }
  const double h_sf = SHCICE * fmin(tsurface, 0.0) - LHFUS;   // new snow: ice at min(Ts, 0 C)
  const bool direct = (!s.topsnow && s.ttop > 0.01);          // Snow_Env: snowfall melts on warm ground
  const double subl = fmax(sublim_kg, 0.0), melt = fmax(melt_kg, 0.0);
  double phys;
  if (direct) {
    phys = melt * h_sf;   // ice in, melt water out at 0 C: the latent heat comes from the ground
  } else {
    const double msf = (Ma - Mb) + subl + melt;                // snow that entered the cells
    const double h_sub = SHCICE * fmin(s.topsnow ? s.ttop : tsurface, 0.0) - LHFUS;
    phys = msf * h_sf - subl * h_sub;                          // melt water leaves at 0 C (h = 0)
  }
  const double delta = (s.hcol + phys) - hfix;
  double C = 0.0;
  for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) {
    if (l->isSnow) { C += fmax(l->liq, 0.0) * SHCLIQ + fmax(l->ice, 0.0) * SHCICE; }
  }
  if (C > 0.0) {
    const double dT = delta / C;
    for (Layer* l = ground->toplayer; l != NULL; l = l->nextl) { if (l->isSnow) { l->tem += dT; } }
  } else {
    // no snow cells (e.g. snowfall that melted on warm ground): the correction (mostly the latent heat
    // of the melted snowfall) is taken from the top ENTH_DIRECT_DEPTH of the soil, shared by heat
    // capacity (equal dT), not from the first cell alone: a thin dry moss cell has almost no heat
    // capacity and would get an unphysical temperature (seen: -2850 C) before the next thermal step.
    std::vector<Layer*> cells; std::vector<double> cap; double ctot = 0.0, ztop = 0.0;
    for (Layer* l = ground->fstsoill; l != NULL && !l->isSnow; l = l->nextl) {
      if (!cells.empty() && ztop >= ENTH_DIRECT_DEPTH) { break; }
      const double c = enthLayerCs(l) + fmax(l->liq, 0.0) * SHCLIQ + fmax(l->ice, 0.0) * SHCICE;
      cells.push_back(l); cap.push_back(c); ctot += c; ztop += fmax(l->dz, 0.0);
    }
    if (ctot > 0.0) {
      for (size_t k = 0; k < cells.size(); ++k) { enthSetT(cells[k], enthLayerH(cells[k]) + delta * cap[k] / ctot); }
    } else {
      Layer* t = (ground->fstsoill != NULL) ? ground->fstsoill : ground->toplayer;
      if (t != NULL) { enthSetT(t, enthLayerH(t) + delta); }
    }
  }
  const double hnew = columnEnthalpyE();
  enth_bk_snow += hnew - s.hcol;             // = phys up to round-off
  enth_rm_snow += hnew - hfix;
}

void Soil_Env::enthRemapSnow(const EnthSnap& s, const double tdrv) {
  const double tsurface = ground->seb_on ? ground->seb_ts : tdrv * ed->d_soid.nfactor;   // as Snow_Env
  enthRemapSnowCore(s, tsurface, ed->d_snw2a.sublim, ed->d_snw2soi.melt);
}
