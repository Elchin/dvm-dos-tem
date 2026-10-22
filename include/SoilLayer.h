/*! \file
 *
 */
#ifndef SOILLAYER_H_
#define SOILLAYER_H_
#include "Layer.h"

#include <math.h>
#include <cmath>
#include <memory>
using namespace std;

class SoilLayer:public Layer {
public:

  SoilLayer();
  virtual ~SoilLayer();

  virtual double getFrzThermCond();//get frozen layer thermal conductivity
  virtual double getUnfThermCond();//get unfrozen layer thermal conductivity
  virtual double getFrzVolHeatCapa();//get frozen layer specific heat capcity
  virtual double getUnfVolHeatCapa();//get unfrozen layer specific heat capacity
  virtual double getMixVolHeatCapa(); //Yuan: all soil components

  double getAlbedoVis();// get albedo of visible radiation
  double getAlbedoNir();// get albedo of Nir radiation

  void derivePhysicalProperty();

  double getMatricPotential();
  double getHydraulicCond();

  // DRAFT (tvc-summer): Kersten number form for unfrozen soil.
  // 0 (default) = legacy Ke = ln(Sr) + 1 (Ke = 0 for Sr <= 1/e = 0.37);
  // 1 = Johansen (1975) fine-grained / CLM / Lawrence & Slater (2008) form
  //     Ke = log10(Sr) + 1 (Ke = 0 only for Sr <= 0.1).
  static int kersten_form;
  // Refs: Johansen 1975 (CRREL transl. 637), doi:10.21236/ADA044002; Lawrence & Slater 2008, doi:10.1007/s00382-007-0278-1; see docs/daily-mode-REFERENCES.md.

  // Option 1 (v5h, OFF by default): moisture-dependent thermal
  // conductivity of moss and organic (fibric/humic) layers following
  // Cote & Konrad (2005) as used in CLASS/CLASSIC (Verseghy 2017; Melton et
  // al. 2019; Elshamy et al. 2020), with Sr relative to the peat retention
  // capacity theta_re (Letts et al. 2000). With kappa = 1 and
  // theta_re = porosity this is the linear dry-saturated interpolation of the
  // JULES moss layer (Chadburn et al. 2015).
  //   k  = Ke*k_sat + (1-Ke)*k_dry,  Ke = kappa*Sr / (1 + (kappa-1)*Sr)
  //   Sr = min(1, (theta_liq + theta_ice) / theta_re)
  //   k_dry = org_tc_dry if > 0, else 0.30*exp(-2.0*porosity) (CLASS organic)
  //   k_sat = existing tcsatunf / tcsatfrz (geometric mean, unchanged)
  // Applies only to moss/organic layers whose top is shallower than
  // org_tc_depth (m). model_settings keys (defaults):
  //   org_tc_scheme 0 (off) | 1 (on); org_tc_kunf 0.6; org_tc_kfrz 0.25;
  //   org_tc_dry -1; org_tc_thre_f 0.275 (moss + fibric); org_tc_thre_h 0.62
  //   (humic); org_tc_depth 10.0
  static int    org_tc_scheme;
  static double org_tc_kunf, org_tc_kfrz, org_tc_dry, org_tc_thre_f,
                org_tc_thre_h, org_tc_depth;
                // Refs: Cote & Konrad 2005, doi:10.1139/t04-106; Letts et al. 2000, doi:10.1080/07055900.2000.9649643; Chadburn et al. 2015, doi:10.5194/gmd-8-1493-2015.
  double getOrgThermCond(const bool frozen_case); // < 0 if not applicable

  // Option 3 (v5h, OFF by default): mineral-layer thermal
  // properties as settable values instead of texture-derived/hard-coded ones.
  //   mineral_tcsolid  > 0: solid conductivity (W/m/K) replaces the
  //       (8.8*sand + 2.92*clay)/(sand+clay) texture formula (CLM4.5/5).
  //   mineral_vhcsolid > 0: solid volumetric heat capacity (J/m3/K) replaces
  //       (2.128*sand + 2.385*clay)/(sand+clay)*1e6 (CLM4.5/5).
  //   mineral_tcdry_form 1: dry conductivity from Johansen (1975)
  //       (0.135*rho_b + 64.7)/(2700 - 0.947*rho_b) (Farouki 1981; CLM;
  //       Lawrence & Slater 2008) instead of the legacy geometric mean
  //       tcsolid^(1-poro)*k_air^poro that overwrites it in
  //       derivePhysicalProperty.
  // Defaults (key absent): -1, -1, 0 = legacy behaviour.
  static double mineral_tcsolid, mineral_vhcsolid;
  static int    mineral_tcdry_form;
  // Refs: Johansen 1975, doi:10.21236/ADA044002; Farouki 1981 (CRREL Monograph 81-1), doi:10.21236/ADA111734; CLM5 Tech Note eqs. 2.6.80, 2.6.82.

};
#endif /*SOILLAYER_H_*/
