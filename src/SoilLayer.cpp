/*! \file
 *
 */
#include "../include/SoilLayer.h"

int SoilLayer::kersten_form = 0; // DRAFT (tvc-summer), set from model_settings.kersten_form
// Option 1 (v5h): see SoilLayer.h; set from model_settings in Cohort.
int    SoilLayer::org_tc_scheme = 0;
double SoilLayer::org_tc_kunf   = 0.6;
double SoilLayer::org_tc_kfrz   = 0.25;
double SoilLayer::org_tc_dry    = -1.0;
double SoilLayer::org_tc_thre_f = 0.275;
double SoilLayer::org_tc_thre_h = 0.62;
double SoilLayer::org_tc_depth  = 10.0;
// Option 3 (v5h): see SoilLayer.h
double SoilLayer::mineral_tcsolid    = -1.0;
double SoilLayer::mineral_vhcsolid   = -1.0;
int    SoilLayer::mineral_tcdry_form = 0;

#include "../include/TEMLogger.h"
extern src::severity_logger< severity_level > glg;

SoilLayer::SoilLayer() {
  BOOST_LOG_SEV(glg, debug) << "==> Creating a SoilLayer object...";
  isSoil = true;
};

SoilLayer::~SoilLayer() {
  BOOST_LOG_SEV(glg, debug) << "--> Deleting a SoilLayer object...";
};

double SoilLayer::getFrzVolHeatCapa() {
  double vhc = vhcsolid * (1-poro) + (liq+ice)/dz *SHCICE;
  return vhc;
};

double SoilLayer::getUnfVolHeatCapa() {
  double vhc= vhcsolid * (1-poro) + (liq+ice)/dz *SHCLIQ;
  return vhc;
};

//Yuan: unfrozen/frozen put together
double SoilLayer::getMixVolHeatCapa() {
  double vhc = vhcsolid * (1-poro) + liq/dz *SHCLIQ+ice/dz *SHCICE;
  return vhc;
};

// Option 1 (v5h): Cote & Konrad (2005) / CLASS organic and moss
// thermal conductivity. Returns a negative value if the option is off or the
// layer is not a moss/organic layer within org_tc_depth.
double SoilLayer::getOrgThermCond(const bool frozen_case) {
  if (org_tc_scheme != 1) { return -1.0; }
  if (!(isMoss || isOrganic)) { return -1.0; }
  if (z >= org_tc_depth) { return -1.0; }
  double thre = isHumic ? org_tc_thre_h : org_tc_thre_f;
  if (thre <= 0.0 || thre > poro) { thre = poro; }
  double sr = (getVolIce() + getVolLiq()) / thre;
  sr = fmin(fmax(sr, 0.0), 1.0);
  double kap = frozen_case ? org_tc_kfrz : org_tc_kunf;
  double ke = kap * sr / (1.0 + (kap - 1.0) * sr);
  double kdry = (org_tc_dry > 0.0) ? org_tc_dry : 0.30 * exp(-2.0 * poro);
  double ksat = frozen_case ? tcsatfrz : tcsatunf;
  return ke * ksat + (1.0 - ke) * kdry;
}

// get frozen layer thermal conductivity
double SoilLayer::getFrzThermCond() {
  if (org_tc_scheme == 1) { // Option 1 (v5h)
    double tco = getOrgThermCond(true);
    if (tco >= 0.0) { return tco; }
  }
  double tc;
  double vice = getVolIce();
  double vliq = getVolLiq();
  double s = (vice + vliq)/poro;
  s = min(s, 1.0);
  double ke= s; // for frozen case

  if(s < 1.e-7) {
    tc = tcdry;
  } else {
    tc = ke *tcsatfrz + (1-ke)*tcdry;
  }

  return tc;
};

// get unfrozen layer thermal conductivity
double SoilLayer::getUnfThermCond() {
  if (org_tc_scheme == 1) { // Option 1 (v5h)
    double tco = getOrgThermCond(false);
    if (tco >= 0.0) { return tco; }
  }
  double tc;
  double vice = getVolIce();
  double vliq = getVolLiq();
  double s = (vice + vliq)/poro;
  s = min(s, 1.0);
  double ke= (kersten_form == 1) ? log10(s) + 1 : log(s) +1; // for unfrozen case (DRAFT option: log10, Johansen fine / CLM)
  ke = fmax(ke, 0.);

  if(s < 1.e-7) {
    tc = tcdry;
  } else {
    tc = ke *tcsatunf + (1-ke)*tcdry;
  }

  if(poro>=0.9 || (poro>=0.8 &&solind ==1)) {
    tc = fmax(tc, tcmin);
  }

  return tc;
};

double SoilLayer::getMatricPotential() {
  double psi;
  double lf = 3.337e5 ;// latent heat of fusion J/kg
  double g =9.8;

  if(tem<0) {
    psi =1000. * lf/g *(tem/(tem+273.16));

    if (psi<-1.e8) {
      psi=-1.e8;
    }
  } else {
    double voliq = getEffVolLiq();
    double ws = fmax(0.01, voliq);
    ws = fmin(1.0, ws);
    psi = psisat * pow(ws, -bsw);

    if (psi<-1.e8) {
      psi=-1.e8;
    }
  }

  return psi;
};

// get albedo of visible radiation
double SoilLayer::getAlbedoVis() {
  double vis;
  double liq1 = getVolLiq();
  double ice1 = getVolIce();
  double delta = 0.11-0.4*(liq1+ice1);
  delta =fmax(0., delta);
  vis = albsatvis + delta;
  vis = min(vis, (double)albdryvis);
  return vis;
};

// get albedo of nir radiation
double SoilLayer::getAlbedoNir() {
  double nir;
  double wat = getVolLiq()+getVolIce();
  double delta = 0.11-0.4*wat;
  delta =fmax(0., delta);
  nir = albsatnir + delta;
  nir = fmin(nir, (double)albdrynir);
  return nir;
};

// derive properties from the assigned property
// called when porosity/thickness is changed
void SoilLayer::derivePhysicalProperty() {
  //hydraulic properties
  minliq = 0.05*poro * DENLIQ * dz;
  maxliq = poro * DENLIQ * dz;
  maxice = poro * DENICE * dz;
  //thermal properties
  tcsatunf= pow((double)tcsolid , (double)1- poro) * pow((double)TCLIQ, (double)poro);
  tcsatfrz= pow((double)tcsolid , (double)1- poro) * pow((double)TCICE, (double)poro);
  tcdry   = pow((double)tcsolid , (double)1- poro) * pow((double)TCAIR, (double)poro);
  tcmin   = tcdry;
  // Option 3 (v5h): Johansen (1975) dry conductivity for mineral
  if (isMineral && mineral_tcdry_form == 1) {
    tcdry = (0.135 * bulkden + 64.7) / (2700. - 0.947 * bulkden);
    tcmin = tcdry;
  }
};

