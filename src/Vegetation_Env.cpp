/*
 * Vegetation_Env.cpp
 *
 * Purpose: Calculating radiation, water processes of canopy
 *
 * History:
 *   June 28, 2011, by F.-M. Yuan:
 *     (1) Recoding based on DOS-TEM's code;
 *
 * Important:
 *   (1) Parameters are read from 'CohortLookup.cpp', and set to
 *       'envpar' (struct:: vegpar_env)
 *
 *   (2) No calibrated Parameters at this stage
 *
 *   (3) The calculation is for ONE PFT only, so when calling it
 *       must be set ONE important index:
 *          pfttype
 *
 *   (4) THREE (3) data pointers must be initialized by calling
 *       corresponding 'set...' methods
 *          chtlu, ed, fd
 *
 */

#include "../include/Vegetation_Env.h"

Vegetation_Env::Vegetation_Env() {
  ipft = MISSING_I;
  // v5l ET-split switches: default 0 = legacy (see Vegetation_Env.h)
  fpc_once_fix = 0;
  nv_snow_bury_fix = 0;
  nv_snow_bury_depth = 0.0;
  rain_int_lai_fix = 0;
  nv_rain_surface = 0;
  nonvasc = 0;
  snow_depth = 0.0;
};

Vegetation_Env::~Vegetation_Env() {
};

//set the parameters from 'CohortLookup'
void Vegetation_Env::initializeParameter() {
  envpar.albvisnir = chtlu->albvisnir[ipft];
  envpar.er        = chtlu->er[ipft];
  envpar.ircoef    = chtlu->ircoef[ipft];
  envpar.iscoef    = chtlu->iscoef[ipft];
  envpar.gl_bl     = chtlu->gl_bl[ipft];
  envpar.gl_c      = chtlu->gl_c[ipft];
  envpar.vpd_close = chtlu->vpd_close[ipft];
  envpar.vpd_open  = chtlu->vpd_open[ipft];
  envpar.glmax     = chtlu->glmax[ipft];
  envpar.ppfd50    = chtlu->ppfd50[ipft];
};

// set the intial monthly LAI from 'CohortLookup'
void Vegetation_Env::initializeState() {
  ed->d_vegs.rwater = chtlu->initvegwater[ipft];
  ed->d_vegs.snow   = chtlu->initvegsnow[ipft];
};

void Vegetation_Env::set_state_from_restartdata(const RestartData & rdata) {
  ed->d_vegs.rwater = rdata.vegwater[ipft];
  ed->d_vegs.snow   = rdata.vegsnow[ipft];
};

//solar radiation (unit: W/m2) on canopy and its energy balance
void Vegetation_Env::updateRadiation(double leaf_area_index, double foliar_projected_cover) {

  double EPAR = 4.55; // average energy for PAR photon (umol/J)
  
  // solar radiation and its energy balance
  ed->d_v2a.swrefl = (ed->d_a2l.nirr * foliar_projected_cover) * envpar.albvisnir; //unit W/m2 (FPC adjusted)
  ed->d_a2v.swdown = (ed->d_a2l.nirr * foliar_projected_cover) - ed->d_v2a.swrefl;
  ed->d_v2g.swthfl  = ed->d_a2v.swdown * exp(-envpar.er * leaf_area_index);
  ed->d_a2v.swinter = ed->d_a2v.swdown - ed->d_v2g.swthfl;
  
  // PAR and its absorption by canopy
  double par = ed->d_a2l.par * foliar_projected_cover;   // fpc adjusted here
  ed->d_a2v.pardown    = par * (1.0 - envpar.albvisnir);
  
  // absorbed PAR: W/m2
  ed->d_a2v.parabsorb  = ed->d_a2v.pardown * (1.0 - exp(-envpar.er * leaf_area_index));
  
  double ppfd50  = envpar.ppfd50 ;//mumol/m2s, ppfd for 0.5 stomatal closure
  // Refs (S1): double FPC scaling introduced by commit 28167031; cover applied once as in LPJ (Gerten et al. 2004, doi:10.1016/j.jhydrol.2003.09.029),
  //   LPJmL4 (Schaphoff et al. 2018, doi:10.5194/gmd-11-1343-2018) and CLM5 Tech Note eq. 2.7.68; FPC = 1-exp(-k LAI), Monsi & Saeki 1953 (transl. doi:10.1093/aob/mci052).
  if (fpc_once_fix == 1 && foliar_projected_cover > 0.) {
    // S1 (v5l): stomatal light factor from the PAR per unit covered area
    double ppfd_c = ed->d_a2v.pardown / foliar_projected_cover * EPAR;
    ed->d_vegd.m_ppfd = ppfd_c / (ppfd_c + ppfd50);
  } else {
    ed->d_vegd.m_ppfd = ed->d_a2v.pardown * EPAR / (ed->d_a2v.pardown * EPAR + ppfd50);
  }
}

//VEGETATION DAILY WATER BALANCE CALCULATION
void Vegetation_Env::updateWaterBalance(const double & daylhr, double leaf_area_index, double foliar_projected_cover) {

  double daylsec = daylhr * 3600.0; // from hour to sec
  double EPAR = 4.55 ;              //an average energy for PAR photon (umol/J)

  // variables calculated in 'atmosphere.cpp'
  double vpd  = ed->d_atmd.vpd;
  double ta   = ed->d_atms.ta;
  double snfl = ed->d_a2l.snfl;  // this is the total atm to land
  double rnfl = ed->d_a2l.rnfl;
  ed->d_a2v.snfl = snfl * foliar_projected_cover;   //note: foliar_projected_cover adjusted here
  ed->d_a2v.rnfl = rnfl * foliar_projected_cover;
  double downpar = ed->d_a2v.pardown;             //note: already foliar_projected_cover adjusted
  ed->d_v2a.evap_moss = 0.;                       // v5l S4 (stays 0 unless nv_rain_surface == 1)

  // v5l ET-split switches for this PFT and day (all false when the switches are 0)
  // Refs (S4): moss as a ground/soil layer (Chadburn et al. 2015, doi:10.5194/gmd-8-1493-2015; Yi et al. 2009, doi:10.1029/2008JG000841); moss/understory evaporation (Heijmans et al. 2004, doi:10.1016/j.agrformet.2003.12.006).
  const bool nv_moss_evap = (nv_rain_surface == 1 && nonvasc > 0);   // S4 mode 1
  const bool nv_to_ground = (nv_rain_surface == 2 && nonvasc > 0);   // S4 mode 2
  // Refs (S2): snow interception is a canopy process (Hedstrom & Pomeroy 1998, doi:10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4);
  //   low vegetation buried by snow (Sturm et al. 2001, doi:10.1175/1520-0442(2001)014<0336:SSIIAT>2.0.CO;2); moss layer in the soil column (Yi et al. 2009, doi:10.1029/2008JG000841).
  const bool nv_buried = (nv_snow_bury_fix == 1 && nonvasc > 0
                          && snow_depth > nv_snow_bury_depth);       // S2
  const bool tran_fpc_once = (fpc_once_fix == 1 && foliar_projected_cover > 0.); // S1

  if(leaf_area_index >0) {
    //precipitation interception and throughfall
    // Refs (S3): p_liq = 0.1 mm per unit (LAI+SAI) (Dickinson et al. 1993, doi:10.5065/D67W6959; CLM5 Tech Note eq. 2.7.12); S_max(LAI) of
    //   von Hoyningen-Huene 1981 (DVWK report, no DOI) as in LISFLOOD (van der Knijff et al. 2010, doi:10.1080/13658810802549154); TVC interception
    //   observations: Zwieback et al. 2019, doi:10.1088/1748-9326/ab1049; Chang et al. 2022, doi:10.1016/j.rse.2021.112747.
    double rstore0 = 0.;   // S3: canopy water stored from previous days (mm)
    if (rain_int_lai_fix == 1) {
      // S3 (v5l): capacity = ircoef [mm per unit LAI] * LAI, minus what is still stored
      double rcap = envpar.ircoef * leaf_area_index;
      rstore0 = fmax(ed->d_vegs.rwater, 0.);
      ed->d_a2v.rinter = fmin(ed->d_a2v.rnfl, fmax(rcap - rstore0, 0.));
    } else {
      ed->d_a2v.rinter = getRainInterception(ed->d_a2v.rnfl, leaf_area_index);
    }
    if (nv_to_ground) {
      ed->d_a2v.rinter = 0.;   // S4 mode 2: rain on moss/lichen reaches the ground surface
    }
    ed->d_a2v.sinter = getSnowInterception(ed->d_a2v.snfl, leaf_area_index);
    if (nv_buried) {
      ed->d_a2v.sinter = 0.;   // S2: moss/lichen under the snowpack intercept no snow
    }
    // water available on the wet canopy today: legacy = today's interception only
    double wetw = ed->d_a2v.rinter;
    if (rain_int_lai_fix == 1) { wetw = rstore0 + ed->d_a2v.rinter; }
    if (nv_to_ground) { wetw = 0.; }   // S4 mode 2: no wet moss "canopy"; stored water drips below
    ed->d_v2g.rthfl = ed->d_a2v.rnfl- ed->d_a2v.rinter;
    ed->d_v2g.sthfl = ed->d_a2v.snfl -ed->d_a2v.sinter;
    //evapotranspiration
    //temperature and pressure correction factor for conductances
    double gcorr = 1.; //pow( (atmsd->ta +273.15)/293.15, 1.75); //* 101300/pa;
    double gl_st = 0.;
    ed->d_vegd.m_vpd=0.;

    if(ed->d_vegd.btran>0) {
      if (tran_fpc_once) {   // S1: PAR per unit covered area
        gl_st = getLeafStomaCond(ta,downpar/foliar_projected_cover*EPAR, vpd, ed->d_vegd.btran,
                                 ed->d_vegd.m_ppfd, ed->d_vegd.m_vpd);
      } else {
        gl_st = getLeafStomaCond(ta,downpar*EPAR, vpd, ed->d_vegd.btran,
                                 ed->d_vegd.m_ppfd, ed->d_vegd.m_vpd);
      }
    }

    gl_st *= gcorr;
    double gl_bl = envpar.gl_bl; //boundary layer conductance (projected area basis) m/s
    gl_bl *= gcorr;
    double gl_c  = envpar.gl_c; // cuticular conductance m/s
    gl_c *=gcorr;
    double gl_t_wv = (gl_bl *gl_st)/(gl_bl+gl_st+gl_c);
    double gl_sh = gl_bl;
    double gl_e_wv = gl_bl;
    double gc_e_wv = gl_e_wv * leaf_area_index;
    double gc_sh =gl_sh * leaf_area_index;
    double gl_t_wv_pet = (gl_bl *envpar.glmax)/(gl_bl+envpar.glmax+gl_c);
    double gl_sh_pet = gl_bl;
    double gl_e_wv_pet = gl_bl;
    double gc_e_wv_pet = gl_e_wv_pet * leaf_area_index;
    double gc_sh_pet =gl_sh_pet * leaf_area_index;
    double sw =ed->d_a2v.swinter;
    double daytimesw = sw;
    double rainsw = sw;
    double rv, rh;
    double rv_pet, rh_pet;
    double vpdpa = vpd; //Pa
    // S1: transpiration radiation per unit covered area (tran *= FPC once, in Cohort)
    double tran_sw = sw;
    if (tran_fpc_once) { tran_sw = sw / foliar_projected_cover; }

    if(wetw>0.) {
      rv = 1./gc_e_wv;
      rh = 1./gc_sh;
      rv_pet = 1./gc_e_wv_pet;
      rh_pet = 1./gc_sh_pet;
      double et1 = getPenMonET(ta,vpdpa, rainsw, rv, rh);
      double et1_pet= getPenMonET(ta,vpdpa, rainsw, rv_pet, rh_pet);
      double dayl1 = wetw/et1;

      if(dayl1>daylsec) {
        if(daylsec>0.) {
          ed->d_v2a.tran = 0.;
          ed->d_v2a.evap = et1 *daylsec;
        } else {
          ed->d_v2a.tran=0.;
          ed->d_v2a.evap=0.;
        }
      } else {
        ed->d_v2a.evap = wetw;
        daylsec -= dayl1;
        rv = 1.0/gl_t_wv;
        rh =1.0/gl_sh;
        double et2;
        if (tran_fpc_once) {
          et2 = getPenMonET(ta, vpdpa, tran_sw, rv, rh);
        } else {
          et2 = getPenMonET(ta, vpdpa, rainsw, rv, rh);
        }
        ed->d_v2a.tran = et2 * daylsec;
      }

      double dayl1_pet = wetw/et1_pet;

      if(dayl1_pet>daylsec) {
        if(daylsec>0.) {
          ed->d_v2a.tran_pet = 0.;
          ed->d_v2a.evap_pet = ed->d_v2a.evap;
        } else {
          ed->d_v2a.tran_pet = 0.;
          ed->d_v2a.evap_pet = 0.;
        }
      } else {
        ed->d_v2a.evap_pet = wetw;
        daylsec -= dayl1_pet;
        rv_pet = 1.0/gl_t_wv_pet;
        rh_pet =1.0/gl_sh_pet;
        double et2_pet= getPenMonET(ta,vpdpa, rainsw, rv_pet, rh_pet);
        ed->d_v2a.tran_pet = et2_pet * daylsec;
      }
    } else { // no interception
      ed->d_v2a.evap = 0.;
      rv = 1.0/gl_t_wv;
      rh =1.0/gl_sh;
      double et3;
      if (tran_fpc_once) {
        et3 = getPenMonET(ta, vpdpa, tran_sw, rv, rh);
      } else {
        et3 = getPenMonET(ta, vpdpa, daytimesw, rv, rh);
      }
      ed->d_v2a.tran = et3 * daylsec;
      ed->d_v2a.evap_pet = 0.;
      rv_pet = 1.0/gl_t_wv_pet;
      rh_pet =1.0/gl_sh_pet;
      double et3_pet= getPenMonET(ta, vpdpa, daytimesw, rv_pet, rh_pet);
      ed->d_v2a.tran_pet = et3_pet * daylsec;
    }

    ed->d_v2a.sublim = getCanopySubl(ed->d_a2v.swdown,ed->d_a2v.sinter, leaf_area_index);
    ed->d_vegd.cc = gc_e_wv;
    ed->d_vegd.rc = 1./ed->d_vegd.cc;
    ed->d_vegs.snow  += (ed->d_a2v.sinter - ed->d_v2a.sublim);
    if (rain_int_lai_fix == 1) {
      // S3: storage = stored + intercepted - evaporated (nothing is lost)
      ed->d_vegs.rwater = rstore0 + ed->d_a2v.rinter - ed->d_v2a.evap;
    } else {
      ed->d_vegs.rwater+= (ed->d_a2v.rinter - ed->d_v2a.evap);
    }
    ed->d_v2g.sdrip = 0.0;
    double maxvegsnow = 0.10*leaf_area_index; //that 0.10 LAI mm snow on vegetation is
                                     //arbitrary - needs more mechanism
                                     //to do this
    if (nv_buried) { maxvegsnow = 0.; }   // S2: stored snow falls onto the snowpack

    if (ed->d_vegs.snow>maxvegsnow) {
      ed->d_v2g.sdrip = ed->d_vegs.snow - maxvegsnow;
    }

    ed->d_vegs.snow -= ed->d_v2g.sdrip;
    ed->d_v2g.rdrip = 0.0;
    double maxvegrain = 0.05*leaf_area_index; //that 0.05 LAI mm rain storage is
                                     //arbitrary - needs more mechanism
                                     //to do this
    if (rain_int_lai_fix == 1) { maxvegrain = envpar.ircoef * leaf_area_index; } // S3 capacity
    if (nv_to_ground) { maxvegrain = 0.; }   // S4 mode 2: stored water drips to the ground

    if (ed->d_vegs.rwater>maxvegrain) {
      ed->d_v2g.rdrip = ed->d_vegs.rwater - maxvegrain;
    }

    ed->d_vegs.rwater -= ed->d_v2g.rdrip;
    if (nv_moss_evap) {
      // S4 mode 1: evaporation of the moss/lichen surface water is moss-surface
      // evaporation (EVAP_MOSS, still part of EET), not canopy evaporation
      ed->d_v2a.evap_moss = ed->d_v2a.evap;
      ed->d_v2a.evap = 0.;
    }
  } else {   //leaf_area_index <=0, i.e., no vegetation?
    ed->d_vegd.cc =0.;
    ed->d_vegd.rc =0.;
    ed->d_a2v.rinter =0.;
    ed->d_a2v.sinter =0.;
    ed->d_v2a.sublim =0.;
    ed->d_v2a.tran   =0.;
    ed->d_v2a.evap   =0.;
    ed->d_v2a.tran_pet =0.;
    ed->d_v2a.evap_pet =0.;
    ed->d_v2a.sublim =0.;
    ed->d_v2g.rdrip = 0;
    ed->d_v2g.sdrip = 0;
    ed->d_v2g.rthfl = 0;
    ed->d_v2g.sthfl = 0;
    ed->d_vegs.snow  =0.;
    ed->d_vegs.rwater=0.;
  }
};

double Vegetation_Env::getPenMonET(const double & ta, const double& vpd,
                                   const double &irad, const double &rv,
                                   const double & rh) {
  double et; // out , mmH2O/m2s
  double CP =1004.64 ; // specific heat capacity of dry air [J/kgK)
  double tk = ta+273.15;
  double pa = 101300;// pressure , Pa
  double rho = 1.292- (0.00428 * ta); // air density  kg/m3
  double EPS=0.6219; // ratio of mole weights
  double SBC= 5.67e-8; //Stefan-boltzmann constant W/m2K4
  /*resistance to raiative heat transfer through air*/
  double rr = rho * CP /(4.0 * SBC * tk* tk*tk);
  /* resistance to convective heat tranfer: rh*/
  /*resistance to latent heat transfer rv*/
  /*combined resistance to convectie and radiative heat transfer,
   * parallel resistances:rhr= (rh*rr)/(rh+rr)*/
  double rhr = (rh*rr)/(rh+rr);
  /*latent heat of vaporization as a function of ta*/
  double lhvap = 2.5023e6 -2430.54 *ta;
  double dt =0.2;
  double t1 = ta+dt;
  double t2 =ta-dt;
  /*saturated vapor pressure at t1, t2*/
  double pvs1 = 610.7 *exp(17.38 *t1/(239.+t1));
  double pvs2 = 610.7 *exp(17.38 *t2/(239.+t2));
  /*slope of pvs vs. T curve at T*/
  double slope = (pvs1-pvs2)/(t1-t2);
  /*evapotranspiration*/
  et = (slope*irad+ rho*CP *vpd/rhr)/((pa * CP *rv)/(lhvap*EPS *rhr)+slope);
  return et/lhvap;
};

double Vegetation_Env::getCanopySubl(const double & rac, const double & sinter,
                                     const double & lai ) {
  double sub;
  double psub; //potential sub
  double snow_int =envpar.iscoef;
  double lamdaw = 2.501e6; // latent heat of vaporization J/kg
  double lf = 3.337e5 ;// latent heat of fusion J/kg
  //rac in unit W/m2
  // change the unit from J/kg to MJ/mm
  double phasechange = (lamdaw+lf)/10e6 ;
  double sub1 = (rac*86400)/phasechange;
  //double sub2 = lai/envpar.all2prj * snow_int;
  double sub2 = lai * snow_int;

  if(sub1 <sub2) {
    psub = sub1;
  } else {
    psub =sub2;
  }

  if(sinter<=psub) {
    sub =sinter;
  } else {
    sub = psub;
  }

  return sub;
};

double Vegetation_Env::getLeafStomaCond(const double & ta,
                                        const double & ppfdabsorb,
                                        const double & vpdin,
                                        const double& btran, double & m_ppfd,
                                        double & m_vpd ) {
  double gl; // leaf conductance m/s, per unit LAI
  double m_psi; // soil water matric potential effect
  double m_co2;
  double m_tmin;
  double m_tot; // total effect of different factos
  double vpd_open = envpar.vpd_open ; //Pa , start of conductance reduction
  double vpd_close =envpar.vpd_close; //Pa , complete conductance reduction
  double glmax = envpar.glmax; //m/s, maximum stomata conductance
  m_ppfd= ppfdabsorb/(ppfdabsorb+envpar.ppfd50);
  double vpd = vpdin; //

  if(vpd <vpd_open) {
    m_vpd = 1.0;
  } else if(vpd>vpd_close) {
    m_vpd = 0.;
  } else {
    m_vpd =(vpd_close-vpd)/(vpd_close-vpd_open);
  }

  m_psi = btran;
  m_co2=1;
  double tmin = ta-5;

  // here ta is average air temperature (degC)  ;
  // 5 is derived from mean air temperature and minimum air temperature
  // difference from station data of Fairbanks
  if(tmin>0) {
    m_tmin=1;
  } else if(tmin<-8) {
    m_tmin=0.;
  } else {
    m_tmin = 1.0 +(0.125*tmin);
  }

  m_tot = m_ppfd * m_vpd * m_co2* m_tmin*m_psi;
  gl =  m_tot * glmax;
  return gl;
};

double Vegetation_Env::getRainInterception(const double & rain,
                                           const double & lai) {
  if (interception_scheme == 0) { return getRainInterceptionLegacy(rain, lai); }
  // interception_scheme 1 (upstream regional_wetland_deepc): "Based on CLM5.0 Technical Note"
  double rinter = rain * envpar.ircoef * (1 - tanh(-0.5 * lai));

  if(rain <= rinter) {
    rinter = rain;
  }

  return rinter;
};

double Vegetation_Env::getSnowInterception(const double & snow,
                                           const double & lai) {
  if (interception_scheme == 0) { return getSnowInterceptionLegacy(snow, lai); }
  // interception_scheme 1 (upstream regional_wetland_deepc): "Based on CLM5.0 Technical Note"
  double sinter = snow * envpar.iscoef * (1 - exp(-0.5 * lai));

  if(snow <= sinter) {
    sinter = snow;
  }

  return sinter;
};

// interception_scheme 0 (default): pre-merge (legacy/v5q) interception, verbatim
double Vegetation_Env::getRainInterceptionLegacy(const double & rain,
                                                 const double & lai) {
  // input: rain mm/day
  double rinter;
  double raincm = rain/10.; // convert to cm
  //20200930 We do not have a citation for these equations
  // and do not know why ircoef from the parameter files is not used.
  double max_int = raincm - ((raincm*0.77-0.05)+0.02 * raincm);
  // may need to add LAI adjustment?
  //
  max_int *= 10. ; // convert back to mm;

  if(rain <= max_int) { // all intercepted
    rinter = rain;
  } else { //partly intercepted
    rinter= max_int;
  }

  return rinter;
};

double Vegetation_Env::getSnowInterceptionLegacy(const double & snow,
                                                 const double & lai) {
  double sinter;
  double psinter; //potential snow interception
  double ISmax = envpar.iscoef ; //this is parameter for snow interception
                                 //  by canopy
  //it should be vegetation specific [mm /LAI /day ],
  //  e.g., for tussock tundra set to 0
  psinter = ISmax * lai;

  if(psinter >= snow) {
    sinter = snow;
  } else {
    sinter = psinter;
  }

  return sinter;
};

void Vegetation_Env::setCohortLookup(CohortLookup * chtlup) {
  chtlu = chtlup;
};

void Vegetation_Env::setEnvData(EnvData* edatap) {
  ed = edatap;
};

void Vegetation_Env::setFirData(FirData* fdp) {
  fd =fdp;
};