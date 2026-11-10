#include "../include/Snow_Env.h"

int Snow_Env::albedo_fix = 0; // Fix B4 (v5k), default off (legacy)

Snow_Env::Snow_Env(): rad_melt(0.) {
};

Snow_Env::~Snow_Env() {
};

void Snow_Env::updateDailyM(const double & tdrv) {
  //update tsurface with nfactor
  double tsurface;
  // DRAFT (daily-structural 3/4): restore the n-factor on the snow surface
  // (was 'tsurface = tdrv; // TESTING ... snow nfactor = 1'). With
  // nfactor(w) = 1 this only matters on days where the summer n-factor is in
  // effect (n-factor is chosen from the previous day's Tair, Cohort.cpp ~503).
  tsurface = tdrv *ed->d_soid.nfactor;
  if (ground->seb_on) { // Option 4 (v5j) SEB
    tsurface = ground->seb_ts;
  }
  //dsmass is the total snowfall during one timestep:
  double dsmass = ed->d_v2g.sthfl + ed->d_v2g.sdrip //these 2 items has already
                                                    //  been adjusted by FPC
                  + (1.0 - cd->m_vegd.fpcsum) * ed->d_a2l.snfl;

  //Note unit converstion: 1 mm H2O = 1 kgH2O/m2
  bool slchg1 = false;
  rad_melt = 0.;

  if (!ground->toplayer->isSnow && ground->toplayer->tem>0.01) {
    //melting all snowfalling ('dsmass')
    // PROPOSAL: pending sub-layer snow ('extramass') melts too; it was discarded before.
    ed->d_snw2soi.melt = dsmass + fmax(ground->snow.extramass, 0.);
    ground->snow.extramass = 0.;
  } else { // otherwise, construct a new snow layer, if any
    //add a new snow layer as the new toplayer
    slchg1 = ground->constructSnowLayers(dsmass, tsurface);

    if(slchg1)  {
      ground->resortGroundLayers();
      initializeNewSnowState(); //for new snow layer when it is a front layer
    }
  }

  // snow-melting/sublimating
  bool slchg2 = false;

  if(ground->toplayer->isSnow) {
    //sublimating of existed surface snow layers accompanying snow surface
    //  energy budget process (Yuan: need further thinking here ???)
    updateDailySurfFlux(ground->toplayer, tsurface);
    // Option C: with thermal_solver = 1 the radiation/SEB melt energy already went into the
    // snow cells in the thermal step (melt water is collected from liq below); no mass removal.
    if (thermal_solver == 1) { rad_melt = 0.; }
    // melting below snow layers, if any
    ed->d_snw2soi.melt = meltSnowLayersAfterT(ground->toplayer);
    // DRAFT: radiation-driven melt (computed in updateDailySurfFlux) goes to
    // the soil (infiltration/runoff), not to the atmosphere.
    ed->d_snw2soi.melt += rad_melt;
  }

  //here 'dsmass' is for snow-melting and sublimating
  // PROPOSAL: melt water was already taken out of the snow layers by meltSnowLayersAfterT()
  // (ice->liq in Stefan, then liq zeroed). Removing it again here double-counted the melt
  // (same issue as the sublimation duplicate noted in updateDailySurfFlux). Only sublimation remains.
  // DRAFT: mass removed from the pack = true sublimation + radiation-driven melt.
  dsmass = ed->d_snw2a.sublim + rad_melt;
  //Note unit converstion: 1 mm H2O = 1 kgH2O/m2

  //removal of melted snow layers
  slchg2 = ground->constructSnowLayers(-dsmass, tsurface);

  if(slchg2)  {
    ground->resortGroundLayers();
  }

  // snow layer adjusting, if needed
  bool slchg3 = ground->combineSnowLayers();

  if(slchg3) {
    ground->resortGroundLayers();
  }

  bool slchg4 = ground->divideSnowLayers();

  if(slchg4) {
    ground->resortGroundLayers();
  }

  // after snow layer adding/removal/adjusting, needs updating thickness
  ground->updateSnowLayerPropertiesDaily();   // each layer
  if (ground->snow_scheme == 1) { ground->stratifySnowLayers(true); } // Schaefer 2009 two-layer profile
  ground->updateSnowHorizon();                // all snow horizon
  ground->checkWaterValidity();
  updateSnowEd(ground->toplayer);
};

void Snow_Env::updateDailySurfFlux( Layer* toplayer, const double & tdrv) {
  if (toplayer->isSnow) {
    double albnir = getAlbedoNir(toplayer->tem);
    double albvis = getAlbedoVis(toplayer->tem);
    double insw =  ed->d_v2g.swthfl * cd->m_vegd.fpcsum
                   + ed->d_a2l.nirr *(1.- cd->m_vegd.fpcsum);
    ed->d_snw2a.swrefl = insw *0.5 * albnir + insw *0.5 * albvis;
    double rn = insw- ed->d_snw2a.swrefl;
    double pradmelt = 0.;
    double sublim;
    if (ground->seb_on) { // Option 4 (v5j) SEB: melt and sublimation from the SEB
      sublim = ground->seb_sublim;
      pradmelt = ground->seb_melt;
    } else {
      sublim = getSublimation(rn, ed->d_snws.swesum, tdrv, pradmelt); //mm SWE/day, or kgSW/day
    }
    Layer * currl=toplayer;
    double actsub = 0.;
    double actmelt = 0.;

    while (currl!=NULL) {
      if (currl->isSnow && (sublim>0. || pradmelt>0.)) {
        double avail = currl->ice;
        double layicermv = fmin(sublim, avail);
        actsub += layicermv;
        sublim -= layicermv;
        avail -= layicermv;
        double laymelt = fmin(pradmelt, avail);
        actmelt += laymelt;
        pradmelt -= laymelt;
        //ice is removed in Ground::constructSnowLayers() (see updateDailyM)
      } else {
        break;
      }

      currl = currl->nextl;
    }

    ed->d_snw2a.sublim = actsub;
    rad_melt = actmelt;
  }

  /*  // the following will cause confusion OR duplicated snowmelting process
      //see Zhuang et al., 2004 D4 not sure whether equation D8 is right, when compared with the one from
      //Brubaker et al., 1996
     // double melt = 2.99*rn/0.2388 *86400/10000.+2.0*tdrv; //mm/day
      double melt =0.26*rn +2.0*tdrv;//mm/day

      if(melt>0 && melt<=sublim){
         ed->d_snw2soi.melt = melt;
         ed->d_snw2a.sublim = sublim-melt;

      }else if(melt<=0){
           ed->d_snw2soi.melt = 0.;
           ed->d_snw2a.sublim = sublim;
      }else if (melt>sublim){
           ed->d_snw2a.sublim = 0.;
         ed->d_snw2soi.melt = sublim;
      }
  */
};

// get albedo of visible radiation on snow
double Snow_Env::getAlbedoVis(const double & tem) {
  double vis;

  if(tem<-10) {
    // Refs (B4): no published source; hard-coded 0.9 (upstream 419a4eb3, fix commit 07f6c8d0); restores consistency with snwalbmax and continuity at -10 C (cf. Warren 1982, doi:10.1029/RG020i001p00067).
    vis = (albedo_fix == 1) ? (double)snowenvpar.albmax : 0.9; // Fix B4 (v5k)
  } else {
    vis = snowenvpar.albmax - (snowenvpar.albmax -snowenvpar.albmin)
          * (tem+10)/10;
    vis = fmin(snowenvpar.albmax, vis);
    vis = fmax(snowenvpar.albmin, vis);
  }

  return vis;
};

// get albedo of invisible radiation
double Snow_Env::getAlbedoNir(const double & tem) {
  double nir;
  nir = getAlbedoVis(tem);
  return nir;
};

// modified - called after 'Stefan.cpp:: processFrontSnowLayers()'
// Here, the 'liq' in snowlayer will be collected as melting water
// and removed from the layer

double  Snow_Env::meltSnowLayersAfterT(Layer* toplayer) {
  double melt=0;   // unit: mm H2O
  Layer* currl ;
  currl= toplayer;

  while(currl!=NULL) {
    if(currl->isSnow) {
      if(currl->frozen==-1) {
        melt += currl->liq;
        currl->liq = 0.;
      } else {
        if(currl->frozen==0) {
          melt += currl->liq; //Note: 'liq' here is the conversion
                              //  of 'ice' due to partial snow-melting
          currl->liq = 0.; // restore 'liq' to zero

          if (currl->tem <0.) {
            currl->frozen=1;
          } else {
            currl->frozen=-1;
          }
        }
      }

      currl=currl->nextl;
    } else {
      break;
    }
  }

  return melt;    //Note: 1kgH2O/m2 = 1 mm H2O
};

// assign double-linked snow horizon data to 'ed'
void Snow_Env::updateSnowEd(Layer *toplayer) {
  ed->d_snws.swesum = 0.;
  ed->d_snws.snowthick = 0.;
  Layer* currl=toplayer;
  int snowind = 0;

  while(currl!=NULL) {
    if(currl->isSnow) {
      ed->d_snws.snwice[snowind] = currl->ice;
      ed->d_snws.snwliq[snowind] = currl->liq;
      ed->d_snws.swe[snowind] = currl->ice;
      ed->d_snws.swesum += currl->ice;
      ed->d_snws.snowthick += currl->dz;
      currl=currl->nextl;
    } else {
      if (snowind>=MAX_SNW_LAY) {
        break;
      }

      ed->d_snws.swe[snowind] = MISSING_D; //in this way, 'swe' will be
                                           //  refreshed for all
      ed->d_snws.snwice[snowind] = MISSING_D; //in this way, 'swe' will be
                                              //  refreshed for all
      ed->d_snws.snwliq[snowind] = MISSING_D; //in this way, 'swe' will be
                                              //  refreshed for all
    }

    snowind++;
  }

  ed->d_snws.extraswe= ground->snow.extramass;
  ed->d_snws.swesum += ground->snow.extramass;

  if (ed->d_snws.swesum>0.) {
    //assuming no evap and reflection from soil surface, if snow layer
    //  still exists or is constructed
    ed->d_soi2a.evap   = 0.; //Yuan: these might not be needed? but no
                             //  harmless here to do so
    ed->d_soi2a.swrefl = 0.;
    ed->d_soi2l.qover  = 0.;
  }
};

void Snow_Env::initializeParameter() {
  snowenvpar.albmax = chtlu->snwalbmax;
  snowenvpar.albmin = chtlu->snwalbmin;
};

void Snow_Env::initializeNewSnowState() {
  Layer *toplayer = ground->toplayer;

  while(toplayer!=NULL) {
    if(toplayer->isSnow) {
      //other 'Layer' properties already updated in 'Ground.cpp',
      //  when constructing this layer
      toplayer->frozen = 1;
      toplayer->maxliq = 0.;
      toplayer->maxice = toplayer->dz*ground->snowdimpar.denmax;
      toplayer->liq = 0.;
      toplayer->ice = toplayer->dz*toplayer->rho;
      break;
    } else {
      break;
    }
  }
};

void Snow_Env::set_state_from_restartdata(const RestartData & rdata) {
  //set the parameters
  double TSsnow[MAX_SNW_LAY];
  double ICEsnow[MAX_SNW_LAY];

  for (int i=0; i<MAX_SNW_LAY; i++) {
    TSsnow[i] = rdata.TSsnow[i];
    ICEsnow[i] = rdata.ICEsnow[i];
  }

  Layer* currl = ground->toplayer;
  int snind = -1;
  ed->d_snws.swesum = 0;
  ed->d_snws.snowthick = 0;

  while(currl!=NULL) {
    if(currl->isSnow) {
      snind ++;
      currl->tem = TSsnow[snind];
      currl->ice = ICEsnow[snind];
      currl->liq = 0.;
      currl->frozen = 1;
      currl->maxliq = 0.;
      currl->maxice = currl->dz*ground->snowdimpar.denmax;
    } else {
      break;
    }

    currl = currl->nextl;
  }

  //
  //ground->checkWaterValidity();
};

double Snow_Env::getSublimation(double const & rn, double const & swe,
                                double const & ta, double & radmelt) {
  // rn: radiation. W/m2
  // swe: snow water equivalent. mm
  // LHVAP: latent heat of vaporization. J/kg
  // LHFUS: latent heat of fusion. J/kg
  // output sublimation. mm/day
  double sub = 0.;
  radmelt = 0.;
  double sabsorb =0.6;     // radiation absorptivity of snow, Zhuang 0.6 ?
  double rabs = rn * sabsorb; //W/m2

  if( swe>0. && ta<=-1) {
    double psub= rabs*86400/(LHVAP + LHFUS);

    if(psub>swe) {
      sub = swe;
    } else {
      sub = psub;
    }
  } else if(swe>0. && ta>=-1) {
    // DRAFT (daily-structural): this is radiation-driven MELT, not sublimation.
    // It used to be returned as 'sub' (removed to the atmosphere and added to
    // EET). Return it via 'radmelt' so the caller routes it to melt water.
    double pmelt= rabs*86400/(LHFUS);

    if(pmelt>swe) {
      radmelt = swe;
    } else {
      radmelt = pmelt;
    }
    sub = 0.;
  } else {
    sub = 0.;
  }

  return sub;
};

void Snow_Env::checkSnowLayersT(Layer* toplayer) {
  Layer* currl = toplayer;
  SnowLayer* sl;

  while(currl!=NULL) {
    if (currl->isSnow) {
      sl = dynamic_cast<SnowLayer*>(currl);

      if(currl->tem>0) {
        currl->tem =-0.01;
      }
    }

    currl=currl->nextl;
  }
};

void Snow_Env::setGround(Ground* grndp) {
  ground = grndp;
};

void Snow_Env::setCohortLookup(CohortLookup* chtlup) {
  chtlu = chtlup;
};

void Snow_Env::setCohortData(CohortData* cdp) {
  cd = cdp;
};

void Snow_Env::setEnvData(EnvData* edp) {
  ed = edp;
};



// Option C (thermal_solver = 1), see Snow_Env.h.
double Snow_Env::radMeltPotential(const double & tsurface) {
  if (ground == NULL || ground->toplayer == NULL || !ground->toplayer->isSnow) { return 0.; }
  Layer* toplayer = ground->toplayer;
  double albnir = getAlbedoNir(toplayer->tem);
  double albvis = getAlbedoVis(toplayer->tem);
  double insw =  ed->d_v2g.swthfl * cd->m_vegd.fpcsum
                 + ed->d_a2l.nirr *(1.- cd->m_vegd.fpcsum);
  double swrefl = insw *0.5 * albnir + insw *0.5 * albvis;
  double rn = insw - swrefl;
  double pradmelt = 0.;
  getSublimation(rn, ed->d_snws.swesum, tsurface, pradmelt);
  return pradmelt;
}
