/**
 * @file FW_Interface.cpp
 * \author Joshua M. Rady
 * Woodwell Climate Reseach Center
 * \date 2024
 *
 * @brief This file contains functions to connect the DVM-DOS-TEM WildFire class to the process based
 * wildfire model components.
 *
 * This is currently a very rough draft, possibly a placeholder.  Everything including the name
 * should be expected to change.
 * 
 * All of the functions requiring direct access to WildFire members have been moved into the class
 * but have been left here pending testing so the new code can be more clearly identified.  There are
 * number of attendant functions, most of which can also be moved into the class but are being left
 * as is for now to make dependancies stark, i.e. they do not need access to member data or
 * functions.
 */

#include "../include/FW_Interface.h"
#include "../include/TEMLogger.h"
#include "../include/TEMUtilityFunctions.h"//For length_of_day().
#include "../include/WildFire.h"
#include "../include/Layer.h"

#include "FireweedCrownFireScottReinhardt.h"
#include "FireweedDeadFuelMoistureFosberg.h"
#include "FireweedFuelTools.h"
#include "FireweedLiveFuelMoistureGSI.h"
#include "FireweedMetUtils.h"
#include "FireweedUtils.h"

#include "../include/GroundFire.h"

#include <cmath>//For fmin() & isnan().
#include <iomanip>
#include <limits>

extern src::severity_logger< severity_level > glg;

//Constants:
const double c2b = 2.0;//The carbon to biomass multiplier for vegetation on a dry basis.
//This could vary by PFT but many models use a single value.
const double gPerKg = 1000;//Move to FireweedUnits.h?

/** Calculate wildfire behavior and effects using the modeled vegetation, fuels, and meteorology.
 *
 * This is the main entry point for the process based wildfire model.
 * The function needs the ecosystem state, meteorology, and time of year as inputs.
 *
 * JMR_NOTE: I attempted multiple ways to pass this information in but ultimately found it necessary to
 * modify WildFire to provide access to most of the necessary inputs.
 * 
 * @param[in] monthIndex The current month as a zero based index.
 *
 * @returns The soil burn depth from ground fire (meters).  Other output are stored in model objects.
 * 
 * @note As discussed in a note in WildFire::burn() the burn depth is stored in the FirData member so the
 * return value is not strictly needed.  We are keeping this code parallel to getBurnOrgSoilthick()
 * for now.
 * @note Used in the process wildfire model only.
 */
double WildFire::ProcessWildfire(const int monthIndex)//Name could change.
{
  BOOST_LOG_SEV(glg, debug) << "Entering WildFire::ProcessWildfire()...";

  double burnDepth = 0.0;//Return value.

  //Gather weather and environmental conditions:---------------------------
  double tempAir = GetAirTemperature();
  //Current humidity is needed for calculating fuel moisture but that code handles it itself.
  //We may also need it for duff moisture soon.
  
  double windSpeed = GetMidflameWindSpeed();
  
  //The percent slope is stored in the CohortData object and also in the wildfire object:
  double slopeSteepness = SlopePctToSteepness(cd->cell_slope);

  double snowDepth_cm = ground->snow.thick * 100.0;//Snow can prevent fire.

  //Shortwave radiation may be needed for more advanced moisture calculations added in the future.

  //Determine the fuel model for the location and set its parameters:------

  //Get the CMT for the grid cell:
  //JMR_NOTE: The CMT number is in the CohortData member of the cohort.  The WildFire object maintains a
  //pointer to it's parent cohort data but it is private.  This code needs to be in the class,
  //a friend or have access to the cohort.
  int theCMTnumber = cd->cmttype;

  //Crosswalk from CMT to the matching fuel model:
  int fuelModelNumber = GetMatchingFuelModel(theCMTnumber);

  //Load the fuel model from the fuel model table file:
  FuelModel fm = GetFuelModelFromCSV(md.fire_fuel_model_file, fuelModelNumber);//Or tab delimited!!!!!
  fm.ConvertUnits(Metric);//Convert to metric units.

  siteFM = fm;//Save the input for use in getAbgVegetationBurntFractionsProcess().

  //If the snow depth is above the specified limit fire does not initiate.  The extinguishement
  //depth may be non-zero.  We propose very small amounts of snow may not prevent a fire.  Also we
  //only have one measure of snow thickness for the whole site.  In reality where there is a small
  //amounts of snow on the ground it is likley that some areas will be bare and more condusive to
  //a fire starting.  The appropriate extinguishment depth is an area that needs more research. 
  //This has nothing to do with snow persisting under snow, AKA overwintering or zombie fires, which
  //are not yet represented in the model.
  //If the fire does not initiate we can skip the fire calculations:
  if ((md.fire_max_snow >= 0.0) && (snowDepth_cm > md.fire_max_snow))
  {
    BOOST_LOG_SEV(glg, info) << "Snow is too deep (" << snowDepth_cm << " cm) for ignition to start a fire.";
    
    //For the remaining fire functions to complete properly we need to:
    // - Store and return the burn depth/thickness (below).
    // - Save the fuel model for use in getAbgVegetationBurntFractionsProcess().  We won't need the
    //loading and fuel moisture information there so we can stop before caculating that.
    siteFM = fm;
    // - Save a BurnupSim object for use in getAbgVegetationBurntFractionsProcess().  Since the fire
    //isn't burning we don't run Burnup but return an empty BurnupSim object with an flag
    //indicating no fire burned.  The flag needs to be checked by subsequent code that uses siteBU
    //since in this case the object's contents is not valid.
    BurnupSim mockOutput;
    siteBU = mockOutput;
  }
  else
  {
    //Determine the surface fuels from the model vegetation and soil states and update the fuel
    //loadings from their default values:
    CohortStatesToFuelLoading(fm, md.fire_moss_as_dead_fuel);

    //Save the fuel loading prior to fire: (will be compared below...)
    std::vector <double> fuelLoadingBefore = fm.w_o_ij;

    //Calculate the fuel bed depth:
    CalculateFuelBedDepth(fm, md.fire_calculate_delta);

    //Calculate fuel moisture:
    //Note: It is better to calculate fuel moisture after calculating fuel loadings since that process
    //might change the fuel sizes.
    std::vector <double> M_f_ij = CalculateFuelMoisture(fm, monthIndex);

    //Add the moisture to the fuel model possibly computing dynamic fuel moisture:
    BOOST_LOG_SEV(glg, debug) << "Apply fuel moisture to fuel model...";
    if (md.fire_dynamic_fuel)
    {
      fm.CalculateDynamicFuelCuring(M_f_ij);
    }
    else
    {
      fm.SetFuelMoisture(M_f_ij);
    }

    //Feed fuels and weather conditions into the surface fire models:------

    //Dump the fuel model if debugging:  This may be temporary?????
    BOOST_LOG_SEV(glg, debug) << "Dumping the fuel model prior to fire:";
    BOOST_LOG_SEV(glg, debug) << fm;

    //First run the Rothermel & Albini spread rate model:
    //It takes a fuel model (and attendant data) and returns the calculation details.
    //M_f_ij does not need to be included since it is added to the fuel model object above.
    BOOST_LOG_SEV(glg, debug) << "Perform surface file spread rate calculations...";
    SpreadCalcs raData = SpreadCalcsRothermelAlbini_Het(fm, windSpeed, slopeSteepness);
                                                        //M_f_ij,//fuel moisture
                                                        //false,//useWindLimit
                                                        //FALSE);//debug

    //Dump the output of the spread rate calculations if debugging:
    //JMR_NOTE: This may be overkill?????
    BOOST_LOG_SEV(glg, debug) << "Dump the spread calculations:";
    BOOST_LOG_SEV(glg, debug) << raData;

    //Simulate the combustion of surface fuels:----------------------------
    BurnupSim burnupOutput = SimulateSurfaceCombustion(fm, raData, tempAir, windSpeed);
    BOOST_LOG_SEV(glg, debug) << "Dump the combustion calculations:";
    BOOST_LOG_SEV(glg, debug) << burnupOutput;
    BOOST_LOG_SEV(glg, debug) << burnupOutput.history;

    siteFM = fm;//Save the input for use in getAbgVegetationBurntFractionsProcess().
    siteBU = burnupOutput;//Save the output for use in getAbgVegetationBurntFractionsProcess().

    if (FireBurned())
    {
      BOOST_LOG_SEV(glg, info) << "Surface fire ignited.";
  
      //Simulate crown fire:-------------------------------------------------
      std::vector <double> crownFireCHPAs = SimulateCrownFire();

      //Simulate ground fire:------------------------------------------------
      //The ground fire simulation uses the soil profile conditions (obtained by the function) and the
      //energy flux from the aboveground fire into the soil surface (RA + Burnup + crown).
      double surfFireHPA = burnupOutput.history.IntegrateFireIntensity();//Surface fire heat per area, kJ/m^2
      double abgFireHPA = surfFireHPA;
      double fireHeatFracToSoil = md.fire_heat_frac_to_soil;//For surface and passive crown fire.

      if (CFB > 0.0)
      {
        //The surface HPA is calculated by both Burnup and the crown fire model.  They should be
        //similar.  Report them so we can be check.
        BOOST_LOG_SEV(glg, debug) << "Surface fire HPA from Burnup: " << surfFireHPA << " kJ/m^2";;
        BOOST_LOG_SEV(glg, debug) << "Surface fire HPA from Scott & Reinhart 2021: "
                                  << crownFireCHPAs[0] << " kJ/m^2";;
        BOOST_LOG_SEV(glg, debug) << "Crown fire HPA from Scott & Reinhart 2021: "
                                  << crownFireCHPAs[1] << " kJ/m^2";;

        abgFireHPA += crownFireCHPAs[1];

        //For 'active' crown fire calculate the heat fraction into the soil:
        if (CFB >= 0.8)
        {
          //We estimated a relationship between crown fire intensity and heat into the soil based on
          //Thompson, Wotton, & Waddington 2015.  More support is needed.
          //We model the relationship above CFB = 0.8 as exponential decay with the parameters
          //a = 0.3753007, k = -1/4792.0683714, a maximum value of 0.35, and a minimum of 0.14.
          fireHeatFracToSoil = 0.3753007 * exp(-totalFireIntensity/4792.0683714) + 0.14;
          fireHeatFracToSoil = fmin(fireHeatFracToSoil, 0.35);
        }
        //For passive crown fire we use the surface fire parameter.
      }

      //Report fire heat with full numerical precission so it can be re-imported:
      BOOST_LOG_SEV(glg, debug) << "Aboveground fire heat per area (HPA): "
                                << std::setprecision(std::numeric_limits<double>::max_digits10)
                                << abgFireHPA << " kJ/m^2";
      //Only a fraction of the heat of the aboveground fire enters the soil:
      BOOST_LOG_SEV(glg, debug) << "Heat fraction to soil: " << fireHeatFracToSoil;
      double fireHeatToSoil = abgFireHPA * fireHeatFracToSoil;//kJ/m^2
      BOOST_LOG_SEV(glg, debug) << "Heat into soil: " << fireHeatToSoil << " kJ/m^2";

      burnDepth = SimulateGroundFire(fireHeatToSoil);
    }
    else
    {
      /*If the surface fire fails to ignite there is no reason to perform the the calculations for
      the subsequent stages.
      Note: Ground fire would be safe if there was no crown fire since it will produce no burn when
      the surface fire HPA is 0.  However, it is not clear if under the right conditions the Scott &
      Reinhart model could predict a crown fire even if Burnup doesn't predict a surface fire.
      Therefore it is safer to skip it.*/
      BOOST_LOG_SEV(glg, info) << "Surface fire failed to ignite.";
    }
  }

  //Litter and soil carbon stocks are updated in WildFire::updateBurntOrgSoil().
  //Moss and herbaceous carbon stocks are updated in WildFire::burnVegetation().

  //Match messaging of getBurnOrgSoilthick():
  BOOST_LOG_SEV(glg, debug) << "Setting the burn thickness in FirData...";
  fd->fire_soid.burnthick = burnDepth;

  BOOST_LOG_SEV(glg, info) << "Final Calculated Organic Burn Thickness: " << burnDepth;
  //This message is only here to have an invariant final line for searching the log:
  BOOST_LOG_SEV(glg, debug) << "Leaving WildFire::ProcessWildfire().";
  return burnDepth;
}

/** Determine the fuel model (number) matching a given community type:
 *
 * Each CMT has a predetermined fuel model assigned to it via the CMT fire parameter file.
 *
 * @param[in] cmt The number of the CMT for this location. [Currently ignored!]
 *
 * @returns The fuel model number (not the index) matching the current CMT.
 *
 * @note Used in the process wildfire model only.
 * 
 * @note The CMT to fuel model crosswalk is specified in the CMT fire parameter file.  Only the
 * values for the current site's CMT are loaded.  Therefore no lookup has to be done here and the
 * CMT, which is implied, could be omitted.  However, this may change in the future with a CMT
 * mapping to a differnt fuel model depening on other site factors like tree cover or drainage type.
 * A separate lookup table might be needed in that case.  We leave the CMT paramter for now.
 * The crosswalk is a work in progress with a placeholder value until fuel models are finalized.
 * We use TU1 = 161 since it is good for testing.  All fuel types are occupied and it is dynamic.
 */
int WildFire::GetMatchingFuelModel(const int cmt) const
{
  int fuelModelNumber;

  if (md.fire_temp_fm == -1)//Use the value supplied by the CMT parameter file:
  {
    fuelModelNumber = firpar.cmt2fm;
  }
  else//If provided override the fuel model with the value in the configuration file:
  {
    fuelModelNumber = md.fire_temp_fm;
    //We don't currently have a way to check that a fuel model number is valid.
  }
  //There is a bare land fuel model by number but it doesn't have parameters.  Do we need to provide
  //a bare land parameter set or can we just signal the calling code that it should skip fire
  //calculations?

  //If no match either throw an error or warn and return a default fuel model.

  return fuelModelNumber;
}

/** Determine fuel loadings based on the cohorts states and store them in the site's fuel model:
 *
 * This process is a theory of what represent fuels in DVM-DOS-TEM.  Once the mapping is defined
 * the fuel loadings are know, since model states are clearly defined.  Since the mapping itself is
 * a theory it represets a place where assumptions could change.  The current fuel mapping is:
 *
 * @par Dead fuels:
 * Dead surface fuels include litter and fallen debris.  In DVM-DOS-TEM these are represented by:
 * - Litter: The rawc component of the first non-moss soil layer.  rawc is part of the 'soil' but
 * is where all aboveground litter-fall ends up.  We assume this is not yet buried.
 * This may include some root cfall as well??????
 * We upscaled the total 'litter' carbon stock to biomass and distribute it among the fuel model's
 * dead fuel size classes.
 * - Large woody debris: Only after fire there is a pool of woody debris (wdebrisc) (and standing
 * dead stems?).
 * We could put all of that in the largest class but some of it might belong in the smaller classes
 * or in a 1000 hour class, which will have a limited effect on fire behavior and is missing from most
 * fuel models.
 * - Possibly live moss (see below).
 *
 * @par Live fuels:
 * Live fuels include graminoids, forbs, shrubs, and moss (see below).  We convert PFT aboveground
 * carbon stocks to biomass and sort into herbaceous (graminoids & forbs) and woody (shrubs) size
 * classes.
 *
 * @par Moss:
 * Moss is tricky.  It is a live fuel that can act like a dead fuel, with highly dynamic moisture.
 * The standard fuel models were developed with a focus on the lower 48 so the authors were not
 * thinking a lot about moss.  Moss is hoever an important driver of fire behavior in the north and
 * is ubiquitous in the DVM-DOS-TEM CMTs.  Figuring out how to handle moss is ongoing work.
 *
 * Given this uncertainty we are starting with the two simplest options, place it into either the
 * fine dead or live herbaceous type.  There are several other options.
 * I think it may make sense to treat it as a custom live fuel or potentially as a dynamic fuel.  It
 * should be given its own SAV if the live herbaceous is too coarse.  The model has a couple of
 * different moss PFTs and an SAV could be given for each.
 *
 * We don't include dead moss as a surface fuel.  We treat that as duff, part of the ground fuels.
 *
 * @param[in,out] fm The fuel model for the site (with default loadings).
 * @param[in] treatMossAsDead Should moss be treated as a fine dead fuel (true) or herbaceous live fuel (false)?
 *
 * @returns Nothing but the fuel model loadding (w_o_ij) is updated on return.
 *
 * @note Used in the process wildfire model only.
 *
 * ToDo:
 * - Record the mapping of stocks to fuels in some way record the value before fire so we can update
 * stocks appropriately afterwards.  We have done this to some extend in the calling code.
 * - Deal with wdebrisc.
 */
void WildFire::CohortStatesToFuelLoading(FuelModel& fm, const bool treatMossAsDead)
{
  BOOST_LOG_SEV(glg, debug) << "Entering WildFire::CohortStatesToFuelLoading()...";

  //Dead fuels:
  double rawC = GetLitterRawC();//Get the total litter carbon.

  //Distribute the rawc among the dead fuel size classes:

  //Get the dead fuel SAVs:
  int numDead = fm.NumDeadClasses();
  std::vector <double> savsDead(fm.SAV_ij.begin(), fm.SAV_ij.begin() + numDead);

  //Get an estimated distribution of dead fuel sizes:
  std::vector <double> distribSAVs, distribWts;//Empty vectors to hold the distribution.
  GetDeadFuelSizeDistribution(fm, distribSAVs, distribWts);

  //Distribute the litter carbon using the distribution:
  //All carbon stocks are in g/m^2 C and need to be converted to kg/m^2 dry biomass.
  std::vector <double> w_o_Dead = DistributeFuel(distribSAVs, distribWts, (rawC * c2b / gPerKg), savsDead);

  //Update the loadings (or do below):
  for (int i = 0; i < numDead; i++)
  {
    fm.w_o_ij[i] = w_o_Dead[i];
  }

  //Live fuels:
  //Standard fuel models have herbaceous and woody classes but in some models both are not actually
  //occupied.  Making sure the fuel classes match the CMT's PFTs is something that should be done
  //in advance but we do checking here for safety.

  //Sort the PFT biomass into the appropriate live fuel classes:
  
  //The fuel model may contain default loadings that should be disregarded:
  int liveHerbIndex = fm.LiveHerbaceousIndex();
  int liveWoodyIndex = fm.LiveWoodyIndex();
  fm.w_o_ij[liveHerbIndex] = 0;
  fm.w_o_ij[liveWoodyIndex] = 0;
  
  for (int pftIndex = 0; pftIndex < NUM_PFT; pftIndex++)
  {
    if (!cd->d_veg.ifwoody[pftIndex])//Or m_veg??????
    {
      //Put all herbaceous PFTs in the herbaceous:
      //Note: There is currently no way to tell graminoids and forbs appart.
      if (cd->d_veg.nonvascular[pftIndex] == 0)
      {
        //If the herbaceous class is not present adding carbon to it will not influence the fire
        //behavior.  If herbaceous fuels are not important in this system we could ignore them.
        //This is not something we can really resolve at run time.  This is a science question that
        //needs to be addressed during fuel model selection.
        if (!fm.LiveHerbaceousPresent())
        {
          BOOST_LOG_SEV(glg, fatal) << "The live herbaceous fuel type is not active in this fuel model.";
        }

        //Include aboveground parts:
        double leafC = bd[pftIndex]->m_vegs.c[I_leaf];
        double stemC = bd[pftIndex]->m_vegs.c[I_stem];
        fm.w_o_ij[liveHerbIndex] += (leafC + stemC) * c2b / gPerKg;//Convert to dry biomass.
      }
      else//Mosses:
      {
        //My best reading is that mosses and lichens are all leaf in DVM-DOS-TEM.  However if they
        //had stems and 'roots', i.e. rhizoids = root C, they should burn too.  Include all
        //compartments now in case for now:
        double mossBiomass = GetPFTBiomass(pftIndex);

        if (treatMossAsDead)
        {
          fm.w_o_ij[0] += mossBiomass;//Assumes fine fuel is first, which is pretty safe.
        }
        else
        {
          //See notes above.
          if (!fm.LiveHerbaceousPresent())
          {
            BOOST_LOG_SEV(glg, fatal) << "The live herbaceous fuel type is not active in this fuel model.";
          }

          fm.w_o_ij[liveHerbIndex] += mossBiomass;
        }
      }
    }
    else//Woody PFTs:
    {
      //Put shrubs in the live woody fuel:
      if (IsShrub(cd->cmttype, pftIndex))
      {
        //If the woody class is not present adding carbon to it will not influence the fire
        //behavior.  See notes for herbaceous fules above.
        if (!fm.LiveHerbaceousPresent())
        {
          BOOST_LOG_SEV(glg, fatal) << "The live woody fuel type is not active in this fuel model.";
        }

        //Include aboveground parts:
        double leafC = bd[pftIndex]->m_vegs.c[I_leaf];
        double stemC = bd[pftIndex]->m_vegs.c[I_stem];
        fm.w_o_ij[liveWoodyIndex] += (leafC + stemC) * c2b / gPerKg;//Convert to dry biomass.
      }
      //Ignore trees.
    }
  }
}

/** FW_MOD: Get the litter carbon for the site.
 *
 * The 'litter' is defined as the rawC componenet of the first non-moss soil layer.  See
 * CohortStatesToFuelLoading() for more information.
 * 
 * @note This started as a hack and could be moved into CohortStatesToFuelLoading().  However, it
 * does make the code a bit more organized.
 * @note Used in the process wildfire model only.
 * @note Now that we have added access to Ground this could be revised to use that.
 */
double WildFire::GetLitterRawC() const
{
  BOOST_LOG_SEV(glg, debug) << "Entering WildFire::GetLitterRawC()...";

  int topFibricIndex;

  BOOST_LOG_SEV(glg, debug) << "Getting litter carbon.";

  //Determine the top non-moss layer (fstshlwl in Ground terminology):
  //The site's Ground object is a private member of the parent cohort's Soil_Bgc object.  If we
  //could access this we wouldn't have to search the soil layers:
  //Layer* topFibric = thisCohort.soilbgc.ground->fstshlwl;//Get the top non-moss layer.
  //double rawC = topFibric->rawc;//Get the total 'litter' carbon mass.

  for (int i = 0; i < cd->m_soil.numsl; i++)//Assumes we are starting at the top layer.
  {
    if (cd->m_soil.type[i] == 1)//Shallow organic / peat ~ I_FIB
    {
      topFibricIndex = i;
      break;
    }
  }

  double topFibricRawC = bdall->m_sois.rawc[topFibricIndex];
  //Rarely this stock can be NaN.  Return zero instead:
  if (std::isnan(topFibricRawC))
  {
    return 0.0;
  }
  else
  {
    return topFibricRawC;
  }
}

/** Get the total biomass for a non-vascular PFT on the site.
 * 
 * @param[in] pftIndex The PFT index.
 *
 * @returns The total dry biomass (not carbon) for the specified PFT (kg/m^2).
 * 
 * @note: This code was moved out of CohortStatesToFuelLoading() to reduce code repetition.  It is
 *        only currently used for non-vascular PFTs but is there is nothing specific to those PFTs
 *        in it.
 * @note Used in the process wildfire model only.
 */
double WildFire::GetPFTBiomass(const int pftIndex) const
{
  double leafC = bd[pftIndex]->m_vegs.c[I_leaf];
  double stemC = bd[pftIndex]->m_vegs.c[I_stem];//Should be 0.
  double rootC = bd[pftIndex]->m_vegs.c[I_root];//Should be 0.
  double pftBiomass = (leafC + stemC + rootC) * c2b / gPerKg;//Convert to dry biomass.
  
  return pftBiomass;
}

/** Get the total biomass for all non-vascular PFTs on the site.
 *
 * @returns The total live dry biomass (not carbon) (kg/m^2).
 *
 * @note Used in the process wildfire model only.
 */
double WildFire::GetNonVascularBiomass() const
{
  double nvBiomass = 0.0;

  for (int pftIndex = 0; pftIndex < NUM_PFT; pftIndex++)
  {
    if (cd->d_veg.nonvascular[pftIndex] != 0)
    {
      nvBiomass += GetPFTBiomass(pftIndex);
    }
  }

  return nvBiomass;
}

/** Get the estimated dead fuel size distribution for the site.
 * The distribution is returned via the vectors passed.
 *
 * This is a temporary stub-ish implementation of the code.  For now we assume fm is a standard fuel
 * model and that the default loadings represent a typical dead fuel size distribution.  This is
 * probably not the case and a better method for estimating the size distribution will be developed
 * in the future using literature or calculations.  Then this function will be updated to some
 * lookup process using the CMT of fuel model number to get the estimate that is calculated or
 * supplied via an input file.
 * 
 * @param[in] fm A fuel model to base the distribution on. [Temporary!!!!!]
 * @param[out] distribSAVs An empty vector to return the SAVs of the distribution in.
 * @param[out] distribWts An empty vector to return the weights of the distribution in.
 * 
 * @returns Nothing.  Parameters are updated instead.
 *
 * @note Used in the process wildfire model only.
 */ 
void GetDeadFuelSizeDistribution(const FuelModel& fm, std::vector <double>& distribSAVs,
                                 std::vector <double>& distribWts)
{
  BOOST_LOG_SEV(glg, debug) << "Entering GetDeadFuelSizeDistribution()...";

  //Copy the appropriate dead members:
  int numDead = fm.NumDeadClasses();
  distribSAVs.assign(fm.SAV_ij.begin(), fm.SAV_ij.begin() + numDead);
  distribWts.assign(fm.w_o_ij.begin(), fm.w_o_ij.begin() + numDead);
}

/** Is this PFT a shrub?
 *
 * Currently there is no shrub or tree flag so we have to use other information to infer
 * shrubbiness.  This is a work in progress and may be replaced by a flag later.
 *
 * The vegetation parameter files include names for most of the PFTS, but only in comments, which
 * are not imported.  These include DecShrub and EvrShrub.  The order of PFTs in a CMT are also not
 * systematic.  The only approach that remains is to check the PFT numbers against their CMTs.  This
 * requires knowledge of the CMT file (cmt_bgcvegetation.txt).  Since the CMTs will change in due
 * course this implementation is not robust and should be replaced soon.
 *
 * Since plants can have different growth forms and dwarf statures, especially in harsh climates, we
 * could try to used stature to help infer when trees are shrubby in a way that is relevant to fuel
 * models.  That might be getting a bit fancy.
 *
 * @param[in] cmtNumber The CMT number.
 * @param[in] pftIdx The index of the PFT to check.
 *
 * @returns True if this PFT is a shrub. 
 */
bool IsShrub(const int cmtNumber, const int pftIdx)
{
  BOOST_LOG_SEV(glg, debug) << "Entering IsShrub()...";

  //I believe all the following PFTs are considered shrubs:
  //A few of these CMT cases could be combined but I'm prioritizing readability over compactness.
  switch (cmtNumber) {
    case 0://Bare ground, contains junk values.
      break;

    case 1://Boreal Black Spruce
    case 2://Boreal White Spruce Forest
    case 3://Boreal Deciduous Forest
      if (pftIdx == 1)//DecidShrub
      {
        return true;
      }
      break;

    case 4://Shrub Tundra (Toolik area)
      if (pftIdx >= 0 && pftIdx <= 3)//Salix, Betula, Decid, EGreen
      {
        return true;
      }
      break;

    case 5://Tussock Tundra
      if (pftIdx >= 0 && pftIdx <= 2)//Betula, Decid, EGreen
      {
        return true;
      }
      break;

    case 6://Wet Sedge Tundra (Toolik area)
      if (pftIdx == 0)//Decid
      {
        return true;
      }
      break;

    case 7://Heath Tundra
      if (pftIdx >= 0 && pftIdx <= 1)//Decid, EGreen
      {
        return true;
      }
      break;

    case 12://Lowland Boreal Wet Shrubland
      if (pftIdx >= 0 && pftIdx <= 1)//DecShrub, EvrShrub
      {
        return true;
      }
      break;

    case 20://EML Shrub Tundra
      if (pftIdx == 0)//Betnan (Betuala nana)
      {
        return true;
      }
      break;

    case 21://EML Tussock Tundra (2022)
      if (pftIdx >= 0 && pftIdx <= 1)//Decidsh, Egreensh 
      {
        return true;
      }
      break;

    case 31://Boreal Bog (Based on PFT order in all parameter files except cmt_firepar.txt, which differs)
      if (pftIdx >= 2 && pftIdx <= 3)//DecShrub, EvrShrub
      {
        return true;
      }
      break;

    case 44://Shrub Tundra (Kougarok)
      if (pftIdx >= 0 && pftIdx <= 3)//Salix, Betula, Decid, EGreen
      {
        return true;
      }
      break;

    case 50://Shrub Tundra TVC
      if (pftIdx >= 0 && pftIdx <= 1)//DecidShrub, EgreenShrub
      {
        return true;
      }
      break;

    case 51://Tussock Tundra TVC
    case 52://Heath-lichen (Trail Valley)
      return false;
      break;

    case 55://Arctic Fen (Daring Creek)
      if (pftIdx == 1)//EgreenShrub
      {
        return true;
      }
      break;

    case 56://Mixed Tundra (Darling Lake)
      if (pftIdx >= 2 && pftIdx <= 3)//EgreenShrub, DecidShrub
      {
        return true;
      }
      break;

    case 57://Heath Tundra (Daring Lake)
      if (pftIdx == 2 || pftIdx == 4)//EgreenShrub, DecidShrub
      {
        return true;
      }
      break;

    case 60://Peat Plateau
      if (pftIdx == 1)//EricShrub
      {
        return true;
      }
      break;

    case 65://Mature Aspen (BOREAS)
    case 66://Mature Jack Pine (BOREAS)
      if (pftIdx == 1)//DecidShrub
      {
        return true;
      }
      break;

    case 67://Mixed Forest (Groundhog Rover)
      return false;
      break;

    case 69://Black Spruce (BOREAS)
      if (pftIdx == 2)//Shrub
      {
        return true;
      }
      break;

    case 70://Dwarf Shrub (Chokurdakh)
      if (pftIdx >= 0 && pftIdx <= 1)//DecidShrub, EgreenShrub
      {
        return true;
      }
      break;

    case 71://Larch Forest (Yakutsk Spasskaya Pad)
      if (pftIdx == 2)//Shrubs
      {
        return true;
      }
      break;

    case 73://Tussock Tundra (Cherskii)
      if (pftIdx == 1)//Shrub
      {
        return true;
      }
      break;

    case 74://Scots Pine (Zotino)
      if (pftIdx == 2)//DwarfShrub
      {
        return true;
      }
      break;

    case 75://Bog-Igarka
      if (pftIdx == 0)//Shrubs
      {
        return true;
      }
      break;

    case 76://Shrub Tundra-Seida
      if (pftIdx >= 0 && pftIdx <= 2)//DecidShrub, Betnan, Salix
      {
        return true;
      }
      break;

    case 82://ScotsPine-Fennosc
      if (pftIdx == 2)//EGreenShrub
      {
        return true;
      }
      break;

    case 90://Heath-Zachenberg
      if (pftIdx == 0)//Dwarfshrub
      {
        return true;
      }
      break;

    case 92://Abisko Bog - bonanza calibration
      if (pftIdx == 1 || pftIdx == 3)//Eshrub, Dshrub
      {
        return true;
      }
      break;

    case 93://Adventalen - imnavait calibration
      if (pftIdx == 0)//Decid
      {
        return true;
      }
      break;

    default:
      BOOST_LOG_SEV(glg, fatal) << "IsShrub() does not know this CMT: " << cmtNumber;
      break;
  }
  return false;
}

/** Is this PFT a tree (and what kind)?
 *
 * Currently there is no tree flag so we have to use other information to infer a PFT is a tree.
 * This may be replaced by a flag or wrap a flag later.  i.e:
 * cd->d_veg.TBD_TREE_FLAG[pftIndex] == ?????)
 *
 * This function has dependacies on a particular version of the parameter file via IsShrub().  
 *
 * @param[in] cmtNumber The CMT number.
 * @param[in] pftIdx The index of the PFT to check.
 *
 * @returns Postive if a tree (1 = conifer, 2 = decidous), 0 if not a tree.
 */
int WildFire::IsTree(const int cmtNumber, const int pftIdx) const
{
  int treeCode = 0;

  if (cd->d_veg.ifwoody[pftIdx] && !IsShrub(cd->cmttype, pftIdx))
  {
    if (!cd->d_veg.ifdeciwoody[pftIdx])
    {
      treeCode = 1;//Conifer
    }
    else
    {
      treeCode = 2;//Decidous
    }
    //Larch might deserve it's own class.
  }

  return treeCode;
}

/** Calculate the fuel bed depth and update it in the fuel model passed.
 *
 * Standard fuel models include a fuel bed depth among their default values.  We need a way to
 * estimate this from the model along with fuel loadings.  DVM-DOS-TEM doesn't give us much to work
 * with.  Ecosystem structural information would be most useful but we don't have that.  For now we
 * will make some assumptions.  Here are some options:
 * 
 * - The fuel bed depth is constant.  This may be reasonable for gramaniod dominated ecosystems,
 * e.g. grasslands, tussock tundra, etc., where the vegetation has a typical maximum height but may
 * be more or less dense.
 * 
 * - The mix of fuels and it structure is fairly constant but the amount is not.  That implies the
 * bulk density of the fuel bed is constant but the height varies with fuel loading.  Fuel bed depth
 * is then calculable from fuel loading and an estimate of they typical bulk density.
 *
 * The reality is likely more complicated but we will build our approach around these for now.
 *
 * @param[in,out] fm The fuel model for the site (with default fuel bed depth and updated fuel loadings).
 * @param[in] dynamic Temporary: Should the fuel bed depth be treated as changing with fuel amounts?
 *
 * @returns Nothing.  The fuel bed depth is updated in the fuel model passed.
 *
 * @note Used in the process wildfire model only.
 */
void CalculateFuelBedDepth(FuelModel& fm, const bool dynamic)
{
  BOOST_LOG_SEV(glg, debug) << "Entering CalculateFuelBedDepth()...";

  //Determine which approach to use based on the fuel model or CMT would be ideal but it will
  //require some research.  For now we use a switch.  If the depth is constant there is nothing to
  //do.  Otherwise calculate the depth based on constant bulk density:
  if (dynamic)
  {
    double totalLoading = std::accumulate(fm.w_o_ij.begin(), fm.w_o_ij.end(), 0.0);
    //This assumes fm.bulkDensity is the standard fuel model default value:
    fm.delta = totalLoading / fm.bulkDensity;//kg/m^2 / kg/m^3 = m
  }
}

/** Get the air temperature at the time of fire.
 *
 * For calculating fire behavior the temperature at the time of fire is most important.  We don't
 * know when in the fire occurs but we would like to calcualte midday to afternoon activity.  This
 * function returns the daily temperature, which is itself interpoated from monthly values.  In the
 * future we may use the min, mean, and max daily values to calculate a better approximation of the
 * temperature at a time of day.  To provide a bit more control for now allow the temperature to be
 * overiddin with a value in the configuration file.  This is limited as it will be applied to all
 * fires across all dates and locations.
 *
 * @returns The air temperature at the time of fire (C).
 */
double WildFire::GetAirTemperature() const
{
  double tempAir;

  if (md.fire_tempair == -300.0)//Use the value from the input stream:
  {
    tempAir = edall->d_atms.ta;//Daily air temp (at surface).
    //Would it be better to get the temperature from climate instead?
    //float tempAir = climate->tair[monthIndex]//Monthly
    //float tempAir = climate->tair_d[dayOfYearIndex];
  }
  else//If provided override with the temperature from the configuration file:
  {
    BOOST_LOG_SEV(glg, info) << "Air temperature obtained from config file.";
    tempAir = md.fire_tempair;
  }

  return tempAir;
}

/** Get the wind speed at midflame height in m/min.
 *
 * The Rothermel Albini spread model takes midflame wind speed, an ill defined quantity representing
 * the mean wind speed in the flaming zone.  An estimate of near surface wind speed is what we need.
 * 2 meter wind speed should be fine.  The input may not match this so we need to compute our best
 * estimate.
 *
 * We need a sub-daily value ~ daily value?
 * The code has to do the following:
 * - Get the daily windspeed from the host model.
 * - Calculate the windspeed at ~2m if the provided wind speed is at another height.
 * - Possibly: Estimate a sub-daily from a daily mean value.  The fire may occur on timescale where
 *   it occurs at a specific time of day.  In this case it would be good to estimate what the wind
 *   speed was at this time.  This may be difficult as wind can be highly variable and may not be
 *   well predicted by diurnal cycles.
 * - Convert to m/min if needed.
 *
 * @returns The wind speed at 2 meters height (m/min).
 *
 * @note: We could also pass in the desired height.  We currently need ~2m and will have 2m wind
 * input so height adjustment is not necessary.  Passing in the time of day would allow us to adjust
 * daily to sub-daily wind (see above).
 *
 * @note Used in the process wildfire model only.
 * @note Windspeed is not yet available in DVM-DOS-TEM, but it should be soon.  This code returns
 *       a stub value that may be overriden in the configuraiton file.
 */
double WildFire::GetMidflameWindSpeed() const
{
  double windSpeed;//Return value.
  
  if (md.fire_windspeed == -1.0)//Use the value from the input stream:
  {
    //Draft:
    //If this is daily how do we know what day of the month it is?
    //The wind speed will be provided as directional components in m/s at 2 meters.
    //double u = edall.d_atms.EasternWindSpeed;//Zonal component U.		Confirm units!!!!!
    //double v = edall.d_atms.NorthernWindSpeed;//Meridional component V.
    //windSpeed = std::sqrt(std::pow(u, 2) + std:pow(v, 2));//Get the vector wind speed.
    //windSpeed *= 60;//Convert from m/s to m/min.
    //This a daily average value.  An afternoon value would probably be better.

    //Temporary stub, return an arbitrary value!!!!!:
    //Summer average wind speeds are ~6 mph in Fairbanks Alaska.
    //6 * ftPerMi / 60 / ftPerM = 160.9344
    BOOST_LOG_SEV(glg, debug) << "Stub value for wind speed used.";
    windSpeed = 160.9344;
  }
  else//If provided override with the windspeed from the configuration file:
  {
    BOOST_LOG_SEV(glg, info) << "Wind speed obtained from config file.";
    windSpeed = md.fire_windspeed;
  }

  return windSpeed;
}

/** Get the wind reduction factor for the site.
 *
 * The wind reduction factor is a function of forest stand density and structure and commonly
 * ranges from 0.1 for a dense fully sheltered stand to 0.3 for sparce partially sheltered stands.
 * 0.4 would be normally be considered unsheltered but is also used as a special value in the crown
 * fire calculations of Rothermel 1991.
 *
 * This is currenltly a stub.  A calculation or input parameters may be added.
 *
 * @returns WRF, the wind reduction factor used to estimate midflame wind speed (U) from open
 *          (6.1 m) windspeed (O), (fraction, midflame wind speed / open wind speed).
 *
 * @note Used in the process wildfire model only.
 * @note This currently only intended for use with forested CMTs with crowning potential.
 */
double WildFire::GetWindReductionFactor() const
{
  BOOST_LOG_SEV(glg, debug) << "Stub value for wind reduction factor used.";
  return 0.15;
}

/** Get the relative humidity at the time of fire.
 *
 * @param dayOfYearIndex The day of the year to get the humidity for.
 *
 * @returns The percent relative humidity at the time of fire (%).
 */
double WildFire::GetRelativeHumidity(const int dayOfYearIndex) const
{
  double rhPct;

  if (md.fire_rh_pct == -1.0)//Use the value from the input stream:
  {
    //Atmospheric pressure is needed for the the Fireweed code method.
    //Sea level barometric pressure will be available soon but is not currently.  Use this along with
    //cd->cell_elevation() to get the pressure at this elevation.
    //p_hPa = Z;

    //Partial pressure of water vapor:
    //This is an input variable.  There has been some discussion recently about the units of this.
    //The docs say it is in hPa but it is used in the code as if it is in Pa
    //(see Climate.cpp calculate_vpd()).
    float P_Pa = climate->vapo_d[dayOfYearIndex];//My reading is that vapo_d is a vector of daily values for the whole year.

    //Saturated vapor pressure:
    //double P_s_hPa = SaturationVaporPressureBuck(tempAir, p_hPa);//Fireweed method returns hPa.
    //This is a computed climate variable.
    float P_s_Pa = climate->svp_d[dayOfYearIndex];//From the code this must be in Pa.

    //Use the pressures to calculate relative humidity:
    rhPct = RHfromVP(P_Pa, P_s_Pa);
  }
  else//If provided override with the windspeed from the configuration file:
  {
    BOOST_LOG_SEV(glg, info) << "Relative humidity obtained from config file.";
    rhPct = md.fire_rh_pct;
  }

  return rhPct;
}

/** Calculate fuel moisture based on recent weather:
 *
 * An alternative to calculating this at the time of fire would be to make fuel moisture a
 * continuously calculated state.  This would moke more sense if litter existed as a distinct stock
 * as well.
 *
 * @param[in] fm The fuel model object for this site.
 * @param[in] monthIndex The current month as a zero based index.
 #
 * @returns M_f_ij, the fuel moisture for all fuel classes.  This is not returned in the fuel model
 * passed in because we don't know if curing is being applied.
 *
 * @note Used in the process wildfire model only.
 */
std::vector <double> WildFire::CalculateFuelMoisture(const FuelModel& fm, const int monthIndex) const
{
  BOOST_LOG_SEV(glg, debug) << "Entering WildFire::CalculateFuelMoisture()...";

  //We get the number of fuel classes here and then assume the number below.  To handle more than
  //the standard five fuel types we need additional tools to undestand what they are.
  std::vector <double> M_f_ij(fm.numClasses, 0);//Return value.

  //Get the current date:
  //There is a month member in the cohorts CohortData but that apparently is not the current month.
  int monthOfYear = monthIndex + 1;

  //The day of month is really up to us.  The middle of the month seems resonable.  We'll use the
  //15th for now and add a calculation later?
  int dayOfMonthIndex = 14;//0 based.

  //The daily data provided in Climate is a vector and should be 0 indexed but the code notes that
  //it currently returns 366 days.
  int dayOfYearIndex = temutil::day_of_year(monthIndex, dayOfMonthIndex);

  //Dead fuel moisture:---------------------------
  BOOST_LOG_SEV(glg, debug) << "Calculating dead fuel moisture:";//Temporary?????

  //Current air temperature:
  double tempAir = GetAirTemperature();

  //Calculate current relative humidity:
  double rhPct = GetRelativeHumidity(dayOfYearIndex);
  BOOST_LOG_SEV(glg, debug) << "Percent relative humidity: " << rhPct;//How many digits do we need?

  //We assume that fires occur in the mid-afternoon.  The exact hour might need to be shared across
  //the entire fire code.  The Fosberg model was originally devise with afternoon (~14:30) values in
  //mind, but the table method used below can take any daylight time.
  int hourOfDay = 15;

  //The percent slope is stored in he CohortData object and also in the wildfire object:
  //Note: Slope is also duplicated in the Wildfire object.
  double slopePct = cd->cell_slope;//Percent
  BOOST_LOG_SEV(glg, debug) << "Percent slope: " << slopePct;

  //Get the aspect:
  double aspect = cd->cell_aspect;//Degrees

  //Shading should take into consideration both cloudiness and canopy cover:
  bool shaded = false;

  //Canopy cover is equivalent to the projected foliage area.  For the Fosberg NWCG method this is
  //estimated visually.  As such it best to consider this the as the overhead tree canopy.  We can't
  //get this exactly but we can approximate it by excluding the herbaceous and moss leaf area:
  double fpcWoody = 0;
  for (int i = 0; i < NUM_PFT; i++)
  {
    if (cd->d_veg.ifwoody[i])
    {
      fpcWoody += cd->d_veg.fpc[i];
    }
  }

  //To combine the cloudiness (expressed as a percentage) and canopy cover (a fraction) we need to
  //invert the the values so the product increases the effect when combined.  We then re-invert:
  double shadeFrac = 1 - ((1 - (climate->cld_d[dayOfYearIndex] / 100)) * (1 - fpcWoody));

  if (shadeFrac > 0.5)
  {
  	shaded = true;
  }

  //Perform the 1hr moisture look up and derive the rest from that:
  BOOST_LOG_SEV(glg, debug) << "Calling FosbergNWCG_1HrFM()";//Temporary?????
  double oneHrFM = FosbergNWCG_1HrFM(md.fire_fosberg_a_file, md.fire_fosberg_b_file,
                                     md.fire_fosberg_c_file, md.fire_fosberg_d_file,
                                     tempAir, rhPct, monthOfYear, hourOfDay,
                                     slopePct, aspect, shaded);//Default values for the rest.

  double tenHrFM = NWCG_10hrFM(oneHrFM);
  double hundredHrFM = NWCG_100hrFM(oneHrFM);

  //Live fuel moisture:---------------------------
  //The GSI based fuel moistures should be averages of the 21 days up to today:
  //Monthly values might work but I suspect minimum daily temp could cause big departures at some
  //times of the year.
  BOOST_LOG_SEV(glg, debug) << "Calculating live fuel moisture:";//Temporary?????

  double herbLFM = 0;
  double woodyLFM = 0;
  for (int i = dayOfYearIndex - 20; i <= dayOfYearIndex; i++)
  {
    //Minimum and maximum daily temperature are not currently available but will be soon.  We
    //use a hack value for now:
    float tempCMin = climate->tair_d[i] - 10;
    
    //Get the VPD:
    float vpdPa = climate->vpd_d[i];
    //Or calculate VPD from mean daily weather values:
    //double p_hPa = ?????;
    //...
    //double vpd_hPa = VPDfromRHBuck(tempAir, rh, p_hPa);

    //temutil::length_of_day gives the day length in hours:
    float dayLengthSec = temutil::length_of_day(lat, dayOfYearIndex) * 60 * 60;

    double gsi = GrowingSeasonIndex(tempCMin, vpdPa, dayLengthSec);
    herbLFM += HerbaceousLiveFuelMoisture(gsi);
    woodyLFM += WoodyLiveFuelMoisture(gsi);
  }
  //Average:
  herbLFM = herbLFM / 21;
  woodyLFM = woodyLFM / 21;


  //Combine the live and dead moisture:-----------
  BOOST_LOG_SEV(glg, debug) << "Setting M_f_ij:";//Temporary?????
  //This makes assumption that the order is that of a standard fuel model.
  //It would be better to inform the numbers using information from the fuel model.  See above.
  //The units need to change from percent moisture to moisture fraction.
  M_f_ij[0] = oneHrFM / 100;
  M_f_ij[1] = tenHrFM / 100;
  M_f_ij[2] = hundredHrFM / 100;
  M_f_ij[3] = herbLFM / 100;
  M_f_ij[4] = woodyLFM / 100;

  return M_f_ij;
}

/** Calculate conifer foliar moisture content.
 *
 * The moisture cotent of conifir needles typcally varies from around 75 - 150%.  It varies by tree
 & species, season, location, elevation, and needle age.  It is not stongly responsive to
 * envioromental conditions.
 * 
 * This function is little more than a stub that returns the a typical value of 100% recommened in
 * the absence of other infromation (Scott & Reinhardt 2001).  It may be improved in the future.
 * The month of year and species information are candidates that could be used to better inform
 * values.
 *
 * Could add @param[in] monthIndex The current month as a zero based index.
 #
 * @returns FMC Foliar moisture content of (conifer) canopy (%, water weight/dry fuel weight x 100).
 *
 * @note Used in the process wildfire model only.
 */
double WildFire::CalculateFoliarMoistureContent() const
{
	BOOST_LOG_SEV(glg, debug) << "Stub value for foliar moisture content used.";
	return 100.0;
}

/** Get the canopy bulk density for the stand.
*
* This returns a value but there is a bit of work to make it robust.
*
* @returns CBD Canopy bulk density, the dry mass of canopy fuel, primarily foliage (needles) and
*              fine branches, per volume of the canopy (kg/m^3).  AKA crown bulk density.
*
* @note Used in the process wildfire model only.
*/
double WildFire::GetCanopyBulkDensty() const
{
  double CBD = 0.0;

  if (firpar.cbd > 0.0)
  {
    CBD = firpar.cbd;
    //It the value is passed in we should do some sanity checks on it.
    //Can we use the canopy fuel / m^2 to see if it needs to be reduced?
    //I haven't come with a way to do the check yet.
  }
  else
  {
    //To compute CBD we need fuel mass and canopy volume.  We can get foliage mass per area and
    //estimate a corresponding fine branch (an bark) biomass from that.  We don't have a way to know
    //the canopy depth in the model at this time.

    BOOST_LOG_SEV(glg, debug) << "Stub value for canopy bulk density used.";
    CBD = 0.15;//An initial stub value for when an input it missing.
  }

  return CBD;
}

/** Get the crown base height for the stand.
*
* @returns CBH Crown base height (m). AKA canopy base height, live crown base height (LCBH), z in
*              original Van Wagner notation.
*
* @note Used in the process wildfire model only.
* @note This currently just returns the parameter value (or a stub) but we may add age adjustment
*       etc. to the output.
* @note CBH is only valid for forests and we could check for that.
*/
double WildFire::GetCrownBaseHeight() const
{
  double CBH = 0.0;
  
  if (firpar.cbh > 0.0)
  {
    CBH = firpar.cbh;
  }
  else
  {
    //In the future a missing value may be treated as an error.
    //For now provide a stub value to allow for testing while values are compiled.
    BOOST_LOG_SEV(glg, debug) << "Stub value for crown base height used.";
    CBH = 5.0;
  }
  
  return CBH;
}

/* Get the canopy fuel load for the stand.
*
* Canopy fuel consist of foliage (needles), some degree of fine branches, as well as some bark
* lichen etc. in some cases.  Crown fire models focus on coniferous crowns so this function returns
* an estimate of coniferous fuel.
*
* @returns CFL, canopy fuel load per area (kg/m^2).  AKA crown fuel load.
*/
double WildFire::GetCanopyFuelLoad() const
{
  double CFL = 0.0;
  /*Estimate the non-foliage fuel as a ratio to the foliage (e.g. 0.5 means 5 g of branches etc.
  burn for every 10 of foliage.  Fine branches directly supply needles so this ratio is an
  expression of a pipe model allometry.  Other crown fuel like bark may not scale so directly.  It
  is also the case that the ratio probalbly increases with fire intensity but that is beyond what we
  have data to support at this time.*/
  const double nonFoliageFuelRatio = 0.5;//Uniformed intial guestimate.  Move to parameter?

  //Loop through the CMT's PFTs.
  for (int pftIndex = 0; pftIndex < NUM_PFT; pftIndex++)
  {
    //Add up canopy fuel for all coniferous tree PFTs:
    if (IsTree(cd->cmttype, pftIndex) == 1)
    {
      //Foilage and non-foilage carbon converted to dry biomass:
      CFL += bd[pftIndex]->m_vegs.c[I_leaf] * nonFoliageFuelRatio * c2b / gPerKg;
    }
    //We assume that decidous leaves may burn but due to their moisture they are heat sinks to
    //neutral at best and don't contrbute to crown fire behavior.
  }

  return CFL;
}

/** Simulate combustion of surface fuels.
 *
 * We simulate surface fuel combustion with the Burnup model of Albini 1995.
 *
 * @param[in] fm The fuel model for the site.
 * @param[in] raData The results of a Rothermel Albini spread model calculation for the site.
 * @param[in] tempAir The ambient air temperature (C).
 * @param[in] windSpeed The wind speed at 2 meters height (m/min).
 *
 * @returns An object containing the simulation results.
 *
 * @note Used in the process wildfire model only.
 */
BurnupSim SimulateSurfaceCombustion(const FuelModel& fm, const SpreadCalcs raData,
                                    const double tempAir, const double windSpeed)//Will be const when added to WildFire.
{
  BOOST_LOG_SEV(glg, debug) << "Entering SimulateSurfaceCombustion()...";

  //Burnup takes a number of parameters:
  //Wind and air temperature are needed (passed in).

  //Fuel properties are a superset of those in a standard fire behavior fuel model:
  //Leave all remaining properties at their default values set in BurnupFM().

  //Duff loading and moisture:
  /*It is unclear exactly what should be used here.  The duff concept in Burnup is certainly not
  synonymous with the deep duff of the boreal.  The shallow duff / live moss layer is probably the
  appropriate analog.  We handle that as a surface fuel so we do not include a duff loading to avoid
  double counting.  We treat smoldering ground fire as a separate stage too.  Duff in Burnup treated
  very simply as an additional source of heat rather than a full blown fuel.  Omitting duff seems
  reasonable but we should revisit this assumption in the future.*/
  double duffLoading = 0.0;
  double duffMoisture = 0.0;

  //Igniting fire properties are needed from the surface fire:
  //The model takes fire intensity and residence time.  These can be estimated from the spread rate
  //calculations and associated functions.

  //Fire intensity of the flaming front must be estimated.  The appropriate calculation to use is
  //the question since Rothermel & Albini and Burnup terminology are not identical.
  //The reaction intensity (I_R) represents the total energy generated by the fire per unit area and
  //time (kJ/m^2/min).  The heat source uses the propagating flux ratio, wind, and slope to give the
  //energy that heats fuel ahead of the fire from I_R.
  double fireIntensity = raData.heatSource / 60.0;//Convert kJ/m^2/min -> kW/m^2
  //kJ/m^2/min = kJ/min/m^2 = kJ/(60 * s)/m^2 =  kJ/s / 60 /m^2 = kW / 60 /m^2 = (kW/m^2) / 60

  double t_r = ResidenceTime(raData.cSAV, Metric) * 60.0;//Convert minutes -> seconds.  Or fm.cSAV

  //Simulation settings:
  //We can experiment with these.  Leave other at their default values.
  const double dT = 1.0;//1 second produces simulations with consistant burnout times. WILDFIRE_PARAMETER
  const int nTimeSteps = 45000;//This is 12.5 hours which may be excessive. WILDFIRE_PARAMETER

  //Call Burnup:
  //The simulation calculates the both consumption of fuels and the time evolution of the fire
  //behavior.  The output currently only contains fire history limited to fire intensity.
  BurnupSim output = BurnupFM(fm, duffLoading, duffMoisture, tempAir, windSpeed, fireIntensity,
                              t_r, dT, nTimeSteps);

  return output;
}

/** Did a wildfire ignite and consume fuels.
 *
 * With the classic wildfire model an ignition event specified via the input files will result in
 * a fire, that is fire effect will be calculated.  However, with the process wildfire model a fire
 * may not 'ignite' if condtion prevent it.  This function can be called by outside code following
 * should_ignite() to know if a fire occured and downstream fire effects need to be calculated.
 *
 * @returns Did a wildfire actually burn.
 *
 * @note Added to extend functionality for the process wildfire model.
 * @note Place after WildFire::should_ignite().
 */
bool WildFire::FireBurned()//FireIgnited()?
{
  //We determine whether a fire was simulated to have happen based on the model state.  The result
  //will only be valid if this function is called internally after WildFire::ProcessWildfire() and
  //externally /after WildFire::burn() has been called.
  if (!md.fire_process_wildfire)
  {
    return true;//This assumes that should_ignite() already returned true.
  }
  else
  {
    //There are several ways a process based wildfire can fail.  See Burnup for more information.
    return (siteBU.burnoutTime > 0.0);
  }
}

//Place after WildFire::getBurnAbgVegetation().
/** Calculate the fate of aboveground vegetation after fire using the process based wildfire model.
 *
 * This function determines the fractions of a given PFT that are combusted, are killed but are not
 * combusted, and that survive fire based on the process based wildfire model.
 *
 * The fractions are returned via class data members.  (See notes in getBurnAbgVegetation)
 *
 * @param[in] pftIndex The PFT index.
 *
 * @returns Nothing.  Class members are updated.
 *
 * @note Used in the process wildfire model only.
 * @note This function needs access to the current fuel model and data output from Burnup.  For now
 *       those are stored in private data members.
 */
void WildFire::getAbgVegetationBurntFractionsProcess(const int pftIndex)//Name is a bit long!
{
  //Reset rates to zero::
  this->r_burn2ag_cn = 0;
  this->r_dead2ag_cn = 0;

  //If for any reason fire didn't actually burn skip most of the calculations:
  if (siteBU.burnoutTime > 0.0)
  {
    //Fuel types may contain one or more PFTs but since we are calculating rates we don't need to
    //separate them out.  The rate for all will apply to each.  Likewise the carbon to biomass
    //conversion cancels out so can be ignored.

    //Get fuel indexes:
    int liveHerbIndex = siteFM.LiveHerbaceousIndex();
    int liveWoodyIndex = siteFM.LiveWoodyIndex();
    int deadHerbIndex = siteFM.DeadHerbaceousIndex();

    //Parallel the PFT decision tree in CohortStatesToFuelLoading():
    if (!cd->d_veg.ifwoody[pftIndex])
    {
      if (cd->d_veg.nonvascular[pftIndex] == 0)//Herbaceous PFTs:
      {
        //Herbs map to one fuel or two if dynamic moisture is on:
        if (siteFM.type != Dynamic)
        {
          if (siteBU.w_o_ij_Initial[liveHerbIndex] > 0.0)
          {
            this->r_burn2ag_cn = siteBU.combustion_ij[liveHerbIndex] / siteBU.w_o_ij_Initial[liveHerbIndex];
          }
        }
        else//Dynamic fuel moisture:
        {
          //Dynamic live moss maps to the live herbaceous fuel and cured herbaceous fuel:
          double totalInitialLoading = siteBU.w_o_ij_Initial[deadHerbIndex] + siteBU.w_o_ij_Initial[liveHerbIndex];
          if (totalInitialLoading > 0.0)
          {
              this->r_burn2ag_cn = (siteBU.combustion_ij[deadHerbIndex] +
                                    siteBU.combustion_ij[liveHerbIndex]) / totalInitialLoading;
          }
        }
      }
      else//Mosses:
      {
        //Moss: Maps to either one (fine) dead fuel OR one live (herbaceous) fuel OR a mix of both:
        if (md.fire_moss_as_dead_fuel)//Treat moss as a dead fuel:
        {
          int deadMossIndex = 0;//Assumes moss is in the finest dead fuel class.
          if (siteBU.w_o_ij_Initial[deadMossIndex] > 0.0)
          {
            this->r_burn2ag_cn = siteBU.combustion_ij[deadMossIndex] / siteBU.w_o_ij_Initial[deadMossIndex];
          }
        }
        else//Treat moss as a live fuel:
        {
          //Considered using aliases:
          //int deadMossIndex = deadHerbIndex;
          //int liveMossIndex = liveHerbIndex;
          if (siteFM.type != Dynamic)
          {
            //Live moss maps to the live herbaceous fuel:
            //Do we also need to check LiveHerbaceousPresent()?
            if (siteBU.w_o_ij_Initial[liveHerbIndex] > 0.0)
            {
              this->r_burn2ag_cn = siteBU.combustion_ij[liveHerbIndex] / siteBU.w_o_ij_Initial[liveHerbIndex];
            }
          }
          else//Dynamic fuel moisture:
          {
            //Dynamic live moss maps to the live herbaceous fuel and cured herbaceous fuel:
            double totalInitialLoading = siteBU.w_o_ij_Initial[deadHerbIndex] + siteBU.w_o_ij_Initial[liveHerbIndex];
            if (totalInitialLoading > 0.0)
            {
              this->r_burn2ag_cn = (siteBU.combustion_ij[deadHerbIndex] +
                                    siteBU.combustion_ij[liveHerbIndex]) / totalInitialLoading;
            }
          }
        }

        this->r_dead2ag_cn = 0;//Currently there is not mortality besides combustion.
      }
    }
    else//Woody PFTs:
    {
      if (IsShrub(cd->cmttype, pftIndex))
      {
        //Shrubs map to the live woody fuel:
        this->r_burn2ag_cn = siteBU.combustion_ij[liveWoodyIndex] / siteBU.w_o_ij_Initial[liveWoodyIndex];
        this->r_dead2ag_cn = 0;//Currently there is not mortality besides combustion.
      }
      else//Trees:
      {
        //For our purposes trees only burn by crown fire, indicated by a CFB greater than 0:
        if (CFB > 0.0)
        {
          //The extent of combustion is a linear function of the crown fraction burned:
          //We assume the maximum tree consumption is 0.35.  This is currently true for all tree
          //PFTs but this could change with an update to the parameter file.  We could get
          //fvcomb_sev5 for this PFT and use that but since the dead pool flux calculation depends
          //on the current value things could break.
          double burntFrac = 0.35 * CFB;

          //Calculate the total mortality (burnt + flux to dead pool) from the burnt fraction:
          //We fit an equation to the relationship between the extisting severity level parameters
          //to match the relationship in a continueous fashion:
          double totalMortFrac = 1.168170 * (1.0 - exp(-5.354249 * burntFrac));

          //We expect conifers to drive crown fire behavior but for mixed forests we assume that
          //decidous trees experiance comparable combustion and mortaity.  This may not be true but
          //we lack infomation to support a better scheme and this is constent with the previous
          //model.
          this->r_burn2ag_cn = burntFrac;
          this->r_dead2ag_cn = totalMortFrac - burntFrac;
        }
        else
        {
          //Note: If CFB = 0 the above calculations should yield zeros, but let's not take chances.
          this->r_burn2ag_cn = 0;
          this->r_dead2ag_cn = 0;
        }
      }
    }
  }

  this->r_live_cn = 1.0 - this->r_burn2ag_cn - this->r_dead2ag_cn;
}

/** Return the change in the litter stock after surface fire.
 * 
 * We map litter to the dead fuels.  The litter is spread over multiple size classes so we need to
 * put them back together to get the total litter change.  There may be other stocks mixed into the
 * dead fuel depending on the mode, specifically moss and / or herbaceous fuels.
 *
 * @returns The fraction of litter burnt in the surface fire.
 *
 * @note Used in the process wildfire model only.
 */
double WildFire::GetLitterBurntFraction() const
{
  BOOST_LOG_SEV(glg, debug) << "Entering WildFire::GetLitterBurntFraction()...";

  double litterBurntFraction = 0.0;//Return value.

  //If for any reason fire didn't actually burn skip the calculation:
  if (siteBU.burnoutTime > 0.0)
  {
    double litterInitial = 0.0;//The total dry litter mass before fire.
    double litterCombusted = 0.0;//The total dry litter mass after fire.

    int numDeadClasses = siteFM.NumDeadClasses();
    for (int i = 0; i < numDeadClasses; i++)
    {
      //Skip any cured dynamic fuels:
      //We don't need to check if this is a dynamic fuel model because DeadHerbaceousIndex() will
      //return -1 if not.
      if (i != siteFM.DeadHerbaceousIndex())
      {
        //If live moss has been combined with the fine dead fuel we need to separate it out:
        if ((i == 0) && md.fire_moss_as_dead_fuel)
        {
          double fineBurntFraction = 0.0;
          if (siteBU.w_o_ij_Initial[i] > 0.0)//Avoid divide by zero.
          {
            fineBurntFraction = siteBU.combustion_ij[i] / siteBU.w_o_ij_Initial[i];
          }

          //All fuel in a size class burns at the same rate:
          double fineLitterInitial = siteBU.w_o_ij_Initial[i] - GetNonVascularBiomass();
          litterInitial += fineLitterInitial;
          litterCombusted += fineLitterInitial * fineBurntFraction;//Just the burnt fine litter.
        }
        else
        {
          litterInitial += siteBU.w_o_ij_Initial[i];
          litterCombusted += siteBU.combustion_ij[i];
        }
      }
    }

    //We could check that litterInitial matches up with what we get from GetLitterRawC() here.
    
    if (litterInitial > 0.0)//Avoid divide by zero.
    {
      litterBurntFraction = litterCombusted / litterInitial;
    }
    
    if (!ValidProportion(litterBurntFraction))
    {
      BOOST_LOG_SEV(glg, fatal) << "Invalid litter burnt fraction calculated: " << litterBurntFraction;
    }
  }

  return litterBurntFraction;
}

/** Simulate crown fire.
 * 
 * @returns A vector of the heat per unit area (HPA) of the surface fire component, the crown fire
 *          component, and of the whole fire in that order (kJ/m^2).  If no crown fire occurs the
 *          vector will be all zeros.  CFB and totalFireIntensity are also updated.
 *
 * @note We return a vector of 0s because is seems a bit safer than returning an empty vector.
 *       However, the calling code will do better to check if CFB > 0 to determine if a crown fire
 *       occured.
 * @note Used in the process wildfire model only.
 */
std::vector <double> WildFire::SimulateCrownFire()
{
  std::vector <double> cHPAs = {0, 0, 0};

  BOOST_LOG_SEV(glg, debug) << "Entering SimulateCrownFire()...";

  //Check to see if the CMT can ecologically have a crown fire:
  //Negative prameter values for CBD and CBH signal that this CMT can't crown.
  if (firpar.cbd > 0.0 && firpar.cbh > 0.0)
  {
    //Gather crown fire model inputs:
    //The site's fuel model is saved in siteFM and reflects the pre-fire states.
    double WRF = GetWindReductionFactor();
    double U = GetMidflameWindSpeed();
    double slopeSteepness = SlopePctToSteepness(cd->cell_slope);//Or pass in?
    double FMC = CalculateFoliarMoistureContent();
    double CBD = GetCanopyBulkDensty();
    double CBH = GetCrownBaseHeight();
    FuelModel fuelModel10 = GetFuelModelFromCSV(md.fire_fuel_model_file, 10);

    BOOST_LOG_SEV(glg, debug) << "WRF = " << WRF;
    //Wind and slope are printed elsewhere, but we could print them here too.
    BOOST_LOG_SEV(glg, debug) << "FMC = " << FMC << " %";
    BOOST_LOG_SEV(glg, debug) << "CBD = " << CBD << " kg/m^3";
    BOOST_LOG_SEV(glg, debug) << "CBH = " << CBH << " m";

    //Calculate CFB:
    CFB = CrownFractionBurned(siteFM, U, WRF, slopeSteepness, CBD, CBH, FMC, fuelModel10, 'U');
    //The crown combustion fraction is stored for use in getAbgVegetationBurntFractionsProcess() to
    //derive combustion and non-combustion mortality.

    if (CFB > 0.0)//A crown fire initiates.
    {
      BOOST_LOG_SEV(glg, info) << "Crown fire occured. CFB = " << CFB;

      //Calculate crown fire heat per area:
      double CFL = GetCanopyFuelLoad();
      BOOST_LOG_SEV(glg, debug) << "CFL = " << CFL << " kg/m^2";
      cHPAs = CrownComponentHPA(siteFM, U, WRF, slopeSteepness, CBD, CBH, FMC, CFL, fuelModel10, 'U');

      //Calculate and store Byram's fire intensity:
      totalFireIntensity = CrownFireIntensity(siteFM, U, WRF, slopeSteepness, CBD, CBH, FMC, CFL,
                                              fuelModel10, 'U');//kW/m
      BOOST_LOG_SEV(glg, debug) << "Total fire intensity = " << totalFireIntensity << " kW/m";
    }
    else
    {
      BOOST_LOG_SEV(glg, info) << "Crown fire did not occur.";
    }
  }

  return cHPAs;
}

/** Simulate ground fire returning the burn depth.
 * 
 * This effectively replaces the functionality of WildFire::getBurnOrgSoilthick() in the original
 * wildfire implementation.
 *
 * @param fireHeatInput Total longwave heat input into the soil from the surface fire (kJ/m^2).
 * @note Burnup produces energy over time so it may be better to link the calculations?
 *
 * @returns The soil burn depth from ground fire (meters).
 *
 * @note Used in the process wildfire model only.
 */
double WildFire::SimulateGroundFire(const double fireHeatInput) const
{
  BOOST_LOG_SEV(glg, debug) << "Entering SimulateGroundFire()...";

  GFProfile gfProfile = GroundFireGetSoilProfile();//Get the soil profile information the model needs.
  double burnDepth = DownwardGroundFire(gfProfile, fireHeatInput, md.fire_gf_heat_loss_factor,
                                        md.fire_gf_surface_pd, md.fire_gf_smolder_pd) / 100.0;//Convert cm to meters.

  BOOST_LOG_SEV(glg, debug) << "Profile after DownwardGroundFire():";
  BOOST_LOG_SEV(glg, debug) << gfProfile;

  return burnDepth;
}

/** Get information about the soil profile and convert it to form that can be used by the ground
 * fire model.
 *
 * The main model's soil representation is made up of layers with varying thickness.  This function
 * creates a representation with layers of even thickness, converting, calculating, and
 * interpolating layer properties needed by the ground fire model.
 *
 * Possible parameters could include the desired layer thickness.
 *
 * We assume that the fire will not continue if it reached mineral soil, bedrock, permafrost, or the
 * water table.  However, permafrost is an annual property and can't be determined at a given month.
 * Likewise the water table is calculated from the profile.  We must rely on temperature and water
 * content to control this behavior.  [More?????]
 * Moss is treated as a surface fuel.
 * Snow needs to be considered.
 *
 * @returns A GFProfile object containing soil properties by layers of even thickness.
 *
 * @note Used in the process wildfire model only.
 */
GFProfile WildFire::GroundFireGetSoilProfile() const
{
  //The temperature of ignition of organic soil (C):
  //Based on literature review by JMR 12/10/2025.
  //This may be broken down by fibric and humic profile in the future if data becomes availalble.
  const double t_ig_organic = 250.0;//WILDFIRE_PARAMETER

  BOOST_LOG_SEV(glg, debug) << "Entering GroundFireGetSoilProfile()...";

  //Only consider the organic horizon(s):
  int numOrgLayers =  ground->organic.shlwnum + ground->organic.deepnum;

  //Create an object to hold the profile data:
  GFProfile gfProfile(numOrgLayers);

  //Walk through the organic soil layers and copy / convert properties to our new profile:
  /*There are two representations of the soil layers, the ground object and edall.  We walk both
  layer representaions simultaneously to get the values we need.  The ground object holds all the
  properties we need but the temperature and moisture values, which seem to be daily values from the
  end of the month, have been problematic.  edall does not have all the physical properties we need
  but is has both monthly and daily temperature and moisture values.*/
  Layer* thisLayer = ground->fstshlwl;
  int layerIndex = MAX_MOS_LAY;//Index of the top organic layer, right below the bottom moss layer.
  for (int i = 0; i < numOrgLayers; i++)
  {
    //Copy data from the source layer to the matching layer:
    gfProfile.thickness_cm[i] = thisLayer->dz * 100.0;//Layer thickness (m -> cm).
    gfProfile.layerDepth[i] = thisLayer->z * 100.0;//Depth at top of layer (m -> cm).
    gfProfile.tempC[i] = edall->m_sois.ts[layerIndex];//Layer temperature in Celcius.
    //Compare values from two soil representations:
    //gfProfile.tempC[i] = thisLayer->tem;//Layer temperature in Celcius.
    BOOST_LOG_SEV(glg, debug) << "Soil layer temp (C): ground Layer->tem = " << thisLayer->tem
                              << ", edall->m_sois.ts = " << edall->m_sois.ts[layerIndex];

    //The humic layers should probably have higher temperatures of ignition but more research is
    //needed.  The values should parameters or be calculated from the level of decay.  They will be
    //interpolated in the ground fire calculation.
    if (thisLayer->isFibric)
    {
      gfProfile.t_ig[i] = t_ig_organic;
      gfProfile.type[i] = "Fibric";
    }
    else if (thisLayer->isHumic)
    {
      gfProfile.t_ig[i] = t_ig_organic;
      gfProfile.type[i] = "Humic";
    }
    else//Same as checking !thisLayer->isOrganic.
    {
      BOOST_LOG_SEV(glg, fatal) << "Layer is not an expected organic type.";
    }

    gfProfile.bulkDensity[i] = thisLayer->bulkden / gPerKg;//Dry soil mass per volume (g/m^3 -> kg/m^3).

    //The soil organic / inorganic fractions are not explicit properties tracked by TEM.  Carbon is
    //used for accounting but this will be less than the total soil organic mass, that which is lost
    //on combustion.  We can estimate the SOM, the complement of the inorganic fraction, from SOC.
    //A SOC to SOM ratio of 0.5 is commonly used for soils without much support.  We collected 
    //literature data for histosols, specificlaly peats, and found support for this value.
    const double SOCtoSOM_Ratio = 0.5;//WILDFIRE_PARAMETER
    double totalSOC = 0.0;//g/m^2(/layer)
    if (i == 0)//We treat the rawc compartment of the top non-moss/fibric/shallow layer as litter:
    {
      totalSOC = thisLayer->soma + thisLayer->sompr + thisLayer->somcr;//g/m^2(/layer)

      //FW_NOTE: If the rawc is not considered here it should also be subracted from the layer mass below!
    }
    else
    {
      totalSOC = thisLayer->rawc + thisLayer->soma + thisLayer->sompr + thisLayer->somcr;//g/m^2(/layer)
    }
    //The soil organic carbon pools are in gC/m^2 per layer.  Use the layer depth to convert to density:
    //(gC/m^2 / m) / g/Kg = kgC/m^3
    double totalSOCdensity = (totalSOC / thisLayer->dz) / gPerKg;//kg/m^3
    //Convert from carbon to organic matter:
    double totalSOMdensity = totalSOCdensity / SOCtoSOM_Ratio;//kg/m^3
    //Or:
    //double totalSOMdensity = totalSOCdensity / thisLayer->cfrac;

    //Even if the soil is purely organic it can't exceed a fraction of 1.0 so cap it:
    double organicFraction = totalSOMdensity / gfProfile.bulkDensity[i];

    //Temporary code to report on the carbon in more detail:
    double totalSOCTemp = thisLayer->rawc + thisLayer->soma + thisLayer->sompr + thisLayer->somcr;//The total reguardless of layer.
    double totalSOCdensity2 = (totalSOCTemp / thisLayer->dz) / gPerKg;//kg/m^3
    double cFraction = totalSOCdensity2 / gfProfile.bulkDensity[i];
    BOOST_LOG_SEV(glg, debug) << "True total carbon fraction of layer = " << cFraction;
    if (organicFraction > 1.0)
    {
      BOOST_LOG_SEV(glg, debug) << "organicFraction = " << organicFraction;
      BOOST_LOG_SEV(glg, debug) << "Layer carbon = " << totalSOCTemp;
      BOOST_LOG_SEV(glg, debug) << "totalSOC = " << totalSOC;
    }

    organicFraction = fmin(organicFraction, 1.0);
    gfProfile.inorganicPct[i] = (1.0 - organicFraction) * 100;//Percent inorganic content on a dry basis (~ ash content).

    //Soil moisture content (%):
    //This is the water mass per volume / dry soil mass per volume * 100%.
    //Ice needs to be added to the ground fire model explicitly but for now we convert it to water.
    //liq is liquid water and ice is ice content in kg/m^2 per layer.  Use the layer depth to
    //convert to water mass per volume:
    //kg/m^2 / m = kg/m^3
    double waterContent = (edall->m_sois.liq[layerIndex] + edall->m_sois.ice[layerIndex]) /
                          thisLayer->dz;//kg/m^3
    //Compare values from two soil representations:
    double waterContentL = (thisLayer->liq + thisLayer->ice) / thisLayer->dz;//kg/m^3, Just for info.
    BOOST_LOG_SEV(glg, debug) << "Soil water content (kg/m^3): ground object = " << waterContentL
                              << ", edall = " << waterContent;

    //kg/m^3 / kg/m^3 * 100% = %
    gfProfile.moistureContentPct[i] = waterContent / gfProfile.bulkDensity[i] * 100.0;//%

    //Layer heat capacity (kJ/kg/K):
    //vhcsolid is volumetric heat capacity (J/m^3/K) for dry compacted soil. Converting to specific
    //heat capacity requires considering the pore fraction and the bulk density:
    //(J/m^3/K * solid fraction) / kg/m^3 * 1000J/KJ = kJ/kg/K
    gfProfile.c_s[i] = (thisLayer->vhcsolid * (1.0 - thisLayer->poro)) / gfProfile.bulkDensity[i] / 1000.0;//kJ/kg/K

    thisLayer = thisLayer->nextl;
    layerIndex += 1;
  }

  BOOST_LOG_SEV(glg, debug) << "Initial translated profile:";
  //BOOST_LOG_SEV(glg, debug) << gfProfile;
  //Print with full numeric resolution to enable exact reimport.  We match the default description
  //line and start the profile on a new line for automated log searches.
  gfProfile.SetPrintMode(-1);
  BOOST_LOG_SEV(glg, debug) << "Ground fire soil profile:" << '\n' << gfProfile;
  gfProfile.SetPrintMode(0);//Restore.

  //If there is a moss layer the top organic layer will not start at depth zero.  Adjust for this:
  gfProfile.Resurface();

  BOOST_LOG_SEV(glg, debug) << "Profile after Resurface():";
  BOOST_LOG_SEV(glg, debug) << gfProfile;

  //Check the profile before interpolating:
  if (!gfProfile.Validate())
  {
    BOOST_LOG_SEV(glg, fatal) << "Translated profile it not valid.";
  }

  //Convert to layers of equal thickness and interpolate the values in the original profile:
  gfProfile.Interpolate(md.fire_gf_layer_thickness);
  //It is not necissary to print the profile here since it will be printed after the ground fire
  //calculation.  The interpolated profile infomation witl not have changed.

  return gfProfile;
}
