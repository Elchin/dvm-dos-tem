/***********************************************************************************************//**
 * @file GroundFire.cpp
 * \author Joshua M. Rady
 * Woodwell Climate Research Center
 * \date 2025
 *
 * @brief This file defines functions to simulate smoldering ground fire using a moderate complexity
 * model where combustion is predicted layer by layer.
 *
 **************************************************************************************************/

#include <algorithm>//For fill().
#include <cmath>
#include <numeric>//For accumulate().
#include <string>

#include "../include/GroundFire.h"//"GroundFire.h"//Note: Difference from source repo.
#include "FireweedMessaging.h"
#include "FireweedUtils.h"

//Constants:----------------------------------------------------------------------------------------
//Properties of water:
const double DelatH_vap = 2257.0;//ΔH_vap = (Latent) heat of vaporization of water (J/g or kJ/kg).
const double c_w = 4.184;//The specific heat capacity of water (kJ/kg/K at 20C).
const double DelatH_fus = 333.55;//ΔH°_fus = (Latent) heat of fusion of water (J/g or kJ/kg).
const double c_i = 2.050;//The specific heat capacity of ice (kJ/kg/K at -10C).

//Functions:----------------------------------------------------------------------------------------

/** Calculate and return the heat of combustion of an organic soil based on its organic content.
 *
 * This is a simple calculation that may not be robust for all organic soils.  It assumes that the
 * inorganic content of the soil has no effect on the energy released during combustion and that the
 * pure organic material has a fixed heat of combustion.  The former may be mostly true but the
 * latter is certainly not.  The energy stored in organic material differs by source and level of
 * decomposition.  This function is biased towards rich organic soils, like peats, where we hope
 * these assmumptions will be less problematic.  The organicHOC parameter can also be changed to
 * account for differences in organic quality if information is available.
 * 
 * Additionally the heat of combustion is generally considered for the case of combustion with excess
 * oxygen.  For flaming combustion this may be appropriate but for smoldering combustion can occur by
 * more than one pathway, with different energies.
 * 
 * The smoldering specific formulation includes the component oxidations:
 * ΔH_c = ΔH peat beta char oxidation + ΔH alpha char oxidation + ΔH beta char oxidation
 * However, we only have those from a few papers so we can only make a first order estimate.  We may
 * be able to improve this function over time with more research.
 * 
 * @param[in] inorganicPct The percent organic content of the soil on the dry basis.
 * @param[in] organicHOC The (theoretical) heat of combustion for a purely organic soil (MJ/kg).
 * We can estimate that from a soil's heat of combustion and it's inorganic / ash content as:
 * soil HOC / (100% - soil % Ash) * 100% = organicHOC
 * Our current mean value calculated from organic soils in the literature (duff & peat) is 21.6
 * MJ/kg.
 *
 * @returns The soil heat of combustion (MJ/kg, ~net fuel heat content). (ΔH°c)
 */
double SoilHeatOfCombustion(const double inorganicPct, const double organicHOC)
{
	if (!InRange(inorganicPct, 0.0, 100.0))
	{
		Stop("Invalid inorganic percentage.");
	}

	double soilHOC = organicHOC * (100.0 - inorganicPct) / 100.0;//MJ/kg
	return soilHOC;
}

/** Return the effective the heat of combustion for smoldering organic soil.
 *
 * The amount of heat released by smoldering combustion is lower than the absolute heat content of
 * the fuel for two reasons.  First, in smoldering the pyrosylates do not burn via flaming and are
 * lost tp the atmosphere (though some may recondense as part of the char?).  Second, combustion may
 * not be complete.  We currently have estimates for the smoldering reduction but don't currently
 * have a general estimate of how complete combustion is.
 *
 * @param[in] inorganicPct The percent organic content of the soil on the dry basis.
 *
 * @returns The effective organic soil heat of combustion (MJ/kg).
 */
double EffectiveSmolderingHOC(const double inorganicPct)//EffectiveSmolderingHeatOfCombustion()
{
	//Frandsen 1991 gives the efficiency of smoldering vs. complete oxydataion at 73%.
	//We are seeking more estimates.
	const double smolderingEfficiency = 0.73;
	//const double combustionEfficiency = X;//May add in the future.

	return SoilHeatOfCombustion(inorganicPct) * smolderingEfficiency;
}

/** Calculate the heat sink for a kg of dry organic soil and its associated soil water (as % MC).
 * This is the sum of the energy required for (melting,) drying, and to further heat the dry soil
 * to ignition.
 *
 * @param[in] moistureContentPct = Percent moisture content (water mass / dry soil mass * 100%).
 * @param[in] t_int Intial soil temperature (C).
 * @param[in] t_ig Temperature of (auto/self-)ignition or combustion (C).
 * @param[in] c_s Specific heat capacity of the soil (kJ/kg/K). [at 20C?]
 *        This is a function of the soil component heat capacities and the temperature, though for
 *        now  we ignore the temperature dependence as is frequently done.
 *
 * @returns The heat sink (kJ/kg as dry soil mass).
 */
double SoilHeatSink(const double moistureContentPct, const double t_int, const double t_ig, const double c_s)
{
	//Here we use degrees C but we could use K.

	//We use the boiling temperature of water at sea level but it would be better to calculate it
	//from the altitude or air pressure:
	const double tBoil_H20 = 100.0;//The boiling point of water.

	double t_H2O = t_int;//The temperature of the soil water.
	//For 1 kg dry soil the soil water (kg) = 1 kg * percent moisture content / 100%:
	double mass_H20 = moistureContentPct / 100.0;

	//If soil ice is present calculate the heat required to melt it (kJ/kg ice):
	double DeltaH_melt = 0.0;
	if (t_int < 0.0)
	{
		//This is the sum of the heat to warm the ice to 0 and the heat to melt it.
		DeltaH_melt = (0.0 - t_int) * c_i + DelatH_fus;
		t_H2O = 0.0;//Warm the water for the next step.
	}
	
	//Calculate the heat of drying (kJ/kg water):

	//Soil water is in free and bound forms.  The bound water fraction requires additional energy to
	//evaporate.  I'm still trying to figure out how to estimate what the bound fraction is.  There may
	//even be more than one bound fraction as hygroscopic and tightly bound water may be different?
	//As a start we could exclude free water above ~100 - 105% moisture content?
	//These terms are being set to zero until I figure this out.
	double frac_bound = 0.0;//The fraction of water in a bound form.
	double DeltaH_des = 0.0;//ΔH_des = Heat of desorption of water (kJ/kg)

	double DeltaH_dry = (tBoil_H20 - t_H2O) * c_w + DelatH_vap + (frac_bound * DeltaH_des);

	//Heat of ignition of dry soil (kJ/kg soil):
	double DeltaH_ig = (t_ig - t_int) * c_s;

	//Total heat required, (kJ/kg moist soil):
	double h_sink = (DeltaH_melt * mass_H20) + (DeltaH_dry * mass_H20) + DeltaH_ig;//1 kg dry soil is implied.
	return h_sink;
}

/** Calculate the distribution of heat over a number of layers, decreasing linearly with depth.
 *
 * The distribution is treated as triangular with the fraction of the distribution in a layer given
 * by the area of a slice of that triangle.  We assume that we are distributing over full layers.
 * The function returns fractional weights adding to one.  The weight order is top to bottom.  The
 * actual heat per layer can be obtained by multiplication.
 *
 * @param[in] numLayers The number of layer for which to compute the distribution.
 *
 * @returns A vector of fractional weights adding to one for each layer from top to bottom.
 */
std::vector <double> HeatDistributionLinear(const int numLayers)
{
	if (numLayers < 1)
	{
		Stop("The number of layer must be at least 1.");
	}

	std::vector <double> layerWeights(numLayers);
	double layerThickFrac = 1.0 / numLayers;//The relative thickness of layers as a fraction of the whole depth.

	for (int i = 1; i <= numLayers; i++)
	{
		//The distribution is a right triangle so the width can be defined by a line of form m*x + b,
		//or m*y + b since it is facing down.  We don't know anything about were the layers start, end,
		//or how deep they are.  We treat everything as relative. For convenience we can make the
		//dimensions of the triangle 2 x 1 so it has an area of 1.   It follows that m = 2.  Since we
		//are starting from the top we have to set the intercept b = 2.  Each layer slice has a depth
		//of 1 / numLayers.
		double m = -2.0;
		double b = 2.0;
		
		//The depths of the top and bottom of the slice:
		double d_top = (i - 1) * layerThickFrac;
		double d_bottom = i * layerThickFrac;
		
		//The area for a layer from depth d1 to d2 = (d2 - d1) * (w1 + w2)/2, where w is the width of
		//the triangular distribution at depth d.  d2 - d1 always = layerThickFrac:
		double w_top = m * d_top + b;
		double w_bottom = m * d_bottom + b;
		
		layerWeights[i - 1] = layerThickFrac * (w_top + w_bottom) / 2.0;
	}

	if (!FloatCompare(std::accumulate(layerWeights.begin(), layerWeights.end(), 0.0), 1.0))
	{
		Stop("Weights do not sum to 1.");
	}

	return layerWeights;
}

/** Calculate the depth of burn for a generic soil profile.
 * 
 * Ground fire is simulated as smoldering combustion that moves progressively in the downward
 * direction.  The soil column is represented as layers of equal thickness.  The soil column is
 * treated as homogeneous in the vertical direction and vertical spread is not considered.
 *
 * Ground fire initiates when energy input from the aboveground fire into the soil drys and heats
 * the top layer to the temperature of ignition.  Currently the top layer is not treated in a
 * special way.  However, given the importance of ignition we may change this in the future to
 * simulate the surface fire input more explicitly in time and with a full heat transfer
 * representation.
 *
 * The fire extinguishes when the energy to dry and heat the next layer to combustion is less than
 * that transfered from combustion above.  The model uses simplified representations of the process
 * of heat transfer and loss.  The calculation occurs layer by layer without the timing of
 * progression simulated directly, i.e. the model used layer-steps rather than time-steps.
 *
 * @param[in,out] soilCol A GFProfile object representing the soil profile for the simulation.
 *                        A profle with layers of equal thickness is expected.
 * @param[in] fireHeatInput Total heat input into the soil from the surface fire (kJ/m^2).
 * @param[in] heatLossFactor The fraction of heat lost from the soil (0-1) .  This intended to
 *                           represent radiant and convective losses from the soil surface in the
 *                           time is take for a layer to combust. The downward efficiency factor
 *                           from Benscoter et al. 2011 = downwardEF. As a heat loss = 0.83.
 * @param[in] surfacePD The maximum heat penetration depth for the soil surface heat transfer
 *                      approximation (cm).
 * @param[in] smolderPD The maximum heat penetration depth for the smoldering heat transfer
 *                      approximation (cm).
 * The following parameters may be temporary:
 * @param[in] surfaceTM The surface heat transfer mode.  0 to distribute the heat evenly, 1 to
 *                      distribute the heat in a linearly decreasing manner with depth.
 * @param[in] smolderTM The smoldering heat transfer mode.  0 to transfer all heat to the next
 *                      layer, 1 to distribute the heat across a depth of smolderPD in a linearly
 *                      decreasing manner with depth.
 *
 * @returns The burn depth (cm).  Additionally the GFProfile passed is updated on return.
 */
double DownwardGroundFire(GFProfile& soilCol, const double fireHeatInput,
                          const double heatLossFactor, const double surfacePD,
                          const double smolderPD, const int surfaceTM, const int smolderTM)
{
	//If there is no organic soil there is nothing to burn:
	//Note: This check may be integrated with validity checking in the future.
	if (soilCol.NumLayers() == 0)
	{
		return 0.0;
	}

	//Validity checking:
	if (!soilCol.Validate())//This is may be overkill since Validate() will be called with interpolation.  However, in the testing setting is useful.
	{
		Stop("DownwardGroundFire(): The soil column is not valid.");
	}

	if (fireHeatInput < 0.0)
	{
		Stop("Invalid value for fireHeatInput.");
	}
	if (!ValidProportion(heatLossFactor))
	{
		Stop("Invalid value for heatLossFactor.");
	}

	//Check that the layers have equal thickness:
	//If the profile was interpolated the layer thicknesses will have been checked but we can't know
	//that.
	if (!soilCol.EqualThickness())
	{
		Stop("DownwardGroundFire(): Soil layers do not all have uniform thickness.");
	}
	double layerThickness_cm = soilCol.thickness_cm[0];

	if (surfacePD < layerThickness_cm)
	{
		Stop("surfacePD is less than a layer thickness.");
	}
	if (smolderPD < layerThickness_cm)
	{
		Stop("smolderPD is less than a layer thickness.");
	}

	//The single layer smoldering heat transfer scheme is still available but is primarily intended
	//for testing so make it clear when it is being used:
	if (smolderTM == 0)
	{
		Warning("Using single layer (pure domino) heat transfer scheme.");
	}

	//Calculate the number of layers to distribute the heat from the surface fire to:
	//The depth the heat from the surface fire is transfered to is provided to the model.  In
	//reality the depth and profile of heat will be a fuction of the heat input pattern over time
	//and soil properties. The fact that this parameter is only used once means that it could be
	//calculated on a simulation by simulation basis, rather than being fixed.
	//The number of layers needs to be a whole number.  We round down but it might be better to
	//require that the layer thickness be a divisor of the surface penetration depth.
	int numSurfLayers = surfacePD / layerThickness_cm;

	//Calcualte the heat distibution as a set of weights for each layer:
	std::vector <double> surfWts(numSurfLayers);
	if (surfaceTM == 0)//Distribute the heat evenly over the layers:
	{
		std::fill(surfWts.begin(), surfWts.end(), 1.0 / numSurfLayers);
	}
	else//Distribute the heat linearly over the layers:
	{
		surfWts = HeatDistributionLinear(numSurfLayers);
	}

	//Distribute the heat from the aboveground fire into the soil surface:
	for (int i = 0; i < numSurfLayers; i++)
	{
		double heatToLayerI = surfWts[i] * fireHeatInput;

		if (heatToLayerI < 0.0)//This check should probably be replaced by checks on the components.
		{
			Stop("Negative heat flux for surface layer" + std::to_string(i));
		}

		if (i < soilCol.NumLayers())//Don't walk off the bottom of the soil column.
		{
			soilCol.heatSource[i] = heatToLayerI;
		}
	}

	//For each layer see if the heat from above (the fire or previous layer) is enough to ignite it:
	for (int i = 0; i < soilCol.NumLayers(); i++)
	{
		//The heat source for this layer was previously calculated and stored:
		double heatSourceToLayer = soilCol.heatSource[i];

		//The (total) heat sink for a layer is determined exclusively by its physical properties and
		//the mass of soil per layer:
		soilCol.heatSink[i] = SoilHeatSink(soilCol.moistureContentPct[i], soilCol.tempC[i],
		                                   soilCol.t_ig[i], soilCol.c_s[i]) *
		                                   soilCol.DrySoilMassKg(i);

		if (heatSourceToLayer >= soilCol.heatSink[i])
		{
			soilCol.burnt[i] = true;

			if (i == soilCol.NumLayers() - 1)//If we have reached the bottom there is nowhere to transfer heat to:
			{
				break;
			}

			//Energy in excess of that needed to dry and heat this layer has to go somewhere.  Some
			//will be lost and some will heat lower layers.
			double layerExcessHeat = heatSourceToLayer - soilCol.heatSink[i];

			//Calculate the heat released from combustion for this layer:
			double layerHeatOfCombustion = EffectiveSmolderingHOC(soilCol.inorganicPct[i]) *
				soilCol.DrySoilMassKg(i) * 1000.0;//MJ/kg -> kJ/kg

			//The total heat available from this layer:
			double totalLayerHeat = layerHeatOfCombustion + layerExcessHeat;

			//Compute the heat transferred to lower layer(s):
			double heatSourceDown = totalLayerHeat * (1.0 - heatLossFactor);//Or propagating heat?

			//Note: Storing these values (heatSourceToLayer, layerExcessHeat, layerHeatOfCombustion,
			//totalLayerHeat, heatSourceDown) could be useful for interpreting these calulations.
			//Some of these were stored in an earier version of the R code.

			//Compute how the heat is distributed:
			if (smolderTM == 0)//Single layer (pure domino) heat transfer scheme:
			{
				//In the simplest representation all the heat is transferred to the next layer:
				soilCol.heatSource[i + 1] = heatSourceDown;
			}
			else//Simple linear heat transfer scheme:
			{

				//Distribute the heat to layers below:
				int distLayers = round(smolderPD / layerThickness_cm);//The number of layers to distribute heat to.
				//The top and bottom layers to distribute heat to:
				int topLayer = i + 1;
				int bottomLayer = i + distLayers;
				std::vector <double> wts = HeatDistributionLinear(distLayers);

				for (int j = topLayer; j <= bottomLayer; j++)
				{
					double heatToLayerJ = wts[j - topLayer] * heatSourceDown;
		
					if (heatToLayerJ < 0.0)
					{
						Stop("Negative heat flux for layer " + std::to_string(j));
					}

					if (j < soilCol.NumLayers())//Don't walk off the bottom of the soil column.
					{
						//Accumulate the heat:
						soilCol.heatSource[j] += heatToLayerJ;
					}
				}
			}
		}
		else
		{
			break;
		}
	}

	return soilCol.GetBurnDepth();
}
