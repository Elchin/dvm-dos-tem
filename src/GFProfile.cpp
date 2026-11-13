/***********************************************************************************************//**
 * @file GFProfile.cpp
 * \author Joshua M. Rady
 * Woodwell Climate Research Center
 * \date 2025
 *
 * @brief This file defines an object used to represent a soil column during simulations of
 * smoldering ground fire.
 *
 **************************************************************************************************/

#include <algorithm>//For max().
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>//For accumulate().
#include <string>

#include "../include/GFProfile.h"//Note: Difference from source repo.
#include "FireweedMessaging.h"
#include "FireweedStringUtils.h"
#include "FireweedUtils.h"

//Public Functions:---------------------------------------------------------------------------------

/** Default constructor
 *
 * This constructor creates an object with all empty vectors.  If the number of layer is know at
 * initialization the alternate constructor may be preferable.
 */
GFProfile::GFProfile()
{
	numLayers = 0;
	printMode = 0;
}

/** Constructor
 *
 * This constructor makes an 'blank' profile of the specified size.  The soil properties are set to
 * invalid values that should be replaced after initialization.
 * 
 * @param[in] The number of profile layers to create.
 */
GFProfile::GFProfile(const int numLayers)
{
	this->numLayers = numLayers;
	printMode = 0;

	//Size the elements and set to invalid values:
	thickness_cm.resize(numLayers, ProfileUndef);
	layerDepth.resize(numLayers, ProfileUndef);
	tempC.resize(numLayers, ProfileUndef);
	t_ig.resize(numLayers, ProfileUndef);
	bulkDensity.resize(numLayers, ProfileUndef);
	inorganicPct.resize(numLayers, ProfileUndef);
	moistureContentPct.resize(numLayers, ProfileUndef);
	c_s.resize(numLayers, ProfileUndef);

	//Outputs:
	burnt.resize(numLayers, false);
	heatSink.resize(numLayers, 0.0);
	heatSource.resize(numLayers, 0.0);

	//Labels:
	type.resize(numLayers);//Fill with empty strings.
}

/** Get the number of layers.
 * 
 * @returns The number of layers in the soil column.
 */
int GFProfile::NumLayers() const
{
	return numLayers;
}

/** Return the dry soil mass for a layer.
 * 
 * @param[in] The index of the layer.
 * 
 * @returns The dry soil mass (kg/m^2/layer).
 */
double GFProfile::DrySoilMassKg(const int layerIndex) const
{
	if (layerIndex < 0 || layerIndex >= numLayers)
	{
		Stop("DrySoilMassKg(): Invalid layer index: " + std::to_string(layerIndex));
	}

	return bulkDensity[layerIndex] * (thickness_cm[layerIndex] / 100);//kg/m^3 * (cm / 100 cm/m) / layer = kg/m^2/layer.
}

/** Interpolate the soil profile.
 *
 * Check the profile and expand it to have layers of equal thickness and interpolate the current
 * values for the new layers, if necessary.
 *
 * @param[in] newLayerThickness The desired layer thickness (cm) to convert to.
 *
 * @returns Nothing.
 *
 * @note Our interpolation process begins by nudging the properties of the original layers to the
 * nearest new layer, on the basis of the layer centers.  The profile is then interpolated from
 * there.  This approach preserves the original values exactly in the profile, but their depths
 * may shift a bit.  The alternative would be to interpolate everywhere.  This prioritise the depths
 * over the exact values, since all layers would be interpolated in this case.  If the new layers
 * are thin the difference between the two approaches should be small.
 *
 * @note An issue for our interpolation is that each layer has a thickness.  If we associate our
 * know values with the middle of the layer this causes some issues.  If the top or bottom layer are
 * thick we don't have any information to interpolate on the outsides.  We fill new layers at the
 * top and bottom with the properties of original top and bottom layers.
 *
 * @note Incoming layers that are thinner than the new layer thickness have to be handled.  The
 * approach used here is to drop layers that are less than half the new layer thickness.  Loss of
 * very thin layers is not likely to have much of an effect on results.  If layer values were
 * interpolated everywhere we could use the values of these thin layers.
 * It is possible for all layers be dropped if they are all thin.  This results in an invalid
 * profile.  Calling code needs to be aware of this.  It might be advantageous to be able to signal
 * that the layers were removed to distinguish from a badly formed profile.
 *       Further examination is needed to determine if thin layers can end up nudged to the same
 * position, which would be a problem
 */
void GFProfile::Interpolate(const double newLayerThickness)
{
	if (newLayerThickness <= 0)
	{
		Stop("Interpolate(): Invalid value for newLayerThickness: " + std::to_string(newLayerThickness));
	}
	//Check that newLayerThickness is some small round value?

	//If the layer thickness or depths are off the interpolation may succeed but the output will be compromised:
	if (!Validate())
	{
		Stop("Interpolate(): The soil column is not valid.");
	}

	//At the end we will have a set of new layers.  The bottom of the deepest layer should be pretty
	//close to that of the current profile but will be a multiple of the new newLayerThickness:
	int initBottomIndex = numLayers - 1;//The position (index) of the bottom layer before interpolation.
	//This should be used carefully since the numbers of layers will change during this function.

	//Check if interpolation is necessary:
	bool needed = false;
	for (int i = 0; i < numLayers; i++)
	{
		if (thickness_cm[i] != newLayerThickness)//We don't give any wiggle room.
		{
			needed = true;
			break;
		}
	}

	if (needed)
	{
		//The first layer should be at zero, otherwise there are problems.
		if (layerDepth[0] != 0.0)
		{
			Stop("The top layer is not at the surface!");
			//we could dump the layer here.
		}

		//Check for thin layers and drop them.  Start from the bottome so the a index of the next
		//layer doesn't change as we go:
		for (int j = initBottomIndex; j >= 0; j--)
		{
			//Keep layers greater than half the new layer thickness:
			if (thickness_cm[j] <= (newLayerThickness / 2))
			{
				Warning("Layer index " + std::to_string(j) + " is thin and is being dropped.");//Change to info?
				DropLayer(j);

				//If all the layers have dropped there is no interpolaition needed:
				//We could record the status here.
				if (numLayers == 0)
				{
					return;
				}

				if (j == 0)//If the top layer is dropped adjust the surface of the profile:
				{
					Resurface();
				}
			}
		}

		//Layers may have been dropped so update the bottom layer index before interpolation:
		initBottomIndex = numLayers - 1;//Maybe change the name since the value may have changed?
		//Store the original bottom of the soil column before making futher changes to it.
		double columnBottom = layerDepth[initBottomIndex] + thickness_cm[initBottomIndex];

		//Nudge the existing layers to the nearest equidistant layer postion:
		//The most we will move a value is by newLayerThickness / 2, which will probably not matter
		//if newLayerThickness is small.  This process may mean that the top and bottom layers are
		//no longer at the top and bottom (see notes).
		//Thin layers can complicate these calcaultions somewhat.
		double addedThickness = 0.0;
		for (int j = 0; j < numLayers; j++)
		{
			//Adjust layer depths due to expansion of thin layers to the minimum layer thickness:
			if (thickness_cm[j] < newLayerThickness)
			{
				addedThickness += newLayerThickness - thickness_cm[j];
				columnBottom += newLayerThickness - thickness_cm[j];
			}
			if (j + 1 < numLayers)
			{
				layerDepth[j + 1] += addedThickness;
			}

			double layerCenter = layerDepth[j] + (thickness_cm[j] / 2);//The center depth of the original layer.

			//Nudge the original layer properties to the new layer that the original center falls in:
			int newLayerIndex = layerCenter / newLayerThickness;//Round down / truncate the division to get the layer.

			//Without the depth adjustement above it is nessisary to check and adjust the postion in
			//some cases.  Leaving until further testing is complete.
			if (newLayerIndex < j)
			{
				//newLayerIndex = j;
				Warning("This layer is being nudged to a shallower postion and expected.");
			}

			layerDepth[j] = newLayerIndex * newLayerThickness;
			thickness_cm[j] = newLayerThickness;
		}

		//Depth calculations (perform after layers are dropped and depths are tweaked):
		//Calculate the adjusted bottom of the column making sure each current layer translates to at least one new layer:
		double newColumnBottom = std::max(std::round(columnBottom / newLayerThickness) * newLayerThickness,
		                                  numLayers * newLayerThickness);
		double newBottomLayerDepth = newColumnBottom - newLayerThickness;

		//Start the interpolation from the bottom of the column.  By doing this any layers added
		//will be below the next layer we have to examine.  That means the layer indexes won't
		//change as we are going, which makes things less confusing.

		//Add layers below the bottom layer if needed:
		//If the bottom of the lowest layer after nudging is where we expect the new bottom to be it
		//is not necessary.
		if (!FloatCompare(layerDepth[numLayers - 1], newBottomLayerDepth))
		{
			int numNewLayers = std::round((newBottomLayerDepth - layerDepth[numLayers - 1]) /
			                               newLayerThickness);
			AddLayersToBottom(numNewLayers);
		}

		//Starting from the original bottom layer see if layers need to be interpolated between them:
		for (int j = initBottomIndex; j > 0; j--)
		{
			//See if adjacent layers are one new layer thickness apart:
			double deltaZ = layerDepth[j] - layerDepth[j - 1];
			//if (deltaZ != newLayerThickness)
			if (!FloatCompare(deltaZ, newLayerThickness))
			{
				//If not add layers to fill in between them:
				//We use round() to avoid inappropriate rounding down when we cast to int. 
				int numNewLayers = std::round((deltaZ / newLayerThickness)) - 1.0;
				
				//The original layers should have been nudged to regular positions and the distance
				//between them should be some exact multiple of the new layer thickness.  Check to
				//be sure:
				if (!FloatCompare((numNewLayers + 1) * newLayerThickness, deltaZ))
				{
					Stop("Distance between layers is not a multiple of the new layer thickness.");
				}
			
				//Interpolate values for the new layers being inserted:
				InterpolateBetween(j - 1, j, numNewLayers);
			}
		}

		//Add layers above the top layer if needed:
		if (layerDepth[0] != 0.0)
		{
			int numNewLayers = std::round(layerDepth[0] / newLayerThickness);
			AddLayersToTop(numNewLayers);
		}

		//Perform quality control:
		if (!Validate(true))
		{
			Print(std::cout);
			//This could be made fatal.
		}
	}
}

/** Interpolate new layers between two soil layers.
 *
 * The function adds new layers between the layers specified and interpolates the soil physical
 * properties evenly across them.  If there is space between the parent layers then the new layers
 * will be given equal thicknesses to fill the gap. Otherwise they will be left blank and a warning
 * will be given. The interpolation does not change the parent layer thicknesses.
 *
 * @param[in] topLayer The top layer index to interpolate from.
 * @param[in] bottomLayer The bottom layer index to interpolate to.
 * @param[in] numNewLayers The number of layers to interpolate between them.
 *
 * @note Specifying the top and bottom layers for interpolation makes this function call clear to
 * read but specifying the bottom is not really necessary since it only makes sense to interpolate
 * between adjacent layers.  We could allow writing over intervening layers but the utility of that
 * is not clear.
 *
 * @note The proper interpolation for some properties, like t_ig, might not be linear and further
 *       thought is needed.  An alternative could be to interpolate some properties and compute
 *       others from those values.
 *
 * @returns Nothing.
 */
void GFProfile::InterpolateBetween(const int topLayer, int bottomLayer, const int numNewLayers)
{
	double newLayerThickness = ProfileUndef;

	if (bottomLayer != (topLayer + 1))
	{
		Stop("InterpolateBetween() can only interpolate between adjacent layers.");
	}

	//Insert new layers for the interpolation.
	InsertBlankLayers(bottomLayer, numNewLayers);
	int firstNewLayer = bottomLayer;//The position of the top new layer.
	bottomLayer = bottomLayer + numNewLayers;//The bottom layer moves down.

	//Check the parent layer spacing:
	double zBetween = layerDepth[bottomLayer] - (layerDepth[topLayer] + thickness_cm[topLayer]);
	if (zBetween == 0.0)//The layers are touching.
	{
		Warning("The parent layers are touching, leaving interpolated layer thicknesses and depths blank.");
	}
	else if (zBetween > 0.0)//There is space between layers.
	{
		newLayerThickness = zBetween / numNewLayers;
	}
	else//Negative values indicate an error in the profile.
	{
		Stop("Layer spacing doesn't make sense.");
	}

	//Interpolate the properties:
	//Calculate the interpolation steps as fractional distances between the two values:
	double stepSize = 1.0 / (numNewLayers + 1);
	double step = stepSize;

	for (int i = 0; i < numNewLayers; i++)
	{
		thickness_cm[firstNewLayer + i] = newLayerThickness;
		layerDepth[firstNewLayer + i] = layerDepth[firstNewLayer + i - 1] + thickness_cm[firstNewLayer + i];
		
		tempC[firstNewLayer + i] = std::lerp(tempC[topLayer], tempC[bottomLayer], step);
		t_ig[firstNewLayer + i] = std::lerp(t_ig[topLayer], t_ig[bottomLayer], step);
		bulkDensity[firstNewLayer + i] = std::lerp(bulkDensity[topLayer],
		                                           bulkDensity[bottomLayer], step);
		inorganicPct[firstNewLayer + i] = std::lerp(inorganicPct[topLayer],
		                                            inorganicPct[bottomLayer], step);
		moistureContentPct[firstNewLayer + i] = std::lerp(moistureContentPct[topLayer],
		                                                  moistureContentPct[bottomLayer], step);
		c_s[firstNewLayer + i] = std::lerp(c_s[topLayer], c_s[bottomLayer], step);

		step += stepSize;
	}
	//Leave burnt, heatSink, and heatSource at their defaults.
	//Leave type blank.
}

/** Add layers to the top of the profile, cloning the properties of the current top layer.
 *
 * This function is used when the surface layer is being sliced into thinner layers and the code
 * makes assumptions consistant with this.  It would also probably be safe to assume that after this
 * the new top layer will be at the surface, but we don't currently check for that here.
 *
 * @param[in] numNewLayers The number of layers to add.
 *
 * @returns Nothing.
 */
void GFProfile::AddLayersToTop(const int numNewLayers)
{
	if (numNewLayers < 1)
	{
		Stop("AddLayersToTop(): Invalid value for numNewLayers: " + std::to_string(numNewLayers));
	}

	if (numLayers == 0)
	{
		Stop("AddLayersToTop(): Can't add layers to empty profile.");
	}

	double layerThickness = thickness_cm[0];
	thickness_cm.insert(thickness_cm.begin(), numNewLayers, layerThickness);
	
	for (int i = numNewLayers - 1; i >= 0; i--)
	{
		layerDepth.insert(layerDepth.begin(), layerThickness * i);
	}

	//Clone physical properties:
	tempC.insert(tempC.begin(), numNewLayers, tempC[0]);
	t_ig.insert(t_ig.begin(), numNewLayers, t_ig[0]);
	bulkDensity.insert(bulkDensity.begin(), numNewLayers, bulkDensity[0]);
	inorganicPct.insert(inorganicPct.begin(), numNewLayers, inorganicPct[0]);
	moistureContentPct.insert(moistureContentPct.begin(), numNewLayers, moistureContentPct[0]);
	c_s.insert(c_s.begin(), numNewLayers, c_s[0]);

	//The profile should only be modified before a simulation so the outputs are at initial values:
	burnt.insert(burnt.begin(), numNewLayers, false);
	heatSink.insert(heatSink.begin(), numNewLayers, 0.0);
	heatSource.insert(heatSource.begin(), numNewLayers, 0.0);

	type.insert(type.begin(), numNewLayers, type[0]);//Copy the label.

	numLayers = numLayers + numNewLayers;
}

/** Add layers to the bottom of the profile, cloning the properties of the current bottom layer.
 *
 * This function is used when the bottom layer is being sliced into thinner layers and the code
 * makes assumptions consistant with this.
 *
 * @param[in] numNewLayers The number of layers to add.
 *
 * @returns Nothing.
*/
void GFProfile::AddLayersToBottom(const int numNewLayers)
{
	if (numNewLayers < 1)
	{
		Stop("AddLayersToBottom(): Invalid value for numNewLayers: " + std::to_string(numNewLayers));
	}

	if (numLayers == 0)
	{
		Stop("AddLayersToBottom(): Can't add layers to empty profile.");
	}

	thickness_cm.insert(thickness_cm.end(), numNewLayers, thickness_cm[numLayers - 1]);
	
	for (int i = 0; i < numNewLayers; i++)
	{
		layerDepth.push_back(layerDepth[numLayers - 1 + i] + thickness_cm[numLayers - 1]);
	}

	//Clone physical properties:
	tempC.insert(tempC.end(), numNewLayers, tempC[numLayers - 1]);
	t_ig.insert(t_ig.end(), numNewLayers, t_ig[numLayers - 1]);
	bulkDensity.insert(bulkDensity.end(), numNewLayers, bulkDensity[numLayers - 1]);
	inorganicPct.insert(inorganicPct.end(), numNewLayers, inorganicPct[numLayers - 1]);
	moistureContentPct.insert(moistureContentPct.end(), numNewLayers, moistureContentPct[numLayers - 1]);
	c_s.insert(c_s.end(), numNewLayers, c_s[numLayers - 1]);

	//The profile should only be modified before a simulation so the outputs are at initial values:
	burnt.insert(burnt.end(), numNewLayers, false);
	heatSink.insert(heatSink.end(), numNewLayers, 0.0);
	heatSource.insert(heatSource.end(), numNewLayers, 0.0);

	type.insert(type.end(), numNewLayers, type[numLayers - 1]);//Copy the label.

	numLayers += numNewLayers;
}

/** Insert new blank layers into the profile. (Properties must be added afterwards.)
 *
 * @param[in] insertAt The index to insert at (i.e. the layer currently at index insertAt will be
 *                     shifted down).
 * @param[in] numNewLayers The number of layers to insert.
 *
 * @returns Nothing.
 */
void GFProfile::InsertBlankLayers(const int insertAt, const int numNewLayers)
{
	if (numNewLayers < 1)
	{
		Stop("InsertBlankLayers(): Invalid value for numNewLayers: " + std::to_string(numNewLayers));
	}

	thickness_cm.insert(thickness_cm.begin() + insertAt, numNewLayers, ProfileUndef);//Or add values?
	layerDepth.insert(layerDepth.begin() + insertAt, numNewLayers, ProfileUndef);//Or add values?

	//Insert and set physical properties to invalid values:
	tempC.insert(tempC.begin() + insertAt, numNewLayers, ProfileUndef);
	t_ig.insert(t_ig.begin() + insertAt, numNewLayers, ProfileUndef);
	bulkDensity.insert(bulkDensity.begin() + insertAt, numNewLayers, ProfileUndef);
	inorganicPct.insert(inorganicPct.begin() + insertAt, numNewLayers, ProfileUndef);
	moistureContentPct.insert(moistureContentPct.begin() + insertAt, numNewLayers, ProfileUndef);
	c_s.insert(c_s.begin() + insertAt, numNewLayers, ProfileUndef);

	//The profile should only be modified before a simulation so the outputs are at initial values:
	burnt.insert(burnt.begin() + insertAt, numNewLayers, false);
	heatSink.insert(heatSink.begin() + insertAt, numNewLayers, 0.0);
	heatSource.insert(heatSource.begin() + insertAt, numNewLayers, 0.0);

	type.insert(type.begin() + insertAt, numNewLayers, std::string());//Empty string.

	numLayers += numNewLayers;
}

/** Revise the profile so that the top layer starts at a depth of zero.
 * 
 * If a soil profile is being translated from another representation where the top layer does not
 * start at depth zero this function can be called to adjust the layer depths so it does.  This
 * should be called before Interpolate() if needed.
 *
 * @returns Nothing.
 */
void GFProfile::Resurface()//Or AdjustDepths()
{
	if (layerDepth[0] != 0.0)
	{
		Msg.Log("The top layer is not at the surface. Adjusting layer depths.");

		double deltaZ = layerDepth[0];
		for (int i = 0; i < numLayers; i++)
		{
			layerDepth[i] -= deltaZ;
		}
	}
}

/** Perform some checks on the profile to see if looks valid.
 *
 * @param[in] uniformLayers For a profile is ready to be used for a ground fire simulation with
 *                          DownwardGroundFire() it should have layers of equal thickness.  Setting
 *                          this to true will add a check for this.
 *
 * @returns Whether the profile passed the checks.
 *
 * @note This code treats a profile without layers as invalid.  This is appropriate in some cases
 * but in others layers can be dropped during interpolation.  We need to add a way to distinguish
 * between the two.
 */
bool GFProfile::Validate(const bool uniformLayers) const
{
	bool valid = true;

	if (numLayers < 1)
	{
		Warning("Invalid value of numLayer.");
		valid = false;
	}

	//Make sure the number of property elements are consistant:
	if (numLayers != thickness_cm.size() ||
	    !SameLengths(thickness_cm, layerDepth, tempC, t_ig, bulkDensity, inorganicPct,
	                 moistureContentPct, c_s, burnt, heatSink, heatSource, type))
	{
		Warning("Some profile elements are not the same length.");
		valid = false;
	}

	//Confirm the depth profile is continuous:
	for (int i = 0; i < numLayers; i++)
	{
		if (thickness_cm[i] <= 0)
		{
			Warning("Invalid thickness for layer index " + std::to_string(i));
			valid = false;
		}

		if (i == 0)
		{
			if (layerDepth[i] != 0.0)
			{
				Warning("The top layer is not at the surface!");
				valid = false;
			}
		}
		else
		{
			//if (! FloatCompare(layerDepth[i], (layerDepth[i - 1] + thickness_cm[i - 1])))
			//Allow a tenth of a millimeter difference.  This is needed due to rounding issues
			//when output from Print() is reimorted.
			if (! FloatCompare(layerDepth[i], (layerDepth[i - 1] + thickness_cm[i - 1]), 0.02))
			{
				Warning("Layer dimensions not consistant for layer index " + std::to_string(i));
				valid = false;
			}
		}
	}

	if (uniformLayers && !EqualThickness())
	{
		Warning("Layers do not all have uniform thickness.");
		valid = false;
	}

	//Check for invalid values:
	if (!InRange(tempC, -60.0, 50.0))//A broad guess range for non-fire temperatures.
	{
		Warning("Unlikely temperature value(s).");
		valid = false;
	}

	if (!InRange(t_ig, 200, 300))//This range is rough and might be a bit wider.
	{
		Warning("Unlikely temperature value(s).");
		valid = false;
	}

	//Histisols can have low bulk densities approaching zero for some peats.  Compacted glacial till
	//on the other hand approaches that of concrete to we need to give a pretty wide range.  Values
	//pulled from figure 4.44 in:
	//Brady, N.C. and Weil, R.R., 2016.
	//The nature and properties of soils (Vol. 15, pg. 163). Columbus, OH: Pearson.
	if (!InRange(bulkDensity, 0.0, 2400.0))
	{
		Warning("Unlikely or impossible bulk density value(s).");
	}

	if (!InRange(inorganicPct, 0.0, 100.0))
	{
		Warning("Invalid percentage(s) for inorganic content.");
		valid = false;
	}

	//if (!InRange(moistureContentPct, 0.0, 1000.0))//This upper limit is poorly informed..
	//{
	//	Warning("Invalid soil moisture content(s).");
	//	valid = false;
	//}
	for (int i = 0; i < numLayers; i++)
	{
		if (moistureContentPct[i] < 0.0)
		{
			Warning("Invalid soil moisture content: " + std::to_string(moistureContentPct[i]));
			valid = false;
		}
		else if (moistureContentPct[i] > 1000.0)//Only warn until a better upper limit is determined:
		{
			Warning("High soil moisture content: " + std::to_string(moistureContentPct[i]));
		}
	}

	//We don't check the calculated properties.

	return valid;
}

/** Check if the soll layers are all the same thickness.
 *
 * @returns Whether the soll layers are equal thickness.
 */
bool GFProfile::EqualThickness() const
{
	if (numLayers == 0)
	{
		Stop("EqualThickness(): The profile has no layers.");
	}
	else if (numLayers > 1)
	{
		for (int l = 1; l < numLayers; l++)
		{
			if (!FloatCompare(thickness_cm[l], thickness_cm[0]))//It would be good if we didn't have to use this.
			{
				return false;
			}
		}
	}
	//If there is only one layer then we return true.

	return true;
}

/** Get the thickness of the whole soil profile.
 *
 * @returns The total thickness of the soil profile (cm).
 */ 
double GFProfile::GetThickness() const
{
	if (numLayers == 0)
	{
		return 0;
	}
	else
	{
		return std::accumulate(thickness_cm.begin(), thickness_cm.end(), 0.0);
		//Or: layerDepth[numLayers] + thickness_cm[numLayers]
	}
}

/** Get the simulated depth of burn, AKA burn thickness.
 *
 * @returns Burn depth (cm).
 * 
 * @note We could add checking that the a simulation has been completed successfully.  If it hasn't
 * a nonsense value will be returned.
 */ 
double GFProfile::GetBurnDepth() const
{
	double burnDepth = 0.0;

	//Find the lowest layer that burned and calculate the depth at the layer bottom:
	for (int l = numLayers - 1; l >= 0; l--)
	{
		if (burnt[l])
		{
			burnDepth = layerDepth[l] + thickness_cm[l];
			break;
		}
	}

	return burnDepth;
}

/** Print the soil profile data to an output stream.
 *
 * @param[in] output The output stream to print to.
 *
 * @returns The ostream so it can be concatenated to.
 */
std::ostream& GFProfile::Print(std::ostream& output) const
{
	output << "Ground fire soil profile:" << std::endl;

	//Print layer properties in table form:
	const int layerWidth = 6;
	const int thickWidth = 13;
	const int depthWidth = 11;
	const int tempWidth = 7;
	const int igWidth = 14;
	const int bdWidth = 13;
	const int inorgWidth = 13;
	const int mcWidth = 19;
	const int csWidth = 14;
	const int burntWidth = 6;
	const int sinkWidth = 10;
	const int srcWidth = 12;

	//If type labels are not provided for any of the layers drop the column:
	int typeWidth = 0;
	for (int i = 0; i < numLayers; i++)
	{
		typeWidth = std::max<int>(typeWidth, type[i].size());
	}
	if (typeWidth)
	{
		typeWidth += 1;
	}

	//Member name header:
	output << std::setw(layerWidth) << " "
		<< std::setw(thickWidth) << "thickness_cm"
		<< std::setw(depthWidth) << "layerDepth"
		<< std::setw(tempWidth) << "TempC"
		<< std::setw(igWidth) << "t_ig"
		<< std::setw(bdWidth) << "bulkDensity"
		<< std::setw(inorgWidth) << "inorganicPct"
		<< std::setw(mcWidth) << "moistureContentPct"
		<< std::setw(csWidth) << "c_s"
		<< std::setw(burntWidth) << "Burnt"
		<< std::setw(sinkWidth) << "heatSink"
		<< std::setw(srcWidth) << "heatSource" << std::endl;

	//Descriptive header:
	output << std::setw(layerWidth) << "Layer"
		<< std::setw(thickWidth) << "Thickness"
		<< std::setw(depthWidth) << "Depth"
		<< std::setw(tempWidth) << "Temp"
		<< std::setw(igWidth) << "Ignition Temp"
		<< std::setw(bdWidth) << "Bulk Density"
		<< std::setw(inorgWidth) << "Inorganic"
		<< std::setw(mcWidth) << "Moisture Content"
		<< std::setw(csWidth) << "Heat Capacity"
		<< std::setw(burntWidth) << "Burnt"
		<< std::setw(sinkWidth) << "Heat Sink"
		<< std::setw(srcWidth) << "Heat Source";

	if (typeWidth)
	{
		output << std::setw(typeWidth) << "Type";
	}
	output<< std::endl;

	//Units header:
	output << std::setw(layerWidth) << "#"
		<< std::setw(thickWidth) << "cm"
		<< std::setw(depthWidth) << "cm"
		<< std::setw(tempWidth) << "C"
		<< std::setw(igWidth) << "C"
		<< std::setw(bdWidth) << "kg/m^3"
		<< std::setw(inorgWidth) << "%"
		<< std::setw(mcWidth) << "%"
		<< std::setw(csWidth) << "kJ/kg/K"
		<< std::setw(burntWidth) << " "
		<< std::setw(sinkWidth) << "kJ/kg"
		<< std::setw(srcWidth) << "kJ/kg";

	if (typeWidth)
	{
		output << std::setw(typeWidth) << " ";
	}
	output << std::endl;

	//Values:
	output << std::boolalpha;//Print burnt as strings.
	for (int l = 0; l < numLayers; l++)
	{
		output << std::setw(layerWidth) << l + 1
			<< std::setw(thickWidth) << std::fixed << std::setprecision(2) << thickness_cm[l]
			<< std::setw(depthWidth) << std::fixed << std::setprecision(2) << layerDepth[l]
			<< std::setw(tempWidth) << std::fixed << std::setprecision(2) << tempC[l]
			<< std::setw(igWidth) << std::fixed << std::setprecision(1) << t_ig[l]
			<< std::setw(bdWidth) << std::fixed << std::setprecision(1) << bulkDensity[l]
			<< std::setw(inorgWidth) << std::fixed << std::setprecision(2) << inorganicPct[l]
			<< std::setw(mcWidth) << std::fixed << std::setprecision(2) << moistureContentPct[l]
			<< std::setw(csWidth) << std::fixed << std::setprecision(3) << c_s[l]
			<< std::setw(burntWidth) << burnt[l]
			<< std::setw(sinkWidth) << std::fixed << std::setprecision(2) << heatSink[l]
			<< std::setw(srcWidth) << std::fixed << std::setprecision(2) << heatSource[l];

		if (typeWidth)
		{
			output << std::setw(typeWidth) << type[l];
		}
		output << std::endl;
	}
	output << std::noboolalpha;
	output.copyfmt(std::ios(nullptr));//Restore the previous print settings.

	return output;
}

/** Print the soil profile to an output stream as a set of delimited data rows for each layer,
 * suitable for data ingestion.
 *
 * @param[in] output The output stream to print to.
 * @param[in] delim The delimiter character.  Defaults to the tab character.
 *
 * @returns The ostream so it can be concatenated to.
 */
std::ostream& GFProfile::PrintDelimited(std::ostream& output, const char delim) const
{
	//Print the header:
	output << "Layer" << delim << "LayerThickness" << delim << "LayerDepth" << delim <<
	          "TempC" << delim << "t_ig" << delim << "BulkDensity" << delim <<
	          "InorganicPct" << delim << "MoistureContentPct" << delim << "c_s" << delim <<
	          "Burnt" << delim << "HeatSink" << delim << "HeatSource" << delim << "Type" <<
	          std::endl;//Or newline?
	//Omit the soil mass?

	//Print the layer data in rows:
	output << std::boolalpha;//Print burnt as strings.
	//Print so numbers are round-trip guaranteed and can be reingested faithfully:
	output << std::setprecision(std::numeric_limits<double>::max_digits10);
	for (int l = 0; l < numLayers; l++)
	{
		output << (l + 1) << delim << thickness_cm[l] << delim << layerDepth[l] << delim <<
		          tempC[l] << delim << t_ig[l] << delim << bulkDensity[l] << delim <<
		          inorganicPct[l] << delim << moistureContentPct[l] << delim << c_s[l] << delim <<
		          burnt[l] << delim << heatSink[l] << delim << heatSource[l] << delim << type[l] <<
		          std::endl;//Or newline?
	}
	output << std::noboolalpha;
	output.copyfmt(std::ios(nullptr));//Restore the previous print settings.

	return output;
}

/** Set the print mode to use with the stream operator. 
 *
 * @param mode Pass 0 or greater for a human readable table or a negative value for tab delimited.
 *
 * @returns Nothing.
 *
 * @note We may make sybolic constants for the value of mode.  In futre might use postive values to
 * aler the number of digits?
 */ 
void GFProfile::SetPrintMode(const int mode)
{
	printMode = mode;
}

//Private Functions:--------------------------------------------------------------------------------

/** Remove a layers from the profile.
 *
 * @param[in] layerIndex The layer index to remove.
 *
 * @returns Nothing.
 *
 * @note Very thin layers may be removed during interpolation.  This leaves the the profile in an
 * intermediate state with gaps in the profile.  The profile will remain invalid until interpolation
 * is complete.
 */
void GFProfile::DropLayer(const int layerIndex)
{
	if (layerIndex < 0 || layerIndex >= numLayers)
	{
		Stop("DropLayer(): Invalid layer index: " + std::to_string(layerIndex));
	}

	thickness_cm.erase(thickness_cm.begin() + layerIndex);
	layerDepth.erase(layerDepth.begin() + layerIndex);
	tempC.erase(tempC.begin() + layerIndex);
	t_ig.erase(t_ig.begin() + layerIndex);
	bulkDensity.erase(bulkDensity.begin() + layerIndex);
	inorganicPct.erase(inorganicPct.begin() + layerIndex);
	moistureContentPct.erase(moistureContentPct.begin() + layerIndex);
	c_s.erase(c_s.begin() + layerIndex);
	burnt.erase(burnt.begin() + layerIndex);
	heatSink.erase(heatSink.begin() + layerIndex);
	heatSource.erase(heatSource.begin() + layerIndex);
	type.erase(type.begin() + layerIndex);

	numLayers -= 1;
}

//External functions:-------------------------------------------------------------------------------

/* Overloaded stream print operator for GFProfile.
 *
 */
std::ostream& operator<<(std::ostream& output, const GFProfile& gfProfile)
{
	if (gfProfile.printMode >= 0)
	{
		gfProfile.Print(output);
	}
	else
	{
		gfProfile.PrintDelimited(output);
	}
	return output;
}
