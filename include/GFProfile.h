/***********************************************************************************************//**
 * @file GFProfile.h
 * \author Joshua M. Rady
 * Woodwell Climate Research Center
 * \date 2025
 *
 * @brief This header file declares an object used to represent the manipulate the soil column used
 * for simulations of smoldering ground fire.
 *
 **************************************************************************************************/

#ifndef GFPROFILE_H
#define GFPROFILE_H

#include <iostream>//Or just <ostream>?
#include <string>
#include <vector>

/** GFProfile Constants
 *
 * ProfileUndef is used to indicate an undefined/empty/missing value for a soil profile property.
 * 0 or -1 would be a fine value for most properties, since it is imposible and therefore clearly
 * invalid.  We use -300 because it is an imposible celcius temperature.
 */
const double ProfileUndef = -300.0;

/* @class GFProfile
 * @brief A data object representing an organic soil profile:
 */
class GFProfile {
	public:
	
	//Column properties: none. 
	
	//Layer properties:
	//The ground fire calculation currently expects equal layer thicknesses but incoming profiles
	//may have layers of varying thicknesses that will then be interpolated.
	std::vector <double> thickness_cm;//Layer thickness (cm).
	std::vector <double> layerDepth;//Depth at top of layer (cm).
	
	std::vector <double> tempC;//Layer temperature in Celcius.
	std::vector <double> t_ig;//The soil temperature of ignition (C).
	std::vector <double> bulkDensity;//Dry soil mass per volume (kg/m^3).
	std::vector <double> inorganicPct;//Percent inorganic content on a dry basis (~ ash content).
	//Soil mostures content is the water mass per volume / dry soil mass per volume:
	std::vector <double> moistureContentPct;//(water mass / layer) / (dry soil mass / layer) * 100%
	std::vector <double> c_s;//Layer heat capacity (kJ/kg/K)
	
	//Calculated ground fire output values:
	std::vector <bool> burnt;
	//The heat source and sink values can be helpful for understanding the results so we store them:
	std::vector <double> heatSink;//kJ/kg
	std::vector <double> heatSource;//kJ/kg

	//An optional field for a layer description, name, or label.  It is not used in calculations and
	//is mainly for print outs.  Not all layers need to be filled out.
	std::vector <std::string> type;

	//Constructors:
	GFProfile();
	GFProfile(const int numLayers);

	int NumLayers() const;
	double DrySoilMassKg(const int layerIndex) const;

	void Interpolate(const double newLayerThickness);
	void InterpolateBetween(const int topLayer, int bottomLayer, const int numNewLayers);
	void AddLayersToTop(const int numNewLayers = 1);
	void AddLayersToBottom(const int numNewLayers = 1);
	void InsertBlankLayers(const int insertAt, int const numNewLayers);
	void Resurface();
	bool Validate(const bool uniformLayers = false) const;
	bool EqualThickness() const;
	double GetThickness() const;
	double GetBurnDepth() const;
	std::ostream& Print(std::ostream& output) const;
	std::ostream& PrintDelimited(std::ostream& output, const char delim = '\t') const;
	void SetPrintMode(const int mode);
	
	private:
	
	//Profile information:
	//The purpose of this variable is to enforce constancy in property lengths.  This makes it clear
	//how many layers there should be and helps avoid problems since the layer properties are
	//public.  To avoid problems the number of layers should only changed using the public
	//functions.  Use NumLayers() to get the count.
	int numLayers;

	int printMode;//Store the printing style for stream operator output.

	void DropLayer(const int layerIndex);

	friend std::ostream& operator<<(std::ostream& output, const GFProfile& gfProfile);
};

//External functions:
std::ostream& operator<<(std::ostream& output, const GFProfile& gfProfile);

#endif //GFPROFILE_H
