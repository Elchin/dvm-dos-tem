/***********************************************************************************************//**
 * @file GroundFire.h
 * \author Joshua M. Rady
 * Woodwell Climate Research Center
 * \date 2025
 *
 * @brief This header file declares functions to simulate smoldering ground fire using a moderate
 * complexity model where combustion is predicted layer by layer.
 *
 **************************************************************************************************/

#ifndef GROUNDFIRE_H
#define GROUNDFIRE_H

#include <vector>

#include "GFProfile.h"

double SoilHeatOfCombustion(const double inorganicPct, const double organicHOC = 21.6);//WILDFIRE_PARAMETER
double EffectiveSmolderingHOC(const double inorganicPct);
double SoilHeatSink(const double moistureContentPct, const double t_int, const double t_ig, const double c_s);
std::vector <double> HeatDistributionLinear(const int numLayers);
double DownwardGroundFire(GFProfile& soilCol, const double fireHeatInput = 0.0,//With no heat input nothing will burn.
                          const double heatLossFactor = 0.83, const double surfacePD = 1.0,
                          const double smolderPD = 5.0, const int surfaceTM = 1,
                          const int smolderTM = 1);

#endif //GROUNDFIRE_H
