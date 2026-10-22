// test_access.h - include model headers with private/protected members
// accessible to the unit tests.
//
// The tests need to set up and inspect objects whose state is private
// (e.g. Climate::derive_daily_secondary_vars, Stefan::rebuildFrontsFromLayerStates,
// Integrator::applySoilSafetyNet). Rather than adding accessors to the model,
// the test translation units re-define the access keywords before including
// the model headers. This does not change any class layout (GCC lays out
// members in declaration order regardless of access) and the model object
// files are compiled normally, so the model code under test is unchanged.
//
// All system/third-party headers that the model headers pull in must be
// included BEFORE the redefinition.
#ifndef UNIT_TEST_ACCESS_H_
#define UNIT_TEST_ACCESS_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <netcdf.h>
#include <json/json.h>
#include <boost/filesystem.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include "../../include/TEMLogger.h"

#define private public
#define protected public
#include "../../include/timeconst.h"
#include "../../include/layerconst.h"
#include "../../include/physicalconst.h"
#include "../../include/Climate.h"
#include "../../include/TEMUtilityFunctions.h"
#include "../../include/Layer.h"
#include "../../include/SoilLayer.h"
#include "../../include/SnowLayer.h"
#include "../../include/MineralLayer.h"
#include "../../include/Ground.h"
#include "../../include/Stefan.h"
#include "../../include/Richards.h"
#include "../../include/EnvData.h"
#include "../../include/BgcData.h"
#include "../../include/CohortData.h"
#include "../../include/Snow_Env.h"
#include "../../include/Soil_Env.h"
#include "../../include/Soil_Bgc.h"
#include "../../include/Integrator.h"
#include "../../include/Vegetation_Env.h"
#include "../../include/ModelData.h"
#undef private
#undef protected

// Free functions defined in Climate.cpp without a header declaration.
float calculate_vpd(const float svp, const float vp);
float calculate_saturated_vapor_pressure(const float tair);
float calculate_clouds(const float girr, const float nirr);
float calculate_par(const float clds, const float nirr);
float calculate_girr(const float lat, const int im);
std::vector<float> calculate_girr_daily(const float lat);
std::vector<float> calculate_daily_prec(const int midx, const float mta, const float mprec);

#endif
