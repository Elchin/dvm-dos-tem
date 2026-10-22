/*! \file
 *
 */
#ifndef SNOWLAYER_H_
#define SNOWLAYER_H_
#include "Layer.h"
#include "physicalconst.h"
#include "parameters.h"

#include <string>
#include <cmath>
using namespace std;

class SnowLayer: public Layer {
public:
  SnowLayer();
  ~SnowLayer();

  void clone(SnowLayer* sl); /*assign same member to another layer*/

  void updateThick();
  void updateDensity(snwpar_dim *snwpar);

  // Option 4 (v5j) SEB: model_settings.seb_snow_tc, 0 = Goodrich
  // 1982 k = 2.9e-6 rho^2 (legacy, default); 1 = Sturm et al. 1997 (k(rho),
  // 156-600 kg/m3; 0.023+0.234 rho below).
  static int snow_tc_form;
  // Refs: Sturm et al. 1997, doi:10.3189/S0022143000002781; legacy Goodrich 1982, doi:10.1139/t82-047.
  virtual double getFrzThermCond();//get frozen layer thermal conductivity
  virtual double getUnfThermCond();//get unfrozen layer thermal conductivity

  virtual double getFrzVolHeatCapa();//get frozen layer specific heat capcity
  virtual double getUnfVolHeatCapa();//get unfrozen layer specific heat capacity

  virtual double getMixVolHeatCapa(); //Yuan:

private:

  double getThermCond5Sturm();
  double getThermCond();

};
#endif /*SNOWLAYER_H_*/
