/*! \file
 *
 */
#ifndef SOILLAYER_H_
#define SOILLAYER_H_
#include "Layer.h"

#include <math.h>
#include <cmath>
#include <memory>
using namespace std;

class SoilLayer:public Layer {
public:

  SoilLayer();
  virtual ~SoilLayer();

  virtual double getFrzThermCond();//get frozen layer thermal conductivity
  virtual double getUnfThermCond();//get unfrozen layer thermal conductivity
  virtual double getFrzVolHeatCapa();//get frozen layer specific heat capcity
  virtual double getUnfVolHeatCapa();//get unfrozen layer specific heat capacity
  virtual double getMixVolHeatCapa(); //Yuan: all soil components

  double getAlbedoVis();// get albedo of visible radiation
  double getAlbedoNir();// get albedo of Nir radiation

  void derivePhysicalProperty();

  double getMatricPotential();
  double getHydraulicCond();

  // DRAFT (tvc-summer): Kersten number form for unfrozen soil.
  // 0 (default) = legacy Ke = ln(Sr) + 1 (Ke = 0 for Sr <= 1/e = 0.37);
  // 1 = Johansen (1975) fine-grained / CLM / Lawrence & Slater (2008) form
  //     Ke = log10(Sr) + 1 (Ke = 0 only for Sr <= 0.1).
  static int kersten_form;

};
#endif /*SOILLAYER_H_*/
