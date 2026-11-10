/*! \file SnowClass.h
 * model_settings.snow_scheme = 1: Sturm et al. (1995) snow classes as implemented in SiBCASA by
 * Schaefer et al. (2009), JGR 114, F02021, doi:10.1029/2008JF001125 (Sec. 2.2-2.4, eqs 1-7, Table 2),
 * with Anderson (1976) thermal-aging + overburden compaction in the CLM form
 * (Oleson et al. 2010, CLM4.0 tech note NCAR/TN-478+STR, Sec. 7.2.6, eqs 7.40-7.47).
 * snow_scheme = 0 (default) = legacy dvmdostem snow (bit-identical; nothing here is called).
 *
 * Pure functions (no model state) so they can be unit tested; Ground calls them.
 */
#ifndef SNOWCLASS_H_
#define SNOWCLASS_H_
#include <string>
#include <vector>

struct SnowClassPar {
  std::string name;        // tundra, taiga, maritime, ephemeral, prairie, alpine
  double fbotmax = 0.;     // max bottom-layer fraction (-)                     Table 2
  double Dmin = 0.;        // depth where a bottom layer starts to form (m)     Table 2
  double Dmax = 0.;        // depth where fbot reaches fbotmax (m)              Table 2
  double rho_topref = 0.;  // top-layer reference density (kg/m3)               Table 2
  double rho_botref = 0.;  // bottom-layer reference density (kg/m3)            Table 2
  double rho_bulkobs = 0.; // observed bulk density, reference month (kg/m3)    Table 2
  double sigma_bulkobs = 1.; // its standard deviation (kg/m3)                 Table 2
  double k_dh = -1.;       // fixed depth-hoar conductivity of the bottom layer (W/m/K); <0 = none.
                           //   tundra 0.18, taiga 0.072 (Sturm et al. 1997, as Schaefer 2009 Sec 2.3 [19])
  bool wind = false;       // wind term drho_w in new-snow density (tundra, prairie; eq 6)
};

namespace snowclass {
const double RHO_MIN = 50.0;    // minimum new-snow density (kg/m3), eq 4
const double RHO_CAP = 600.0;   // upper bound on any snow layer density (kg/m3) (validity of Sturm 1997 k)

// eqs 4-7: new-snow density from air temperature (C) and wind speed (m/s)
double drhoT(double tair);
double drhoW(const SnowClassPar& p, double u);
double rhoNew(const SnowClassPar& p, double tair, double u);
// eqs 1-2: bottom-layer fraction vs bulk depth D (m)
double fbot(const SnowClassPar& p, double D);
// eq 3: density difference between the layers (kg/m3, >= 0)
double drho(const SnowClassPar& p, double D, double rho_bulk);
// two-layer profile that conserves the bulk mass: rho_top*(1-fb)*D + rho_bot*fb*D = rho_bulk*D.
// sign of (rho_top - rho_bot) follows (rho_topref - rho_botref). drho is reduced if needed so both
// densities stay >= RHO_MIN (and <= RHO_CAP).
void twoLayer(const SnowClassPar& p, double D, double rho_bulk, double& fb, double& rho_top, double& rho_bot);
// Anderson (1976) / CLM4 compaction rate of one layer (1/s, <= 0): thermal aging CR1 + overburden CR2.
//   rho: layer density (kg/m3); tem: layer T (C); burden: mass above the layer mid-point (kg/m2);
//   wet: liquid water present (c2 = 2 in CR1)
double compactionRate(double rho, double tem, double burden, bool wet);
// new thickness after dt seconds (hourly exponential sub-steps; mass constant, rho capped at RHO_CAP)
double compactDz(double dz, double ice, double tem, double burden, bool wet, double dt);
// Conservative remap in mass coordinate: given the old layer masses m_old (top->bottom) and temperatures,
// and the new masses m_new (same sum), return new temperatures (mass-weighted averages of the old
// material in each new layer). Conserves sum(m*T) to round-off.
std::vector<double> remapT(const std::vector<double>& m_old, const std::vector<double>& t_old,
                           const std::vector<double>& m_new);
// Target masses for layers of thickness dz (top->bottom) under the two-layer profile; sum = total mass
// exactly (residual put in the largest layer). Also returns, per layer, the fraction of its thickness
// in the bottom (depth-hoar) layer.
void targetMasses(const SnowClassPar& p, const std::vector<double>& dz, double M,
                  std::vector<double>& m_new, std::vector<double>& fdh);
// Table 2 defaults for a class name; false if unknown
bool defaultClass(const std::string& name, SnowClassPar& p);
// Load the class table + CMT->class map from a JSON file. cmtcode like "CMT04". override (if non-empty)
// forces the class. Returns false and sets err on failure.
bool load(const std::string& file, const std::string& cmtcode, const std::string& override_class,
          SnowClassPar& out, std::string& err);
// series (thickness-weighted harmonic) mix of depth-hoar k and the density-based k for a layer whose
// thickness fraction fdh lies in the depth-hoar layer
double seriesK(double fdh, double k_dh, double k_rho);
}
#endif
