/*! \file EnthalpySolver.h
 * Implicit, globally convergent enthalpy solver (Stage B library; used by the model
 * through Soil_Env_enthalpy.cpp when model_settings.thermal_solver = 1, Stage C) for 1-D snow/soil/rock heat conduction with phase change.
 * Design: /workspace/tem-dev/enthalpy_solver_design.md (enthalpy/stageB).
 *
 * Discretisation: cell-centred finite volumes (one cell per dvmdostem Layer),
 * backward Euler, H_i(T_i^{n+1}) - H_i^n = dt*(face fluxes + S_i), i.e.
 * H(T) + A T = b with A a symmetric tridiagonal M-matrix (Tubini et al. 2021,
 * Eqs. 12-16, doi:10.5194/tc-15-2541-2021).
 * Nonlinear solvers:
 *  - nested Newton (Casulli & Zanolli 2010, doi:10.1137/100786320; Tubini et
 *    al. 2021 Sect. 3.2 / App. B): Jordan split of C_a = dH/dT at T* = 0 degC.
 *  - Katzenelson piecewise-affine solver (Hoetten et al. 2026, ETNA 65,
 *    doi:10.1553/etna_vol65s553) for the free-water curve (sfcc 0), used as a
 *    cross-check.
 * Enthalpy reference: liquid water at 0 degC (H = 0), the same definition as
 * Soil_Env::layerEnthalpy (energy_diag).
 */
#ifndef ENTHALPYSOLVER_H_
#define ENTHALPYSOLVER_H_

#include <vector>

namespace enthalpy {

enum CellKind { SOIL = 0, SNOW = 1, ROCK = 2 };

// Freezing characteristic
//  0: free water, piecewise-affine H(T): phi = phimin below -eps, linear H on [-eps, 0]
//  1: Schaefer & Jafarov 2016 power law, identical to Soil_Bgc::sjLayerPhi for
//     T < -eps (max with phimin = liqmin/W), linear phi ramp on [-eps, 0]
enum Sfcc { SFCC_FREE = 0, SFCC_SJ2016 = 1 };

struct Cell {
  int kind;        // CellKind
  double dz;       // m
  double Cs;       // solid heat capacity per area, J m-2 K-1 (0 for snow)
  double W;        // total water (liq + ice), kg m-2
  double liqmin;   // residual liquid water, kg m-2 (legacy minliq)
  int sfcc;        // Sfcc (snow/rock: forced free water)
  double eps;      // ramp width, degC (> 0)
  double clay, silt, sand; // texture as passed to sjLayerPhi (percent or fraction)
  double forg;     // organic fraction for sjLayerPhi (1 = organic)
  double lam_u;    // unfrozen conductivity W m-1 K-1
  double lam_f;    // frozen conductivity W m-1 K-1
  Cell();
};

// ---- freezing curve / enthalpy of one cell (pure functions) ----
double sjPhi(double t, double forg, double fc, double fsi, double fsa); // == Soil_Bgc::sjLayerPhi
double phimin(const Cell& c);               // residual liquid fraction
double H(const Cell& c, double T);          // J m-2
double Ca(const Cell& c, double T);         // dH/dT, J m-2 K-1 (left derivative at kinks)
double liq(const Cell& c, double T);        // kg m-2, consistent with H(T)
double invertH(const Cell& c, double h);    // T with H(T) = h (monotone; bisection-safeguarded Newton)
double conductivity(const Cell& c, double T); // lam_u^(1-f) lam_f^f, f = ice/W
bool hasLatent(const Cell& c);

// ---- Stage C helpers (pure functions, used by Soil_Env_enthalpy.cpp) ----
// General enthalpy of a cell state that may be off the curve (T, liq, ice), J m-2.
double Hstate(const Cell& c, double T, double liq, double ice);
// True if (T, liq) is on the curve of c (|liq - liq(c,T)| <= tol_kg).
bool onCurve(const Cell& c, double T, double liqv, double tol_kg);
// Put a cell on its curve keeping its enthalpy h (J m-2): returns T, sets liq/ice (sum c.W).
double equilibrate(const Cell& c, double h, double& liqv, double& icev);
// Merge two cells conservatively: temperature of the merged cell for summed heat capacities
// C = Cs + liq c_w + ice c_i and summed liq/ice (keeps Hstate exactly).
double mergeT(double Cu, double Tu, double Cl, double Tl);

// Front reconstruction from a soil profile (top to bottom; soil cells only).
struct FrontIn {
  std::vector<double> z, dz, T, ice, ice_eps, W;  // layer top (m), thickness, temperature, ice,
                                                   // ice at T = -eps (curve), total water
  std::vector<int> latent;                         // 1 if the cell has latent capacity
  double Tabove;   // temperature above the first soil cell (snow cell or surface)
  double Tbelow;   // temperature below the last soil cell (rock)
  double eps;      // ramp width (degC)
  int maxfronts;   // MAX_NUM_FNT
};
struct FrontOut {
  std::vector<double> fz;     // front depths, top to bottom (m)
  std::vector<int> ftype;     // 1 = freezing front (frozen above), -1 = thawing front (thawed above)
  std::vector<int> frozen;    // per cell: 1 frozen, -1 thawed, 0 holds a front
  std::vector<double> frozenfrac; // per cell: frozen thickness fraction
};
// Cells are classified frozen (T < -eps), thawed (no ice) or ramp (latent heat partly released).
// A run of ramp cells holds one front (or a lens pair) placed so that the frozen thickness equals
// sum_k f_k dz_k, f_k = ice_k / ice_k(-eps); a 0 degC crossing between two non-ramp neighbours is
// placed by linear interpolation between their mid-points. Fronts are capped at maxfronts by
// removing the closest adjacent pairs. Cell states follow the legacy front convention.
void reconstructFronts(const FrontIn& in, FrontOut& out);

struct BC {
  int top_type;    // 0 Dirichlet T_s at the surface, 1 Robin: conductance cn_top to T_s
  double Ts;       // degC
  double cn_top;   // W m-2 K-1 (legacy virtual layer: 100)
  double q_bot;    // bottom heat flux INTO the column, W m-2 (0 = legacy zero flux)
  BC();
};

struct Options {
  int face_mode;   // 0 harmonic (series) mean [default], 1 Tubini max rule (fallback)
  double tol;      // max |cell residual|, J m-2
  int max_outer, max_inner;
  int max_katz;    // max Katzenelson segment switches
  Options();
};

struct Stats {
  int outer, inner;   // NCZ outer iterations, inner Newton iterations (= tridiagonal solves)
  int solves;         // tridiagonal solves (NCZ: inner; Katzenelson: linear solves)
  bool converged;
  double res_max;     // max |cell residual| at exit, J m-2
  double dH;          // sum H^{n+1} - H^n, J m-2
  double gtop;        // top flux into column over the step, J m-2
  double gbot;        // bottom flux into column over the step, J m-2
  double src;         // sum of sources over the step, J m-2
  double closure;     // dH - (gtop + gbot + src), J m-2
  Stats();
};

class ColumnSolver {
public:
  Options opt;
  ColumnSolver();
  // One backward-Euler step of length dt [s]. T: in = T^n, out = T^{n+1}.
  // S: optional per-cell source W m-2 (NULL = none). Returns stats.converged.
  bool stepNCZ(const std::vector<Cell>& cells, std::vector<double>& T, double dt,
               const BC& bc, const std::vector<double>* S, Stats& st);
  // Katzenelson on the same discrete system; requires every latent cell sfcc 0.
  bool stepKatz(const std::vector<Cell>& cells, std::vector<double>& T, double dt,
                const BC& bc, const std::vector<double>* S, Stats& st);
  // face conductances of the last assembled step (W m-2 K-1); K[0] = top
  const std::vector<double>& faces() const { return K; }

private:
  std::vector<double> K, lam, lo, di, up, rhs, b, Hn, w1, w2, x, v, r;
  int n;
  void assemble(const std::vector<Cell>& cells, const std::vector<double>& T,
                double dt, const BC& bc, const std::vector<double>* S);
  void applyA(const std::vector<double>& T, std::vector<double>& out) const; // A T
  void thomas(const std::vector<double>& d, std::vector<double>& sol);       // (A + diag d) sol = rhs
  void finish(const std::vector<Cell>& cells, const std::vector<double>& T,
              double dt, const BC& bc, const std::vector<double>* S, Stats& st);
  double rfloor(const std::vector<double>& AT, const std::vector<double>& T) const;
  double dt_, ktop_;
};

} // namespace enthalpy

#endif
