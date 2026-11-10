/*! \file SnowClass.cpp  -- see include/SnowClass.h (model_settings.snow_scheme = 1). */
#include "../include/SnowClass.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <json/json.h>

namespace snowclass {

// eq 5 (Judson & Doesken 2000 fit): 0 below -15 C, 1.7 (Ta+15)^1.5 at or above.
double drhoT(double tair) {
  if (!(tair >= -15.0)) { return 0.0; }
  return 1.7 * pow(tair + 15.0, 1.5);
}

// eq 6 (Sato et al. 2008): tundra and prairie only; 0 below 2 m/s, linear above.
double drhoW(const SnowClassPar& p, double u) {
  if (!p.wind || !(u >= 2.0)) { return 0.0; }
  return (u - 2.0) * (p.rho_topref - RHO_MIN) / 3.0;
}

// eqs 4 and 7
double rhoNew(const SnowClassPar& p, double tair, double u) {
  double r = RHO_MIN + drhoT(tair) + drhoW(p, u);
  r = fmin(r, p.rho_topref);
  return fmax(r, RHO_MIN);
}

// eqs 1-2
double fbot(const SnowClassPar& p, double D) {
  if (!(D > 0.0) || p.fbotmax <= 0.0) { return 0.0; }
  const double Dhalf = 0.5 * (p.Dmax + p.Dmin);
  const double Dslope = 10.0 / (p.Dmax - p.Dmin);
  return p.fbotmax / (1.0 + exp(Dslope * (Dhalf - D)));
}

// eq 3
double drho(const SnowClassPar& p, double D, double rho_bulk) {
  if (p.fbotmax <= 0.0) { return 0.0; }
  const double s = (rho_bulk - p.rho_bulkobs) / p.sigma_bulkobs;
  return fabs(p.rho_topref - p.rho_botref) * (fbot(p, D) / p.fbotmax) * exp(-s * s);
}

void twoLayer(const SnowClassPar& p, double D, double rho_bulk, double& fb, double& rho_top, double& rho_bot) {
  fb = fbot(p, D);
  double dr = drho(p, D, rho_bulk);
  const double sgn = (p.rho_topref >= p.rho_botref) ? 1.0 : -1.0;
  // rho_top = rb + sgn*dr*fb ; rho_bot = rb - sgn*dr*(1-fb)  (bulk mass conserved for any dr)
  double lo, hi;   // limits so that both stay within [RHO_MIN, RHO_CAP]
  if (sgn > 0) { lo = (1.0 - fb) > 0 ? (rho_bulk - RHO_MIN) / (1.0 - fb) : dr; hi = fb > 0 ? (RHO_CAP - rho_bulk) / fb : dr; }
  else         { lo = fb > 0 ? (rho_bulk - RHO_MIN) / fb : dr; hi = (1.0 - fb) > 0 ? (RHO_CAP - rho_bulk) / (1.0 - fb) : dr; }
  dr = fmax(0.0, fmin(dr, fmin(lo, hi)));
  rho_top = rho_bulk + sgn * dr * fb;
  rho_bot = rho_bulk - sgn * dr * (1.0 - fb);
}

// Anderson (1976) in the CLM4.0 form (Oleson et al. 2010, eqs 7.42-7.46):
//   CR1 = -c3 c1 c2 exp(-c4 (Tf - T)), c3 = 2.777e-6 1/s, c4 = 0.04 1/K,
//         c1 = 1 (rho <= 100) else exp(-0.046 (rho - 100)), c2 = 2 if wet else 1;
//   CR2 = -Ps / eta, eta = eta0 exp(c5 (Tf - T) + c6 rho), eta0 = 9e5 kg s m-2, c5 = 0.08 1/K,
//         c6 = 0.023 m3/kg, Ps = mass above the layer + half of the layer's own mass (kg m-2; as in
//         CLM SnowHydrologyMod, Ps in kg m-2 with eta0 in kg s m-2).
double compactionRate(double rho, double tem, double burden, bool wet) {
  const double td = fmax(0.0, -tem);            // Tf - T (K), >= 0
  double c1 = (rho > 100.0) ? exp(-0.046 * (rho - 100.0)) : 1.0;
  const double c2 = wet ? 2.0 : 1.0;
  const double cr1 = -2.777e-6 * c1 * c2 * exp(-0.04 * td);
  const double eta = 9.0e5 * exp(0.08 * td + 0.023 * rho);
  const double cr2 = -fmax(burden, 0.0) / eta;
  return cr1 + cr2;
}

double compactDz(double dz, double ice, double tem, double burden, bool wet, double dt) {
  if (!(dz > 0.0) || !(ice > 0.0)) { return dz; }
  const int nsub = 24; const double h = dt / nsub;
  double d = dz;
  const double dmin = ice / RHO_CAP;
  for (int k = 0; k < nsub; ++k) {
    const double rho = ice / d;
    if (rho >= RHO_CAP) { break; }
    d *= exp(compactionRate(rho, tem, burden + 0.5 * ice, wet) * h);
    if (d < dmin) { d = dmin; break; }
  }
  return d;
}

std::vector<double> remapT(const std::vector<double>& m_old, const std::vector<double>& t_old,
                           const std::vector<double>& m_new) {
  const size_t n = m_old.size(), nn = m_new.size();
  std::vector<double> t(nn, 0.0);
  std::vector<double> co(n + 1, 0.0); for (size_t k = 0; k < n; ++k) { co[k + 1] = co[k] + m_old[k]; }
  double a = 0.0;
  for (size_t j = 0; j < nn; ++j) {
    const double b = a + m_new[j];
    double h = 0.0, w = 0.0;
    for (size_t k = 0; k < n; ++k) {
      const double lo = fmax(a, co[k]), hi = fmin(b, co[k + 1]);
      if (hi > lo) { h += (hi - lo) * t_old[k]; w += hi - lo; }
    }
    if (w > 0.0) { t[j] = h / w; }
    else { t[j] = (j < n) ? t_old[j] : (n ? t_old[n - 1] : 0.0); }
    a = b;
  }
  return t;
}

void targetMasses(const SnowClassPar& p, const std::vector<double>& dz, double M,
                  std::vector<double>& m_new, std::vector<double>& fdh) {
  const size_t n = dz.size();
  m_new.assign(n, 0.0); fdh.assign(n, 0.0);
  double D = 0.0; for (double d : dz) { D += d; }
  if (!(D > 0.0) || !(M > 0.0)) { return; }
  double fb, rt, rb;
  twoLayer(p, D, M / D, fb, rt, rb);
  const double zi = (1.0 - fb) * D;   // depth from top of the top/bottom interface
  double a = 0.0, sum = 0.0; size_t imax = 0;
  for (size_t i = 0; i < n; ++i) {
    const double b = a + dz[i];
    const double top = fmax(0.0, fmin(b, zi) - a);
    const double bot = dz[i] - top;
    m_new[i] = rt * top + rb * bot;
    fdh[i] = (dz[i] > 0.0) ? fmax(0.0, fmin(1.0, bot / dz[i])) : 0.0;
    sum += m_new[i];
    if (m_new[i] > m_new[imax]) { imax = i; }
    a = b;
  }
  m_new[imax] += M - sum;   // round-off residual: SWE conserved exactly
}

double seriesK(double fdh, double k_dh, double k_rho) {
  if (!(fdh > 0.0) || !(k_dh > 0.0)) { return k_rho; }
  if (fdh >= 1.0) { return k_dh; }
  return 1.0 / (fdh / k_dh + (1.0 - fdh) / k_rho);
}

bool defaultClass(const std::string& name, SnowClassPar& p) {
  // Schaefer et al. 2009 Table 2; k_dh Sec 2.3 [19]; wind term Sec 2.3 [20]
  struct R { const char* n; double fbm, dmin, dmax, rt, rb, ro, so, kdh; bool w; };
  static const R t[] = {
    {"tundra",    0.3, 0.0, 0.1, 350, 100, 261.108, 68.206, 0.18,  true},
    {"taiga",     0.7, 0.0, 0.7, 300, 200, 213.363, 49.462, 0.072, false},
    {"maritime",  0.8, 0.2, 0.8, 200, 300, 286.501, 80.979, -1.0,  false},
    {"ephemeral", 0.8, 0.2, 0.8, 250, 350, 321.289, 82.321, -1.0,  false},
    {"prairie",   0.8, 0.2, 0.8, 350, 250, 272.423, 61.057, -1.0,  true},
    {"alpine",    0.8, 0.2, 0.7, 200, 300, 267.796, 56.85,  -1.0,  false}};
  for (const R& r : t) {
    if (name == r.n) {
      p.name = r.n; p.fbotmax = r.fbm; p.Dmin = r.dmin; p.Dmax = r.dmax; p.rho_topref = r.rt;
      p.rho_botref = r.rb; p.rho_bulkobs = r.ro; p.sigma_bulkobs = r.so; p.k_dh = r.kdh; p.wind = r.w;
      return true;
    }
  }
  return false;
}

bool load(const std::string& file, const std::string& cmtcode, const std::string& override_class,
          SnowClassPar& out, std::string& err) {
  std::ifstream f(file.c_str());
  if (!f) { err = "cannot open snow class file " + file; return false; }
  Json::Value root; Json::CharReaderBuilder b; std::string e;
  if (!Json::parseFromStream(b, f, &root, &e)) { err = "parse error in " + file + ": " + e; return false; }
  std::string cls = override_class;
  if (cls.empty()) {
    const Json::Value& m = root["cmt_class"];
    if (m.isObject() && m.isMember(cmtcode)) { cls = m[cmtcode].asString(); }
    else { cls = root.get("default_class", "").asString(); }
  }
  if (cls.empty()) { err = "no snow class for " + cmtcode + " in " + file; return false; }
  SnowClassPar p;
  if (!defaultClass(cls, p)) { p = SnowClassPar(); p.name = cls; }
  const Json::Value& c = root["classes"][cls];
  if (!c.isObject() && p.fbotmax == 0.0) { err = "unknown snow class '" + cls + "'"; return false; }
  if (c.isObject()) {
    p.fbotmax = c.get("fbotmax", p.fbotmax).asDouble();
    p.Dmin = c.get("Dmin", p.Dmin).asDouble();
    p.Dmax = c.get("Dmax", p.Dmax).asDouble();
    p.rho_topref = c.get("rho_topref", p.rho_topref).asDouble();
    p.rho_botref = c.get("rho_botref", p.rho_botref).asDouble();
    p.rho_bulkobs = c.get("rho_bulkobs", p.rho_bulkobs).asDouble();
    p.sigma_bulkobs = c.get("sigma_bulkobs", p.sigma_bulkobs).asDouble();
    p.k_dh = c.get("k_dh", p.k_dh).asDouble();
    p.wind = c.get("wind", p.wind).asBool();
  }
  if (!(p.Dmax > p.Dmin) || !(p.sigma_bulkobs > 0.0) || p.fbotmax < 0.0 || p.fbotmax > 1.0) {
    err = "invalid parameters for snow class '" + cls + "'"; return false;
  }
  out = p;
  return true;
}

} // namespace snowclass
