#ifndef RHMOISTFROZEN_OVERRIDE_H_
#define RHMOISTFROZEN_OVERRIDE_H_
// model_settings.rhmoistfrozen_override: optional per-CMT override of the cmt_bgcsoil.txt parameter
// rhmoistfrozen (frozen-soil RH moisture multiplier, Soil_Bgc bgcpar.rhmoistfrozen), without editing parameters/.
//   absent, null or {}            -> no override (default; bit-identical)
//   {"70": 0.0, "CMT66": 0.05}    -> CMT70 uses 0.0, CMT66 uses 0.05; every other CMT keeps its file value
// Keys: CMT number as a string, digits only ("70", "7", "07") or with the "CMT" prefix ("CMT70"); 0..99.
// Values: finite numbers in [0, 1]. Anything else (non-object, bad key, duplicate CMT, non-numeric or
// out-of-range value) is FATAL in ModelData. Applied overrides are logged on stderr per cohort initialisation.
// Motivation: upstream efa2118a changed CMT70 rhmoistfrozen 0.0 -> 1.0, which respires the Gydan permafrost
// mineral C during equilibrium (114 -> 7 kg/m2); see docs/WORK_NOTES_rhmoistfrozen_override.md.
#include <json/json.h>
#include <map>
#include <string>
#include <cmath>
#include <cctype>

namespace rmfo {

// Parse a CMT key ("70", "07", "CMT70") into 0..99; false if malformed.
inline bool parse_key(const std::string& k, int& cmt) {
  std::string d = k;
  if (d.size() > 3 && d.compare(0, 3, "CMT") == 0) d = d.substr(3);
  if (d.empty() || d.size() > 2) return false;
  for (char c : d) if (!std::isdigit(static_cast<unsigned char>(c))) return false;
  cmt = std::stoi(d);
  return cmt >= 0 && cmt <= 99;
}

// Parse model_settings["rhmoistfrozen_override"] (pass the value, or a null Json::Value when absent).
// On success fills out (possibly empty) and returns true; on error sets err and returns false.
inline bool parse(const Json::Value& v, std::map<int, double>& out, std::string& err) {
  out.clear();
  if (v.isNull()) return true;
  if (!v.isObject()) { err = "model_settings.rhmoistfrozen_override must be an object {\"<CMT>\": value}"; return false; }
  for (const std::string& k : v.getMemberNames()) {
    int cmt = -1;
    if (!parse_key(k, cmt)) { err = "model_settings.rhmoistfrozen_override: bad CMT key \"" + k + "\" (use e.g. \"70\" or \"CMT70\")"; return false; }
    const Json::Value& x = v[k];
    if (!x.isNumeric() || x.isBool()) { err = "model_settings.rhmoistfrozen_override[\"" + k + "\"] must be a number"; return false; }
    const double d = x.asDouble();
    if (!std::isfinite(d) || d < 0.0 || d > 1.0) { err = "model_settings.rhmoistfrozen_override[\"" + k + "\"] must be in [0, 1]"; return false; }
    if (out.count(cmt)) { err = "model_settings.rhmoistfrozen_override: CMT " + std::to_string(cmt) + " given twice"; return false; }
    out[cmt] = d;
  }
  return true;
}

// Override value for this CMT number, if any.
inline bool lookup(const std::map<int, double>& m, int cmt, double& v) {
  auto it = m.find(cmt);
  if (it == m.end()) return false;
  v = it->second; return true;
}

} // namespace rmfo
#endif
