// nc_fixtures.h - tiny synthetic NetCDF files for the daily-climate loader
// and CO2-indexing tests (no real model inputs needed).
#ifndef UNIT_NC_FIXTURES_H_
#define UNIT_NC_FIXTURES_H_
#include <netcdf.h>
#include <string>
#include <vector>
#include <stdexcept>

namespace ut {

inline void ncchk(int st) { if (st != NC_NOERR) throw std::runtime_error(nc_strerror(st)); }

struct DailyFileSpec {
  int ndays = 730;
  std::string units = "days since 1950-01-01";   // first data year = 1950
  std::string calendar = "365_day";
  float tair = -5.0f, prec = 1.0f, nirr = 100.0f, vapo = 4.0f;
  int bad_index = -1;            // if >= 0: put bad_value into bad_var at this step
  std::string bad_var;
  float bad_value = 0.0f;
  float lat = 68.0f, lon = -133.0f;
};

// Write a daily projected-climate-like file with dims (time, Y, X) = (n, 1, 1).
inline void write_daily_climate(const std::string& fn, const DailyFileSpec& s) {
  int nc, dt, dy, dx, vt, vlat, vlon;
  ncchk(nc_create(fn.c_str(), NC_CLOBBER | NC_NETCDF4, &nc));
  ncchk(nc_def_dim(nc, "time", s.ndays, &dt));
  ncchk(nc_def_dim(nc, "Y", 1, &dy));
  ncchk(nc_def_dim(nc, "X", 1, &dx));
  ncchk(nc_def_var(nc, "time", NC_INT, 1, &dt, &vt));
  ncchk(nc_put_att_text(nc, vt, "units", s.units.size(), s.units.c_str()));
  ncchk(nc_put_att_text(nc, vt, "calendar", s.calendar.size(), s.calendar.c_str()));
  int d2[2] = {dy, dx};
  ncchk(nc_def_var(nc, "lat", NC_FLOAT, 2, d2, &vlat));
  ncchk(nc_def_var(nc, "lon", NC_FLOAT, 2, d2, &vlon));
  const char* names[4] = {"tair", "precip", "nirr", "vapor_press"};
  float vals[4] = {s.tair, s.prec, s.nirr, s.vapo};
  int vids[4]; int d3[3] = {dt, dy, dx};
  for (int k = 0; k < 4; ++k) ncchk(nc_def_var(nc, names[k], NC_FLOAT, 3, d3, &vids[k]));
  ncchk(nc_enddef(nc));
  std::vector<int> t(s.ndays);
  for (int i = 0; i < s.ndays; ++i) t[i] = i;
  ncchk(nc_put_var_int(nc, vt, t.data()));
  ncchk(nc_put_var_float(nc, vlat, &s.lat));
  ncchk(nc_put_var_float(nc, vlon, &s.lon));
  for (int k = 0; k < 4; ++k) {
    std::vector<float> v(s.ndays, vals[k]);
    // a weak seasonal signal so that the series is not constant
    if (k == 0) for (int i = 0; i < s.ndays; ++i) v[i] = s.tair + 15.0f * (float)std::sin(2.0 * M_PI * ((i % 365) - 110) / 365.0);
    if (s.bad_index >= 0 && s.bad_var == names[k]) v[s.bad_index] = s.bad_value;
    ncchk(nc_put_var_float(nc, vids[k], v.data()));
  }
  ncchk(nc_close(nc));
}

// Write a co2 file with dimension 'year' and variables 'year', 'co2'.
inline void write_co2(const std::string& fn, int y0, int n, bool gap = false) {
  int nc, dy, vy, vc;
  ncchk(nc_create(fn.c_str(), NC_CLOBBER | NC_NETCDF4, &nc));
  ncchk(nc_def_dim(nc, "year", n, &dy));
  ncchk(nc_def_var(nc, "year", NC_INT, 1, &dy, &vy));
  ncchk(nc_def_var(nc, "co2", NC_FLOAT, 1, &dy, &vc));
  ncchk(nc_enddef(nc));
  std::vector<int> yr(n); std::vector<float> c(n);
  for (int i = 0; i < n; ++i) { yr[i] = y0 + i + ((gap && i >= n / 2) ? 1 : 0); c[i] = 300.0f + (float)(yr[i] - 1900); }
  ncchk(nc_put_var_int(nc, vy, yr.data()));
  ncchk(nc_put_var_float(nc, vc, c.data()));
  ncchk(nc_close(nc));
}

inline std::string tmpname(const std::string& tag) {
  char buf[256];
  snprintf(buf, sizeof(buf), "/tmp/dvm_unittest_%d_%s.nc", (int)getpid(), tag.c_str());
  return buf;
}

} // namespace ut
#endif
