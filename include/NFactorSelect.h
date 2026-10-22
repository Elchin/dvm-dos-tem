#ifndef NFACTORSELECT_H_
#define NFACTORSELECT_H_
// Daily n-factor selection used by Cohort::updateMonthly_Env (pure helpers,
// factored out unchanged so they can be unit tested; see tests/unit).
namespace nfactor_sel {

// Value of the n-factor selector before the first day of each month has set it
// (only its sign matters: > 0 selects the summer n-factor). Fix (b) (v5i).
const double TDRV_MONTH_START = 1.0;

// Air temperature used to choose the n-factor for today:
//   sameday == 1 -> today's Tair; otherwise (legacy) the previous day's Tair
//   (tdrv_prev; TDRV_MONTH_START on the first day of a month).
inline double selector_tair(int sameday, double tair_today, double tdrv_prev) {
  return (sameday == 1) ? tair_today : tdrv_prev;
}

// Summer n-factor unless the selector temperature is freezing (<= 0 C).
inline double select_nfactor(double tsel, double nfactor_s, double nfactor_w) {
  double nf = nfactor_s;
  if (tsel <= 0.0) {
    nf = nfactor_w;
  }
  return nf;
}

} // namespace nfactor_sel
#endif /*NFACTORSELECT_H_*/
