#ifndef DAILYWATEROUT_H_
#define DAILYWATEROUT_H_
// v5p-daily-water-fix (output only): daily NetCDF values of TRANSPIRATION, EET and PET.
//
// The monthly outputs are cohort totals accumulated from the daily cohort ('edall') values:
//   EnvData::atm_endOfDay  m_l2a.eet += d_l2a.eet; m_l2a.pet += d_l2a.pet  (also stores daily_eet/daily_pet)
//   EnvData::veg_endOfDay  m_v2a.tran += d_v2a.tran
// Before this fix the daily outputs did not follow that bookkeeping:
//   * EET/PET daily = sum over ALL PFTs of ed[ip].daily_eet/daily_pet. Each per-PFT value
//     (Cohort::getEd4land_daily) already contains the whole-cohort ground terms (soil
//     evaporation, snow sublimation, copied to every PFT by assignGroundEd2pfts_daily), so the
//     ground part was counted once per active PFT; PFTs with vegcov = 0 contributed stale values.
//   * TRANSPIRATION daily = one value per month (edall->d_v2a.tran of the month's last day) and
//     one value written per year.
// Now the daily values are the cohort ('edall') daily values, so the days of a month sum to the
// monthly output exactly. Nothing here is read by the model; no state or flux changes.
// Refs: no published source; rationale: a daily output of an accumulated flux must be the term
//   that is summed into the monthly output (output bookkeeping only).
#include <vector>
#include "EnvData.h"

namespace dwout {

// store the day's cohort transpiration (call after the day's last change of d_v2a.tran,
// i.e. together with EnvData::veg_endOfDay, which adds the same value to m_v2a.tran)
inline void record_tran(EnvData& e, const int dayidx) {
  if (dayidx >= 0 && dayidx < 31) { e.daily_tran[dayidx] = e.d_v2a.tran; }
}

// append the first dinm days of a month's daily array to an output buffer
inline void append_days(const double (&arr)[31], const int dinm, std::vector<double>& out) {
  const int n = (dinm < 0) ? 0 : ((dinm > 31) ? 31 : dinm);
  for (int id = 0; id < n; ++id) { out.push_back(arr[id]); }
}

} // namespace dwout

#endif /* DAILYWATEROUT_H_ */
