#ifndef TKNETC_H_
#define TKNETC_H_
/** TKNETC: combined total ecosystem carbon loss diagnostic (gC/m2/month,
 *  positive = carbon leaving the ecosystem).
 *
 *  Exact definition implemented here and in Runner::output_netCDF:
 *
 *    TKNETC = NEE_month + TKLATC + TKRESPC - TKREGROWC   (TKREGROWC only if regrow on, else 0)
 *
 *  where NEE_month is the model's monthly net ecosystem exchange in the
 *  loss-positive sign convention used by the RTS test closures
 *  (tests/rts/test_rts_ch4.py, test_rts_regression.py), i.e. the negative of
 *  the cal-json NEE (GPP - respiration):
 *
 *    NEE_month = RG + RM + soil_h Heterotrophic_CO2_and_CH4 + RHDWD - GPP
 *
 *  CH4 handling (thermokarst_ch4, scheme a in RTSMethane.h): soil C leaving
 *  as CO2/CH4 is NOT RHSOM + TKCH4EMIS. In Cohort::updateMonthly_RTSMethane
 *    TKRHCO2   = RHSOM - prod + oxid
 *    TKCH4EMIS = prod - oxid          (emis = prod - oxid, RTSMethane::compute)
 *  so TKRHCO2 + TKCH4EMIS = RHSOM exactly (algebraically). When CH4 is on
 *  this helper therefore uses (TKRHCO2 + TKCH4EMIS) as the soil term: CH4
 *  carbon is counted exactly once, neither double-counted (RHSOM+TKCH4EMIS
 *  would count the emitted CH4 C twice, once inside RHSOM and once as CH4)
 *  nor dropped (TKRHCO2 alone would omit it). When CH4 is off the soil term
 *  is RHSOM (TKRHCO2/TKCH4EMIS are then zero diagnostics and must not be
 *  used). TKLATCH4 is only a diagnostic split of TKLATC and is never added.
 *  TKCH4POOL is a state (0 in scheme a), part of ECOC in closures, not a flux.
 *
 *  With ECOC = VEGC+SOC+DEADC+DWDC+TKCH4POOL this gives, month by month,
 *    dECOC = -TKNETC   (up to the suite's <0.01 gC/m2/month closure),
 *  i.e. cumulative TKNETC = -(ECOC_end - ECOC_start).
 */
struct TKNetC {
  static double soilCarbonLoss(bool ch4_on, double rhsom, double rh_co2, double ch4_emis) {
    return ch4_on ? (rh_co2 + ch4_emis) : rhsom;
  }
  static double compute(bool ch4_on, bool regrow_on,
                        double gpp, double rg, double rm,
                        double rhsom, double rhdwd, double rh_co2, double ch4_emis,
                        double latc, double respc, double regrowc) {
    const double nee_month = rg + rm + soilCarbonLoss(ch4_on, rhsom, rh_co2, ch4_emis) + rhdwd - gpp;
    return nee_month + latc + respc - (regrow_on ? regrowc : 0.0);
  }
};
#endif
