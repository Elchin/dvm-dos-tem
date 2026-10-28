#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

/* \file
 * This class is used to implement
 * Runge-Kutta Method for integration of TEM equations
 * (see ATPress42005a, Numerical Recipes in C++,
 *  The art of scientific computing, 2nd edition
 *  W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery,
 *  2005, Cambridge University Press)
 *  Page 712-727
 *
 * Low-order classical Runge-Kutta Formulas with step size control
 *   and their application to some heat transfer problems
 * by Erwin Feblberg
 * NASA TR R-315, 1969
 * REPORT TL 521 A3312 no.315 in GI-IARC 1st floor
 * see also http://www.ece.uwaterloo.ca/~ece204/TheBook/14IVPs/rkf45/
 */

#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;

#include "Soil_Bgc.h"
#include "Vegetation_Bgc.h"
#include "temconst.h"

class Integrator {
public :
  Integrator();
  ~Integrator();

  // After change the keys here, you must remember to update the
  // NUMEQ_VEG and NUM_VEG_STATE in temconst.h !!!
  enum vegvarkey {
    I_VEGC = 0,              // index starts from zero
    I_STRN = NUM_PFT_PART,
    I_LABN = 2*NUM_PFT_PART,
    I_DEADC, I_DEADN,        // 0 ~ NUM_VEG_STATE-1 for veg state variables

    // NUM_VEG_STATE is the starting index for flux variables
    I_INGPP = NUM_VEG_STATE,
    I_INNPP = NUM_VEG_STATE + 1 * NUM_PFT_PART,
    I_GPP   = NUM_VEG_STATE + 2 * NUM_PFT_PART,
    I_NPP   = NUM_VEG_STATE + 3 * NUM_PFT_PART,
    I_RM    = NUM_VEG_STATE + 4 * NUM_PFT_PART,
    I_RG    = NUM_VEG_STATE + 5 * NUM_PFT_PART,
    I_LTRC  = NUM_VEG_STATE + 6 * NUM_PFT_PART,
    I_SNUP  = NUM_VEG_STATE + 7 * NUM_PFT_PART,
    I_NMBOL = NUM_VEG_STATE + 8 * NUM_PFT_PART,
    I_NRSRB = NUM_VEG_STATE + 9 * NUM_PFT_PART,

    //11*NUM_PFT_PART veg C/N flux variables
    I_LTRN  = NUM_VEG_STATE + 10 * NUM_PFT_PART,

    //because indexed from zero, so here is the sum of all above
    I_INNUP = NUM_VEG_STATE + 11 * NUM_PFT_PART,
    I_LNUP
  };

  // After change the keys here, you must remember to update the
  // NUMEQ_SOI and MAXSTATE in temconst.h !!
  enum soivarkey {
    I_L_RAWC  = 0,
    I_L_SOMA  = MAX_SOI_LAY,
    I_L_SOMPR = 2 * MAX_SOI_LAY,
    I_L_SOMCR = 3 * MAX_SOI_LAY, // (4*LAYERS) these 4 lines -
                                 // Layered Soil C state variables
    I_L_ORGN  = 4 * MAX_SOI_LAY,
    I_L_AVLN  = 5 * MAX_SOI_LAY, // (4*LAYERS) these 2 lines -
                                 // Layered Soil N state variables

    I_WDEBRISC= 6 * MAX_SOI_LAY, // because indexed from zero, so here
                                 // is the sum of all above
    I_WDEBRISN,

    I_L_RH_RAW  = NUM_SOI_STATE,
    I_L_RH_SOMA = NUM_SOI_STATE + MAX_SOI_LAY,
    I_L_RH_SOMPR= NUM_SOI_STATE + 2 * MAX_SOI_LAY,
    I_L_RH_SOMCR= NUM_SOI_STATE + 3 * MAX_SOI_LAY,
    I_L_NMIN    = NUM_SOI_STATE + 4 * MAX_SOI_LAY,
    // 6*MAX_SOI_LAY soil C/N flux variables
    I_L_NIMMOB  = NUM_SOI_STATE + 5*MAX_SOI_LAY,

    I_RH_WD = NUM_SOI_STATE + 6 * MAX_SOI_LAY, // Because indexed from zero,
                                             // so here is the sum of all above

    I_ORGNLOSS, I_AVLNLOSS

  };


  char predstr_veg[NUMEQ_VEG][11];
  char predstr_soi[NUMEQ_SOI][11];

  bool vegbgc;
  bool soibgc;

  void setBgcData(BgcData *bdp);
  void setSoil_Bgc(Soil_Bgc * soib);
  void setVegetation_Bgc(Vegetation_Bgc * vegb);

  void updateMonthlyVbgc();
  void updateMonthlySbgc(const int &numsl);

  void report_array(float array[]);

  // Integrator safety net (soil): when adapt() gives up ("blackhole": the
  // state is reset to the start of the run of tiny steps at fraction t0 of the
  // month), the soil inputs for the rest of the month (litterfall C/N, woody
  // debris C/N, N deposition) are still delivered so that C is conserved.
  int blackhole_count;    // cumulative number of adapt() give-ups
  int safetynet_count;    // cumulative number of soil safety-net deliveries
  bool safetynet_fired;   // set by the last updateMonthlySbgc()
  // Fix C1 (v5k), OFF by default: model_settings.soil_nmin_tol_fix. The soil
  // step-size control (boundcon) tests the error of the monthly net-N-
  // mineralisation accumulator NMIN against ptol*|NMIN|. NMIN is a signed net
  // flux; when mineralisation ~ immobilisation it passes through ~0, every
  // doubled step is rejected and adapt() crawls at dt ~ 2^-19 (264k steps in
  // TVC_B Warming nfeed Nov 2012 with B1+B2+B4), and the float state then loses
  // the tiny C decrements to rounding (+2.74 gC unbooked). 1 = scale the NMIN
  // tolerance by the gross flux/pool size ptol*max(|NMIN|, |NIMMOB|, |AVLN|).
  static int nmin_tol_fix;
  // Refs (C1): no published source for the fix; standard restored: step-size error scaled by a magnitude that cannot pass through zero (Press et al. 2007, Numerical Recipes 3rd ed., sec. 17.2; Hairer, Norsett & Wanner 1993, doi:10.1007/978-3-540-78862-1).
  // Fix L1 (v5m), OFF by default: model_settings.n_safetynet_conserve.
  // When adapt() gives up, the soil state is reset to month fraction t0 and
  // applySoilSafetyNet() delivers the remaining inputs, but the root N
  // extraction for the remaining (1 - t0) of the month is not applied, although
  // the vegetation integrator has already added the WHOLE month's N uptake to
  // the plants (Cohort: veg is integrated before soil). Each event therefore
  // creates (1 - t0) * sum(rtnextract) g N/m2 (N leak; verified: N closure is
  // exactly 0 in every year without an event). 1 = debit that N from the soil:
  // first from ORGN of the layers the roots extracted from (in proportion to
  // rtnextract, capped at each layer's ORGN), then from ORGN of any layer (in
  // proportion to the remaining ORGN), then from AVLN; anything still unpaid
  // (all donor pools empty) is carried and debited at the start of the next
  // soil month. Also delivers the remaining N deposition exactly where
  // Soil_Bgc::deltastate() would (rooted top layer, or any top layer when
  // n_rootless_fix == 2). Total N is then conserved to float round-off.
  // Carbon: no leak (the veg does not remove soil C; all litter C delivered).
  static int n_safetynet_conserve;
  // Refs (L1): no published source for the fix (code: Integrator::applySoilSafetyNet, be8a62fa); principle restored: closed N mass balance of the TEM / DOS-TEM N cycle, plant N uptake debited from the soil N pools (McGuire et al. 1992, doi:10.1029/92GB00219; Raich et al. 1991, doi:10.2307/1941899; Euskirchen et al. 2009, doi:10.1890/08-0806.1; Yuan et al. 2012, doi:10.1890/11-1957.1); donor = ORGN first: boreal/arctic plants take up organic N (Nasholm et al. 1998, doi:10.1038/31921; Schimel & Bennett 2004, doi:10.1890/03-8002); clipping/resetting a stiff ODE state is not mass-conserving unless the removed mass is booked (Shampine et al. 2005, doi:10.1016/j.amc.2004.12.011); unclosed N biases long C-N spin-ups (Thornton & Rosenbloom 2005, doi:10.1016/j.ecolmodel.2005.04.008).
  // Fix G2 (v5n), OFF by default: model_settings.soil_giveup_retry.
  // 1 = event-only retry: the soil month is first integrated exactly as in
  // legacy; only if adapt() gives up (100 steps forced at the minimum step
  // size), the month is integrated again from its start with the G1 extraction
  // limiter (Soil_Bgc::nx_limiter_active) for this month only. If the retry
  // gives up too, the safety net is applied as before. Months without a
  // give-up are bit-identical to legacy (same code path, same arithmetic).
  static int soil_giveup_retry;
  // Refs (G1/G2): plant N uptake is a saturating (Michaelis-Menten) function of soil available N in TEM (Raich et al. 1991, doi:10.2307/1941899; McGuire et al. 1992, doi:10.1029/92GB00219); demand-based uptake is down-regulated when soil mineral N cannot supply it (Thornton et al. 2007, doi:10.1029/2006GB002868); arctic/boreal plants also take up organic N (Nasholm et al. 1998, doi:10.1038/31921; Schimel & Bennett 2004, doi:10.1890/03-8002); a zero-order sink drives a non-negative ODE state below zero in finite time, which no step-size reduction can avoid, positivity needs a state-dependent sink (Shampine et al. 2005, doi:10.1016/j.amc.2004.12.011; Burchard et al. 2003, doi:10.1016/S0168-9274(03)00101-6).
  // Diagnostics (v5n), OFF by default: model_settings.soil_giveup_diag.
  // 1 = record why adapt() gave up (first step forced at dt_min: the negative
  // pool or the error-test variable, with the soil N state of that layer) and
  // count rkf45 calls; reported on stderr only (outputs unchanged).
  static int soil_giveup_diag;
  int retry_count;          // G2: soil months retried with the limiter
  int retry_fail_count;     // G2: retries that still gave up
  bool retry_fired;         // G2: set by the last updateMonthlySbgc()
  int forced_months;        // diag: months with >= 1 step forced at dt_min
  long long rkf45_calls;    // diag: cumulative rkf45 calls (soil integrator)
  int last_negpool;         // diag: last pool index rkf45 found negative
  int diag_test;            // diag: step test before forcing (0 = negative pool, >1 = error var)
  int diag_pool;            // diag: failing y index (negative pool or error var)
  float diag_t0, diag_dt;   // diag: month fraction and dt at the first forced step
  float diag_y[NUMEQ];      // diag: state at the first forced step
  double diag_rtn[MAX_SOI_LAY], diag_netnmin[MAX_SOI_LAY], diag_nimmob[MAX_SOI_LAY];
  double diag_avlnlost, diag_avln0[MAX_SOI_LAY];
  int last_reject_test, last_reject_pool; // diag: last rejected step of this month (test, y index)
  int diag_prevtest, diag_prevpool;       // diag: copies at the first forced step
  double diag_nscale, diag_totnetnmin, diag_totnextract, diag_totdzavln; // diag: Soil_Bgc::deltan budget
  double safetynet_ndebit;  // L1: N debited from soil by the last event (g/m2)
  double safetynet_nneed;   // L1: N the last event had to debit (g/m2)
  double n_debt;            // L1: N still owed (donors empty), paid next month
  double n_debit_total;     // L1: cumulative N debited (g/m2), diagnostic
  // L1 helper (pure function, unit-tested): remove up to 'need' (>0) from
  // pool[0..n-1], first in proportion to weights w (if w != NULL; only entries
  // with w > 0 and pool > 0), then the remainder in proportion to the
  // remaining positive pool. Pools never go below 0; negative pools are left
  // unchanged. Returns the amount actually removed (sum of float changes).
  static double debitProportional(double need, float pool[], const int n,
                                  const double w[]);
  float safetynet_t0;     // month fraction already integrated when it gave up
  double safetynet_c;     // C delivered by the last safety net (g/m2)
  double safetynet_n;     // N delivered by the last safety net (g/m2)

private:

  void applySoilSafetyNet();
  void recordGiveupDiag(const float pstate[], const int numeq, const int ptest,
                        const float ptime, const float pdt);  // v5n diagnostics
  double conserveSafetyNetN(const double frac);  // L1 (v5m)
  void payNDebt();                                 // L1 (v5m)
  float bh_t0;            // time at which adapt() saved oldstate

  int numsl;  // actual number of soil layers

  void c2ystate_veg(float y[]);
  void c2ystate_soi(float y[]);
  int adapt(float pstate[], const int &numeq);
  void y2cstate_veg(float y[]);
  void y2cstate_soi(float y[]);
  void y2cflux_veg(float y[]);
  void y2cflux_soi(float y[]);

  bool rkf45(const int& numeq, float pstate[], float& pdt);
  void delta(float pstate[], float pdstate[]);

  void dc2yflux_veg(float pdstate[]);
  void dc2yflux_soi(float pdstate[]);
  void dc2ystate_veg(float pdstate[]);
  void dc2ystate_soi(float pdstate[]);
  void y2tcstate_veg(float pdstate[]);// t stands for temporary
  void y2tcstate_soi(float pdstate[]);// t stands for temporary

  int checkPools();
  void step(const int& numeq, float pstate[], float pdstate[],
            float ptstate[], float& pdt);
  int boundcon( float ptstate[], float err[], float& ptol) ;

  int blackhol;
  float inittol;
  int maxit;
  int maxitmon;
  int retry;

  int syint; // << -- Appears to be unused??
  int test;

  float y[NUMEQ];
  float dum4[NUMEQ];
  float error[NUMEQ];

  float dum5[NUMEQ];
  float dumy[NUMEQ];

  float ydum[NUMEQ];
  float yprime[NUMEQ];
  float rk45[NUMEQ];

  float f11[NUMEQ];
  float f2[NUMEQ];
  float f13[NUMEQ];
  float f3[NUMEQ];
  float f14[NUMEQ];
  float f4[NUMEQ];
  float f15[NUMEQ];
  float f5[NUMEQ];
  float f16[NUMEQ];
  float f6[NUMEQ];

  static float  a1;
  static float  a3, a31, a32;
  static float  a4, a41, a42, a43;
  static float  a5, a51, a52, a53, a54;
  static float  b1, b3, b4, b5;
  static float  b6, b61, b62, b63, b64, b65;

  Soil_Bgc * ssl;
  Vegetation_Bgc * veg;
  BgcData * bd;

};
#endif /*INTEGRATOR_H_*/


