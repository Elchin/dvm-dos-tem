# T3: process-based MeHg in the dvmdostem Hg module: design (no code, no runs)

Branch feat/hg-t2 at 75fccfb1. Baseline: **B2F** (config/hg_t2_overlay_B2F.json = R1B + hg_t2_ev_diff 1, hg_t2_ev_diff_thaw 1,
hg_fd 0.215; with the MeHg-on preset hg_t2_mehg_preset 1 = M4V1). Scoring: WORK_NOTES §22-23, 9 scored bands 9/9, MeHg
primary 7/7 (TVC 0.751 %, HPC 0.093 %, ratio 8.08, Toolik fen 2.22 %; May-Oct MeHg export HPC 36.3, TVC 23.2, fen 40.6 ng/m2).
Written 2026-11-26. **No compute:** only code reading and a few seconds of post-processing of existing B2F outputs
(an2/B2F_*.json). 0 core-h; ledger unchanged at 12.8608. No commits, no pushes. Box copy: /workspace/hg/t3/.

## 0. Summary
- Today's MeHg is a **T3-lite index**: per layer, MeHg relaxes to P/kd with P = k_meth x layer THg x A x Q10 term, capped at
  the layer Hg. It is moved with water and exported (M2-M4), but **it is never taken out of Hg(II)**. Methylation and
  demethylation move no Hg mass, and every Hg(II) process (sorption, evasion, DOC-limited export, thaw release) still sees
  the full layer THg. The rate law has no explicit redox, sulfate or trophic-class term. Its only constraints are three
  %MeHg classes and an export band.
- **Recommended approach (option B below): MeHg as a mass-conserving sub-pool of the layer Hg.**
  1. THg = Hg(II) + MeHg in the existing pools. Methylation moves Hg(II) to MeHg, and demethylation moves it back to Hg(II)
     (oxidative, the default). An optional reductive share goes to Hg(0).
  2. The Hg(II) processes act on THg - MeHg (a "closure" switch).
  3. MeHg gets its own soil-water partitioning (Kd offset from Hg(II)) and its own export. It keeps today's exemption from
     the Hg:DOC cap.
  4. The rate law keeps the calibrated M4V1 driver (anaerobic C turnover x saturation x freeze factor) as its core. It
     adds multiplicative factors for trophic/fen class and sulfate that default to 1. Literature ranges are given for
     each, but none of them is fitted.
  5. Fen/wetland enters through the CMT. dvmdostem already has CMT06 wet sedge, CMT31 Boreal Bog, CMT32 Boreal Fen and a
     wetland CMT list. A class map in model_settings picks the driver and an optional fen factor. Landscape mixing is
     done in post-processing. There is no fen-fraction input.
- **No new prognostic state.** s.mehg is already in the restart (v4), and partitioning is equilibrium. Everything is
  default off and bit-identical when off. With hg_t3 1 at its defaults the run must reproduce B2F bitwise (an
  index-equivalence test). Only the closure, the Kd offset and the new factors change results.
- **Two read-only findings shape the design** (§3):
  - Gross methylation in the index is 8-28x the MeHg stock per year (TVC 1838 ug/m2/yr gross on about 93 ug/m2). So even
    a 0.1 % reductive demethylation share would add 1.8 ug/m2/yr Hg(0) at TVC, 16 % of its total evasion. Demethylation
    products must default to Hg(II).
  - The literature soil MeHg Kd is about 1 log unit below Hg(II) (Allison & Allison 2005 medians 2.8 vs 3.8). Applied
    literally, it could raise MeHg export up to about 10x and break the HPC/TVC export band. The Kd offset must be scanned
    and chosen by the export bands, not imposed.
- **Cost:** about 0.74 core-h expected, capped at 1.0, in 5 phases (§6). About 7-9 developer-days. The Toolik
  drainage-contrast test needs a copied run-mask (n=12 cells), which needs your OK.

## 1. Current state (include/HgTracer.h at 75fccfb1 unless noted)
| piece | where | what it does |
|---|---|---|
| index parameters | l.261-272 | hg_mehg, hg_mehg_driver; k_meth_ch4 / k_meth_proxy; k_demeth 3 /month at t_ref (about 0.1 /d); q10_meth 1, q10_demeth 2 |
| T2 MeHg switches | l.345-374, parse l.547-592, validation l.637-675 | t2_k_meth_proxy / _ch4; mehg_subst (M1); frz, frzd (M2/M3); transport, mobfac, lag (M2); satd, zrun, frzmob (M4); mehg_preset (M5) |
| MeHg-on preset (= M4V1) | `t2_mehg_preset_` l.451-453 | lag 1, transport 1, frz 1 (cutoff 0 C), frzd 1, satd 1, k_proxy 99.9, k_ch4 6079 |
| index update | `mehg_month` l.1892-1915 | P = km x substrate x max(0, A) x q10_meth^((T - t_ref)/10) x freeze factor; kd = k_demeth x q10_demeth^(...); MeHg <- Meq + (MeHg - Meq) e^-kd, Meq = P/kd; clamp to [0, layer Hg]. **No Hg moves between species.** |
| driver A | Cohort.cpp (a_turn) | CH4 production / layer C when cd.ch4_on (wetland CMT with the CH4 module), else the saturated fraction (daily water table with satd) x RH / layer C |
| MeHg transport | `diss_transport_` l.1201-1243; retention l.1131-1140; balance l.1267-1270 | each daily dissolved transfer F carries F x mobfac x MeHg / THg; exports remove MeHg (HG_MEHGEXP); HG_MEHGBAL closes |
| DOC cap | l.1220 | the Hg:DOC export cap (DOC step 2) applies to THg - MeHg: "MeHg keeps its Kd path" (already a partial species split) |
| geometry, fire | l.1578-1624, l.1749-1797, l.1853 | MeHg follows its layer's Hg losses proportionally and is remapped with the layers |
| restart | l.1919-1942 | s.mehg[MAX_SOI_LAY] is part of hg_state (v4) |
| wetland gate | Cohort.cpp l.96-99; ModelData.cpp l.974 | cd.is_wetland = CMT in wetland.cmts; ch4_on = is_wetland && wetland.ch4 |

Site set-up (read from cs2/B2F and sites2/B2F_*):
| case | CMT | drainage_class | wetland block | MeHg driver | target class |
|---|---|---|---|---|---|
| TVC | 51 (Tussock Tundra TVC) | 1 (poor) | absent (off) | proxy | bog 0.7 % (our assumption; lit_7questions Q1) |
| HPC | 01 (Boreal Black Spruce) | 1 (poor) | absent (off) | proxy | permafrost plateau 0.1 % |
| Toolik fen | 06 (Wet sedge tundra) | cell (0,0) = 1 | enabled, ch4 on, cmts [6, 31, 32, 55, 61, 75, 92, 93] | CH4 | fen 2.0 % |
| Toolik proxy | 06 | 1 | wetland disabled | proxy | none (sensitivity) |

The toolik06 input grid is 10x10 (28 well-drained and 72 poorly drained cells, all CMT06). The run-mask enables one cell.

## 2. Requirements and options
**What T3 must add beyond T3-lite:**
1. Mass conservation between species: d(Hg(II))/dt includes -M + D_ox, d(MeHg)/dt = M - D, and d(Hg(0)) includes D_red.
2. Hg(II) processes see Hg(II), not THg.
3. MeHg partitioning and export that are its own (Kd, DOC link) rather than "same share as Hg(II)".
4. A rate law with explicit, literature-bounded factors for temperature, redox/saturation, electron donor (labile or
   anaerobic C), sulfate and trophic/fen class. Each is switchable, and the ones that are unconstrained default to 1.
5. A fen/wetland class that the CMT structure can carry, without site-specific tuning.
6. B2F 9/9 and MeHg 7/7 must hold with the recommended set.

### A) Index-plus (no species split)
Keep the index and add only the new rate factors and a MeHg Kd for mobility. This is cheap and changes no Hg(II) numbers,
but it does not meet requirement 1-2, so it is still T3-lite. Useful only as the "closure 0" sub-mode of B.

### B) MeHg as a mass-conserving sub-pool of the layer Hg (recommended)
- **State:** unchanged. s.soil + s.diss stay THg per layer; s.mehg is the MeHg part of it (0 <= MeHg <= THg, as now). Hg(II) =
  THg - MeHg. Dissolved and sorbed MeHg come from equilibrium partitioning each step. No restart change.
- **Methylation / demethylation:** the existing exponential-relaxation update stays. With hg_t3_closure 1 the month's net
  change dM = MeHg_new - MeHg_old is booked as a species transfer inside THg (THg unchanged), and the gross terms are
  diagnosed: M = P (gross), D = M - dM. With hg_t3_fred > 0 the reductive share f_red x D leaves THg as Hg(0) to the soil
  Hg(0) pathway (the B2F escape probability P_l applies). The default is f_red = 0, so THg stays exactly conserved.
- **Closure (hg_t3_closure 1):** these processes use Hg(II) = THg - MeHg:
  - decay-release re-emission pot and its escape probability;
  - Kd sorption / dissolved share for Hg(II) (MeHg gets its own share, below);
  - the DOC-limited Hg(II) export (already done at l.1220);
  - thaw release;
  - deep routing / retention of Hg(II) (MeHg retention follows its own mobile share).
  Fire keeps today's rule: MeHg is combusted with the layer Hg, which is conservative.
- **Partitioning:** phi_MeHg,l = liq / (liq + Kd_MeHg M_l), with log10 Kd_MeHg = log10 Kd_Hg(II),layer + hg_t3_dlogkd_mehg.
  - Daily transfers carry F_MeHg = F_w x C_MeHg,diss, where F_w is the water flux, and Hg(II) transfers carry F_w x
    C_Hg(II),diss.
  - THg export becomes the sum of the two species.
  - With dlogkd 0 and closure 0 this reduces exactly to today's mobfac 1 rule.
  - DOC link: the DOC step-1 scaled-Kd coupling applies to MeHg in the same form. The Hg:DOC cap stays Hg(II)-only.
- **Rate law (per layer, month; daily factors where the code already has them):**
  P_l = k_m,driver x S_l x A_l x Q10_m^((T - t_ref)/10) x f_frz,l x g_class x f_SO4 x f_avail,l
  D_l = k_d x Q10_d^((T - t_ref)/10) x MeHg_l
  - S_l = Hg(II)_l with closure, else THg_l (today).
  - A_l = today's driver: CH4 production / C, or sat_d x RH / C.
  - f_frz = today's daily freeze factor.
  - g_class = trophic/fen factor from the CMT class map (default 1).
  - f_SO4 = optional sulfate factor (default 1, §4.4).
  - f_avail = (1 - w_a) + w_a phi_Hg(II),l / phi_ref, a bioavailability blend. The default w_a = 0 means total Hg(II).
- **Cost of mismatch:** none for off; on, Hg outputs shift by about the MeHg share of THg in the layers where each process
  acts (§3.2).

### C) Full separate Hg(II) / MeHg species arrays
Duplicate every pool and every transport, sorption, fire and geometry routine for MeHg, as a new restart field set (v5).
This is the cleanest, but the same physics as B at monthly equilibrium partitioning. It costs about 2-3x the development
time, a restart-format change and its regression suite. It is only worth doing if MeHg ever needs kinetic sorption (a
non-equilibrium MeHg sorbed pool). Deferred.

## 3. Read-only estimates (B2F outputs, last 20 yr; 0 core-h)
### 3.1 How fast the index cycles
| case | active-layer Hg (mg/m2) | MeHg stock (ug/m2, JJAS %) | gross MEHGPROD (ug/m2/yr) | turnover (1/yr) | total evasion (ug/m2/yr) | Hg(0) from f_red = 0.1 % (ug/m2/yr) |
|---|---|---|---|---|---|---|
| TVC | 12.37 | 93 | 1838 | 20 | 11.2 | 1.84 (+16 % of evasion) |
| HPC | 16.01 | 15 | 423 | 28 | 7.0 | 0.42 (+6 %) |
| Toolik fen (CH4) | 28.76 | 639 | 4961 | 7.8 | 3.3 | 4.96 (+150 %) |
| Toolik proxy | 28.06 | 33 | 284 | 8.6 | 3.2 | 0.28 |

- Gross production is about 40-720 times the deposition (6.7-7.9 ug/m2/yr) [derived]. This is a fast internal cycle that %MeHg constrains
  only through the ratio k_m,eff / k_d. Gross rates are not observed at our sites.
- **Consequence:** any demethylation product other than Hg(II) is a large new Hg(0) source.
  - Marvin-DiPasquale et al. 2000 found oxidative demethylation (product CO2, so Hg(II)) dominant at low-contamination sites
    with k_deg <= 0.2 /d. The reductive mer pathway (CH4 + Hg(0)) dominated only at mine-impacted sites.
  - Oremland et al. 1991 found anaerobic demethylation products to be mainly CO2.
  - So f_red = 0 is both the literature default and the only value compatible with the B2F evasion fit. A non-zero f_red
    would need a joint refit of fd.
### 3.2 Closure effect on Hg(II)
- MeHg / THg in the active layer (JJAS) is 0.75 % (TVC), 0.09 % (HPC) and 2.2 % (Toolik fen).
- The closure lowers the Hg(II) seen by sorption, evasion and export by about that share in the producing layers. The
  surface 10-20 cm, where most re-emission happens, carries up to about 1.9 % at TVC (JJAS profile).
- **Expected:** TVC Ro and evasion change by <= 1-2 %. The TVC Ro margin is 0.008 of 0.41 (about 2 %), so the closure may
  need a small fd refit (one-box, as B2F did). HPC is essentially unaffected. Toolik is report-only.
### 3.3 MeHg Kd and export
- Literature soil/soil-water log Kd (Allison & Allison 2005, EPA/600/R-05/074, Table 1/3):
  - Hg(II): median 3.8 (2.2-5.8, n = 17);
  - MeHg: median 2.8 (mean 2.7, 1.3-4.8, n = 11).
  - Lyon et al. 1997 give soil-water means of 3.3e3-6.0e4 L/kg for Hg(II) and 2.0e1-6.7e3 for MeHg.
  - So the median offset is about -1 log unit, with a very wide spread.
- Organic soils are different: MeHg binds reduced-S (thiol) groups in soil and stream organic matter (Qian et al. 2002;
  Skyllberg 2008). So in peat and organic horizons, the offset to the model's Hg(II) organic Kd (log 3.356) may be much
  smaller than the mineral-soil median.
- **Export sensitivity [derived]:** where Kd M >> liq, the dissolved share scales as 1/Kd. An offset of -1 raises the mobile
  MeHg about 10x. MeHg export would scale less than 10x where it is limited by production or the flow path, but HPC (36.3 of
  the 50 band) and the fen (40.6) have little headroom. **The export band, not the literature median, has to pick
  dlogkd.** The scan in §6 is {0, -0.5, -1.0}.
- MeHg:THg export in B2F is 2.7 % (TVC), 3.7 % (HPC) and 5.2 % (fen), from May-Oct MeHg against annual THg dissolved
  export. That is within or at the top of the 1-5 % order of the Mackenzie and Scotty Creek observations
  (lit_7questions Q6). Raising MeHg mobility would push the fen above it.

## 4. Process laws and parameter ranges (all citations checked, §9)
### 4.1 Temperature
- Keep Q10_m = 1 on top of the C-turnover driver, which already carries a temperature response. Keep Q10_d = 2.
- Literature: k_m rose about 5x from 5 to 20 C (Q10 about 3, derived) and k_d showed no significant temperature response in
  boreal wetland soils (Sagisaka Méndez & Mitchell 2026). k_d varied only 2.4x along a chronosequence where k_m varied 26-53x
  (Hu et al. 2020).
- **Known conflict:** in M1, Q10_d = 1 could not keep the 7:1 TVC/HPC contrast with one k (ratio 14.6-16.8). So Q10_d 2 is
  a structural fit, not a literature value. V22/V32 (Q10_m 2-3, Q10_d 2) are kept as sensitivities, and T3 does not reopen
  this.
### 4.2 Redox / water table / saturation
- Keep today's representation:
  - the daily saturation factor (satd) x the freeze cutoff for the proxy;
  - CH4 production, which implies anaerobic conditions, for wetland CMTs.
- Rationale: enhanced methylation follows the change to a saturated, collapsed or minerotrophic state, not thaw itself, and
  "redox conditions alone do not drive" it (Tarbier et al. 2021).
- No new redox variable: dvmdostem has no Eh, Fe or sulfide state.
### 4.3 Electron donor / labile C
- Keep A = anaerobic C turnover. Methylation correlated with CO2 + CH4 production (Sagisaka Méndez & Mitchell 2026) and with
  CH4 and Fe(II) production, and labile C fuelled early MeHg production (Yang et al. 2016).
- In peat mesocosms, organic C alone did **not** raise MeHg production. Sulfate alone raised it to about 200 pg/L/d, and
  sulfate + C to about 500 pg/L/d (Mitchell, Branfireun & Kolka 2008). So a DOC/labile-C multiplier without sulfate is not
  supported. **No separate DOC factor.**
### 4.4 Sulfate (optional, default off)
- **Evidence:**
  - Sulfate stimulates methylation in freshwater sediments (Gilmour, Henry & Mitchell 1992) and in a boreal poor fen
    (Branfireun et al. 1999).
  - In an experimental wetland, about 4x the annual sulfate load raised porewater MeHg from about 0.5 to 1.63 ng/L (3x)
    two weeks after the first addition, and the MeHg flux from the wetland about 2.4x. Later additions did not raise
    porewater MeHg further (Jeremiason et al. 2006).
  - At high sulfate, sulfide inhibits methylation: neutral HgS0 falls as sulfide rises (Benoit et al. 1999). The optimum is
    0.2-0.5 mM sulfate (Gilmour & Henry 1991, as summarised by Ullrich et al. 2001).
  - Methylation became sulfate-limited at 25 C (Sagisaka Méndez & Mitchell 2026).
- **Proposed form:** f_SO4 = [S/(K_S + S)] x [K_I/(K_I + S)] / f_ref, normalised to 1 at a reference S_ref. S is a
  **class or site constant** (porewater sulfate, mM), because there is no S cycle and no sulfate input. K_S and K_I are
  chosen so the optimum falls at 0.2-0.5 mM.
- **Default hg_t3_so4 0 (f_SO4 = 1).** Fitting is not possible: no porewater sulfate data at TVC, HPC or Toolik were found.
  The switch exists so scenario work (for example, minerotrophic fens fed by groundwater sulfate) can use it.
### 4.5 Substrate availability (optional, default w_a = 0)
- The porewater Hg(II) concentration is decisive (Liem-Nguyen et al. 2021), and k_m differs > 100x between Hg(II) forms
  (Jonsson et al. 2012). But in M1 the fully dissolved substrate (subst 1 or 2) failed by an order of magnitude (ratio
  0.8-1.0 or 91, k x740).
- So f_avail is a blend with weight w_a: it defaults to 0, and only small w_a (<= 0.2) is a sensitivity.
### 4.6 Demethylation
- k_d = 3 /month at t_ref (about 0.1 /d) stays. It is within <= 0.2 /d at low-Hg sites (Marvin-DiPasquale et al. 2000) and
  <= 0.5 /d without an Hg source (Helmrich et al. 2022). The production-weighted k_d,eff is 0.05 /d.
- Products: Hg(II) (f_red = 0, §3.1).
- Photodemethylation (3.69e-3 m2/E, Lehnherr & St. Louis 2009) belongs to open water. dvmdostem has no stream or pond, so it
  is out of scope (report-only note).
### 4.7 Fen / trophic class factor (optional, default 1)
- Tracer k_m rises plateau < bog < fen: 8e-4, 4e-3, 0.059 and 0.146 /d (Thompson et al. 2025). The hgcA transcripts were
  highest in fens and absent in palsa / shallow bog (McDaniel et al. 2020).
- The model already makes the contrast through its drivers: the TVC/HPC ratio is 8.1, and the fen comes from the CH4
  driver. So g_class stays 1 for every class unless a calibration shows a driver cannot reach a class target. If it is
  needed, it should be **one** fen factor, bounded by the tracer ratios (poor fen / bog = 15x), and not a per-site factor.

## 5. Wetland / fen class with the CMT structure
**What exists:**
- one CMT per cell (veg_class input);
- the wetland CMT list (wetland.cmts) and the CH4 module (layer CH4 production, water table) for those CMTs;
- drainage_class 0/1 (input);
- daily layer T, liquid water, ice, water table and RH;
- the DOC step-1 diagnostic DOC;
- CMT31 Boreal Bog and CMT32 Boreal Fen parameterisations in parameters/ (not used at our sites).

**What is missing:**
- a sub-grid wetland/fen fraction per cell;
- trophic status or minerotrophy and groundwater chemistry (wetland.groundwater is false);
- sulfate, Fe, sulfide and pH;
- MeHg deposition;
- site soil %MeHg or stream MeHg at TVC, HPC or Toolik (none found; the targets are Thompson's NW-Canada classes).

**Proposal (W1 + W3):**
- **W1 class map** (model_settings `hg_t3_class`: {"fen": [6, 32, ...], "bog": [31, 51?], "upland": [...]}, default absent
  = no class, g_class 1). The class picks g_class and an optional class sulfate S. The driver is still chosen by ch4_on as
  now.
  - Mapping TVC's CMT51 to "bog" only reproduces today's assumption. It is the TVC target problem flagged in the T2 design,
    and T3 does not solve it without TVC %MeHg data.
  - Recommendation: keep TVC on the bog target but label it "assumed class" in reports.
- **W2 drainage proxy** (rejected as a class): TVC and HPC are both drainage_class 1, so it does not separate them. It is
  still useful as a **within-Toolik test**.
  - A copied run-mask enabling 6 well-drained (0) and 6 poorly drained (1) CMT06 cells gives a fen-vs-upland contrast
    under one climate, with the CH4 driver.
  - Report-only plausibility: wetland vs upland MeHg yields 26-79x (1.84-5.55 vs 0.07 mg/ha/yr, St. Louis et al. 1994),
    fen-drainage vs bog water MeHg 4.5-14.5x (Gordon et al. 2016), collapse fen vs plateau %MeHg 13x (Tarbier et al. 2021).
- **W3 landscape mixing in post-processing:** landscape MeHg export = sum over classes of f_class x export_class from
  separate class runs. There is no code change, and f_class comes from a map when one is available (missing for TVC/HPC).
  This is the honest way to compare with catchment yields (Thompson et al. 2023: 5.2-36 ng/m2 per open-water season in
  peatland catchments).
- **W4 sub-grid fen tile inside a cohort** (rejected for T3): it would need a second soil column per cell.

## 6. Switches, tests, phases
### 6.1 Switches (all default off; need hg_mehg 1, hg_t2 1, hg_t2_diss 1, hg_transport 1, hg_t2_mehg_transport 1)
| switch | values (default) | meaning |
|---|---|---|
| hg_t3 | 0/1 (0) | master; with all sub-switches at default it must equal B2F bitwise |
| hg_t3_closure | 0/1 (0) | Hg(II) processes use THg - MeHg; species transfer booked |
| hg_t3_dlogkd_mehg | -3..+1 (0) | log10 Kd_MeHg - log10 Kd_Hg(II) of the layer; 0 = today's mobfac 1 |
| hg_t3_fred | 0..1 (0) | reductive demethylation share to Hg(0) (needs closure 1) |
| hg_t3_wavail | 0..1 (0) | bioavailable-substrate weight w_a |
| hg_t3_so4 | 0/1 (0), with hg_t3_so4_ks, _ki, _sref and class S | sulfate factor (§4.4) |
| hg_t3_class | object (absent) | CMT class map, g_class per class (1) |

New outputs (spec rows invalid by default):
- MEHG_DEMETH: gross demethylation by layer and total;
- MEHG_DISS: dissolved MeHg (ug/m2) and MEHG_PWC: porewater MeHg (ng/L) by layer;
- HG_MEHGEXP: the real removal (already exists);
- HG_T3RED: Hg(0) from reductive demethylation;
- HG_T3BAL: species balance error.

The restart stays v4 and parameters/ are untouched. All values go through model_settings, with a new overlay config
(hg_t2_overlay_T3x.json) written next to B2F.

### 6.2 Tests
- **Unit (test_hg_tracer.cpp):**
  - parse, ranges and dependencies;
  - species conservation (THg constant under methylation and demethylation with f_red 0; THg loss = f_red D otherwise);
  - closure algebra;
  - partition algebra (dlogkd 0 gives today's mobile share exactly);
  - export of two species sums to THg export;
  - f_SO4 optimum and inhibition, and normalisation to 1 at S_ref;
  - g_class 1 identity;
  - frozen layer: no methylation;
  - restart round trip.
- **Model identities (8 run sets: TVC / HPC spin + main, Toolik proxy / CH4 sp + tr):**
  - hg_t3 absent vs B2F: 0 differing variables.
  - hg_t3 1 with defaults vs B2F: 0 differing variables, or <= 1e-12 relative if the refactor changes operation order.
    This is the index-equivalence test.
  - Hg-off production 8-case PASS.
  - closure 0 + any dlogkd / factors: 0 non-MeHg differences (the M-phase no-feedback test still holds).
  - C / N / water bitwise equal in every T3 run (Hg has no feedback on C).
- **Acceptance for a recommended T3 set:**
  - 9/9 scored bands and MeHg 7/7.
  - HG_BALERR and HG_T3BAL <= 1e-8 ug/m2.
  - k_d,eff <= 0.5 /d, and k_m,eff below the class tracer k_m.
  - MeHg:THg export reported against 1-5 %.
  - If a band fails, report it and do not force the fit (as in M2 D3).
### 6.3 Phases (test runs n <= 75, <= 4 cores; builds -j1/-j2; transient-only from copied restarts where possible)
Unit costs from LEDGER.txt: a full 4-case set (TVC/HPC spin + main, Toolik proxy + CH4 sp + tr) is about 0.05 core-h (B2F
0.048). One Toolik CH4 cell is 58 s, about 0.016 core-h. A 75-yr W6 segment is about 0.007.

| phase | content | runs | core-h |
|---|---|---|---|
| T3-P0 | code + unit tests; identities (absent, index-equivalence, Hg-off prod) | 2 sets + prod | 0.17 |
| T3-P1 | closure 1 (1 set); dlogkd {-0.5, -1.0} with closure (2 sets); linear k rescale + fd one-box refit if Ro moves; 1 verification set | 4-5 sets | 0.25 |
| T3-P2 | sensitivities: f_red 1e-3 (1 set; predicted analytically first, §3.1); fen g_class bounded or SO4 x3 at Toolik CH4 only (2 Toolik cases) | about 1.5 sets | 0.08 |
| T3-P3 | Toolik drainage contrast: copied run-mask, 6 + 6 CMT06 cells, CH4 driver, sp + tr, recommended set | 12 cells | 0.20 |
| T3-P4 | warming: TVC / HPC ctl + w8 from copied 2024 restarts; Toolik fen w8 | 5 segments | 0.04 |
| docs | WORK_NOTES §24, figures, overlay, ledger line "PHASE T3" | 0 | 0 |
| **total** | | | **about 0.74 expected; cap 1.0** |

Stop rules:
- If index-equivalence fails, stop and fix before P1.
- If closure 1 moves a scored band out, the fd refit is the only allowed fix. If that fails, report it and leave closure as
  a sensitivity.
- If no dlogkd <= 0 passes the export bands, keep dlogkd 0 and report the conflict with the literature Kd.

## 7. Design choices for Elchin (recommendation in bold)
- **T1 Pool structure:** **B, MeHg as a mass-conserving sub-pool of the layer Hg (no restart change)**; or C, separate
  species arrays (v5 restart, 2-3x the work); or A, index-plus (not real T3).
- **T2 Closure:** **build it and adopt it in the T3 set if 9/9 holds (a one-box fd refit allowed).** Keep closure 0 as the
  no-feedback mode for MeHg-only studies.
- **T3 Rate-law core:** **keep the calibrated M4V1 driver (A x sat x freeze, THg or Hg(II) substrate) and add only factors
  that default to 1.** Alternative: re-derive a new rate law and recalibrate from scratch (M1 showed the literature-pure
  choices break the class contrast).
- **T4 Temperature:** **keep Q10_m 1 (on top of the C driver) and Q10_d 2, with the literature conflict documented**;
  V22/V32 stay as sensitivities.
- **T5 Demethylation products:** **100 % oxidative, to Hg(II) (f_red 0)**, with f_red only as a sensitivity because 0.1 % of
  gross production is already 16 % of TVC evasion.
- **T6 Sulfate:** **switch built, default off, not fitted** (no sulfate data at the sites); or omit it entirely.
- **T7 MeHg partitioning:** **own Kd as an offset from Hg(II); scan {0, -0.5, -1.0} and let the export bands pick; MeHg
  stays outside the Hg:DOC cap and uses the DOC step-1 scaled-Kd coupling.** Alternative: impose the Allison median offset
  (-1.0) and accept possible export-band failures.
- **T8 Wetland / fen:** **W1 CMT class map (g_class 1 unless needed; at most one bounded fen factor) + W3 landscape mixing
  in post-processing + the Toolik drainage-contrast test as a report-only plausibility check.** TVC stays on the assumed
  bog class. W4 (sub-grid fen tile) is deferred.
- **T9 Calibration:** **fit only k_m,proxy and k_m,CH4 (linear rescale), and pick dlogkd from the scan.** Fix everything
  else from the literature or at 1. New report-only checks: fen/upland yield ratio, MeHg:THg export, porewater MeHg, k_m,eff,
  k_d,eff.
- **T10 Out of scope:** **photodemethylation (no surface water), MeHg deposition (0), kinetic MeHg sorption, and turning the
  CH4 module on at TVC/HPC** (that is a C-configuration change and needs its own approval).
- **T11 Identity bar:** **off = B2F bitwise; hg_t3 1 at defaults = B2F bitwise (index-equivalence)**, and closure-on runs
  are not identical by design.
- **T12 Budget and inputs:** **0.74 core-h expected, cap 1.0, own ledger line "PHASE T3"; a copied 12-cell Toolik
  run-mask in our own tree for P3** (cal-cmt06 stays read-only and no input data change).

## 8. Calibration targets, data and identifiability
**Targets kept (primary, scored as now):**
- TVC 0.7 % ±20 %, HPC 0.1 % ±20 %, ratio 4.67-10.5, Toolik fen 2.0 % ±20 % (Thompson et al. 2025 classes; Jun-Sep active
  layer);
- May-Oct MeHg export 2-50 ng/m2 at HPC, TVC and the Toolik fen.

**Secondary / report-only (new in T3):**
- 0-50 cm order plateau < bog < fen (Tarbier et al. 2021: top-50 cm 0.48 / 1.7 / 6.1 %);
- MeHg:THg export 1-5 % (Emmerton et al. 2013; Thompson et al. 2023);
- dissolved share of exported MeHg, about half (Emmerton et al. 2013: 44-51 %);
- porewater MeHg of order 0.5 ng/L in a boreal peatland at baseline (Jeremiason et al. 2006; not an Arctic site);
- the Toolik drainage contrast against the wetland/upland yield ratio of 26-79x (St. Louis et al. 1994) and the 4.5-14.5x
  fen/bog water MeHg (Gordon et al. 2016);
- fen export against wetland yields of 184-555 ng/m2/yr (St. Louis et al. 1994, boreal Ontario). Our fen band ceiling of
  50 ng/m2 per season is below that, which is worth noting rather than fixing;
- k_m,eff below the tracer k_m, and k_d,eff <= 0.2-0.5 /d.

**Identifiability:**
- At steady state %MeHg fixes only k_m,eff / k_d. k_d and Q10_d are separated only by the seasonal cycle and the TVC/HPC
  contrast, and gross rates are unobserved. So any process that depends on gross fluxes (f_red, MeHg-specific evasion) is
  unconstrained.
- k_m and g_class are collinear: one site per class means a class factor is a site factor. Hence g_class 1.
- f_SO4 and f_avail are collinear with k_m and have no data. Hence default off.
- dlogkd is constrained only by the export bands, together with production. Two bands (HPC and the fen) close to their
  ceilings make it effectively one-sided (<= 0 offsets cost export headroom).
- Three sites against 2-3 fitted parameters is the maximum that is defensible.

**Data gaps:**
- no soil %MeHg, porewater MeHg or stream MeHg at TVC, HPC or Toolik;
- no porewater sulfate or groundwater chemistry;
- no fen fraction for the TVC / HPC landscapes;
- the TVC class mapping is an assumption;
- no Arctic tundra k_m / k_d with simultaneous temperature response;
- the organic-soil MeHg Kd is not separated from mineral soils in the compilations.

## 9. References (DOIs resolved at doi.org or checked against Crossref / publisher / PubMed, 24-26 Nov 2026)
- Allison, J. D. & Allison, T. L. (2005) Partition coefficients for metals in surface water, soil, and waste. US EPA,
  EPA/600/R-05/074. https://semspub.epa.gov/work/01/466653.pdf (soil log Kd Hg(II) 3.8 [2.2-5.8], MeHg 2.8 [1.3-4.8])
- Benoit, J. M., Gilmour, C. C., Mason, R. P. & Heyes, A. (1999) Sulfide controls on mercury speciation and bioavailability to
  methylating bacteria in sediment pore waters. Environ. Sci. Technol. 33, 951-957. https://doi.org/10.1021/es9808200
- Branfireun, B. A. et al. (1999) In situ sulphate stimulation of mercury methylation in a boreal peatland. Global
  Biogeochem. Cycles 13, 743-750. https://doi.org/10.1029/1999GB900033
- Emmerton, C. A. et al. (2013) Mercury export to the Arctic Ocean from the Mackenzie River, Canada. Environ. Sci. Technol.
  47, 7644-7654. https://doi.org/10.1021/es400715r
- Gilmour, C. C. & Henry, E. A. (1991) Mercury methylation in aquatic systems affected by acid deposition. Environ. Pollut.
  71, 131-169. https://doi.org/10.1016/0269-7491(91)90031-Q (sulfate optimum as summarised by Ullrich et al. 2001)
- Gilmour, C. C., Henry, E. A. & Mitchell, R. (1992) Sulfate stimulation of mercury methylation in freshwater sediments.
  Environ. Sci. Technol. 26, 2281-2287. https://doi.org/10.1021/es00035a029
- Gordon, J. et al. (2016) Mercury and methylmercury biogeochemistry in a thawing permafrost wetland complex, Northwest
  Territories, Canada. Hydrol. Process. 30, 3627-3638. https://doi.org/10.1002/hyp.10911 (fen vs bog drainage MeHg 4.5-14.5x)
- Helmrich, S. et al. (2022, online 2021) Critical review of mercury methylation and methylmercury demethylation rate
  constants in aquatic sediments for biogeochemical modeling. Crit. Rev. Environ. Sci. Technol. 52, 4353-4378.
  https://doi.org/10.1080/10643389.2021.2013073
- Hu, H. et al. (2020) Shifts in mercury methylation across a peatland chronosequence: from sulfate reduction to
  methanogenesis and syntrophy. J. Hazard. Mater. 387, 121967. https://doi.org/10.1016/j.jhazmat.2019.121967
- Jeremiason, J. D. et al. (2006) Sulfate addition increases methylmercury production in an experimental wetland. Environ.
  Sci. Technol. 40, 3800-3806. https://doi.org/10.1021/es0524144 (PMID 16830545)
- Jonsson, S. et al. (2012) Mercury methylation rates for geochemically relevant Hg(II) species in sediments. Environ. Sci.
  Technol. 46, 11653-11659. https://doi.org/10.1021/es3015327
- Lehnherr, I. & St. Louis, V. L. (2009) Importance of ultraviolet radiation in the photodemethylation of methylmercury in
  freshwater ecosystems. Environ. Sci. Technol. 43, 5692-5698. https://doi.org/10.1021/es9002923
- Liem-Nguyen, V. et al. (2021) Methylmercury formation in boreal wetlands in relation to chemical speciation of mercury(II)
  and concentration of low molecular mass thiols. Sci. Total Environ. 755, 142666. https://doi.org/10.1016/j.scitotenv.2020.142666
- Lyon, B. F. et al. (1997) Calculation of soil-water and benthic sediment partition coefficients for mercury. Chemosphere 35,
  791-808. https://doi.org/10.1016/S0045-6535(97)00200-2
- Marvin-DiPasquale, M. et al. (2000) Methyl-mercury degradation pathways: a comparison among three mercury-impacted
  ecosystems. Environ. Sci. Technol. 34, 4908-4916. https://doi.org/10.1021/es0013125
- McDaniel, E. A. et al. (2020) Expanded phylogenetic diversity and metabolic flexibility of mercury-methylating
  microorganisms. mSystems 5, e00299-20. https://doi.org/10.1128/msystems.00299-20
- Mitchell, C. P. J., Branfireun, B. A. & Kolka, R. K. (2008) Assessing sulfate and carbon controls on net methylmercury
  production in peatlands: an in situ mesocosm approach. Appl. Geochem. 23, 503-518.
  https://doi.org/10.1016/j.apgeochem.2007.12.020
- Oremland, R. S., Culbertson, C. W. & Winfrey, M. R. (1991) Methylmercury decomposition in sediments and bacterial
  cultures. Appl. Environ. Microbiol. 57, 130-137. https://doi.org/10.1128/aem.57.1.130-137.1991
- Qian, J. et al. (2002) Bonding of methyl mercury to reduced sulfur groups in soil and stream organic matter. Geochim.
  Cosmochim. Acta 66, 3873-3885. https://doi.org/10.1016/S0016-7037(02)00974-2
- Sagisaka Méndez, S. & Mitchell, C. P. J. (2026) Effects of temperature on mercury methylation and demethylation in boreal
  wetland soils. J. Geophys. Res. Biogeosci. 131, e2025JG009302. https://doi.org/10.1029/2025JG009302
- Skyllberg, U. (2008) Competition among thiols and inorganic sulfides and polysulfides for Hg and MeHg in wetland soils and
  sediments under suboxic conditions. J. Geophys. Res. Biogeosci. 113, G00C03. https://doi.org/10.1029/2008JG000745
- St. Louis, V. L. et al. (1994) Importance of wetlands as sources of methyl mercury to boreal forest ecosystems. Can. J.
  Fish. Aquat. Sci. 51, 1065-1076. https://doi.org/10.1139/f94-106
- Tarbier, B. et al. (2021) Permafrost thaw increases methylmercury formation in subarctic Fennoscandia. Environ. Sci.
  Technol. 55, 6710-6717. https://doi.org/10.1021/acs.est.0c04108
- Thompson, L. M. et al. (2023) Concentrations and yields of mercury, methylmercury, and dissolved organic carbon from
  contrasting catchments in the discontinuous permafrost region. Water Resour. Res. 59, e2023WR034848.
  https://doi.org/10.1029/2023WR034848
- Thompson, L. M. et al. (2025) Production of methylmercury in peatlands following permafrost thaw increases along a trophic
  gradient. Environ. Sci. Technol. 59, 19457-19467. https://doi.org/10.1021/acs.est.5c04510
- Ullrich, S. M., Tanton, T. W. & Abdrashitova, S. A. (2001) Mercury in the aquatic environment: a review of factors
  affecting methylation. Crit. Rev. Environ. Sci. Technol. 31, 241-293. https://doi.org/10.1080/20016491089226
- Yang, Z. et al. (2016) Warming increases methylmercury production in an Arctic soil. Environ. Pollut. 214, 504-509.
  https://doi.org/10.1016/j.envpol.2016.04.069

Numbers marked [derived] are my arithmetic on the quoted values. The sources marked earlier in lit_7questions.md (24 Nov)
were reused with their checks. Login-walled full texts were not read; values come from abstracts, the EPA report or public
pages.
