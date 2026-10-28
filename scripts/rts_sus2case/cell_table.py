#!/usr/bin/env python3
"""Per-cell RTS susceptibility / inventory table for the Gydan 10x10 grid (phase 1).

One row per land cell of spatial10x10/inputs/grid.json (90 cells, 0.5 deg,
centres 67.75..72.25N, 75.75..80.25E). The 0.5-deg blocks of the aggregated
Jia file line up with the grid (block edges offset ~0.0016 deg, negligible), so
cells are matched by centre coordinates (|d| < 0.01 deg), no regridding.

Columns (all fractions relative to the TOTAL 0.5-deg cell area, as in aggregate.py):
  jia_{mean,gt050,gt075}_ssp5_{2030,2050,2070,2090}   Jia et al. 2026 SSP5-8.5
      cell-mean susceptibility (0-1) and area fractions with s>0.5 / s>0.75
  jia_{mean,gt050,gt075}_current   Jia Current map, REFERENCE ONLY: Current is
      built from a different predictor set and must never be differenced
      against an SSP period (grid-mean s 0.76 vs ~0.5 for the SSPs).
  valid_fraction_ssp5_YYYY, valid_fraction_current   Russia + valid pixels / cell
  russia_land_fraction, cell_area_km2   from the aggregation file
  nesterova_count   Nesterova et al. 2025 West Siberia RTS v0.1 points
                    (doi:10.1594/PANGAEA.974406; ws_count.npy from fit_ws.py)
  jia_inv_count     Jia et al. 2026 training inventory (RTS dataset.xlsx), 100 m
                    de-duplicated (russia_rts_inventory_count_0p5deg.nc, calib.py)
  arts_count, arts_area_ha   ARTS v1.0 positive polygons (Yang et al. 2025,
                    doi:10.18739/A2PK0738B), by representative point
  darts_fobs        DARTS v1.3 footprint fraction of valid area (darts_fit.py),
                    NaN where DARTS covers <50% of the cell ("where available")
  darts_n_footprints  DARTS footprint-union polygons by representative point
                    (NaN where darts_fobs is NaN)
  jia_nan_ssp5 / jia_nan_current   1 if any of mean/gt050/gt075 is NaN
  case_a_F_{mean,gt050}_YYYY, case_b_* : provisional derived columns, see
                    derive_cases(); recomputed in phase 2.
NaN policy: values are written as NaN (never silently 0); flags list them.

Usage (VM, env with geopandas + xarray):
  cell_table.py --grid .../spatial10x10/inputs/grid.json --sus-dir ~/box_archive/rts-susceptibility -o table.csv
"""
import argparse, json, os
import numpy as np

PERIODS = (2030, 2050, 2070, 2090)
ANCHOR_HA = 1.44      # ARTS v1.0 Russia mean RTS area (union footprints), Gydan anchor
MF_CENTRAL, MF_LOW, MF_HIGH = 7.5e-4, 1.3e-4, 2.0e-3   # rts_postproc JIA2026_CALIBRATED_MAX_FRACTION (+range)


def match_index(axis, values, tol=0.01):
    """Index of each value in a coordinate axis (nearest; error if farther than tol)."""
    axis = np.asarray(axis, float); idx = []
    for v in np.atleast_1d(values):
        k = int(np.argmin(np.abs(axis - v)))
        if abs(axis[k] - v) > tol:
            raise ValueError(f'{v} not on axis (nearest {axis[k]})')
        idx.append(k)
    return np.array(idx)


def point_cell(lat, lon, top_edge, west_edge=-180.0, d=0.5):
    """(row, col) of points on the aggregation grid (rows from the top edge southwards)."""
    return (np.floor((top_edge - np.asarray(lat)) / d).astype(int),
            np.floor((np.asarray(lon) - west_edge) / d).astype(int))


def linear_F(s, mf, thr=0.5):
    """rts_postproc susceptibility_to_target_fraction: mf*(s-thr)/(1-thr) if s>thr else 0; NaN stays NaN."""
    s = np.asarray(s, float)
    return np.where(np.isnan(s), np.nan, np.where(s > thr, mf * (s - thr) / (1 - thr), 0.0))


def smooth3x3(field, land):
    """3x3 nan-mean over land cells of a (ny, nx) field (edges use the available neighbours)."""
    f = np.where(land & np.isfinite(field), field, np.nan); ny, nx = f.shape; out = np.full_like(f, np.nan)
    for j in range(ny):
        for i in range(nx):
            if not land[j, i]: continue
            w = f[max(j - 1, 0):j + 2, max(i - 1, 0):i + 2]
            if np.isfinite(w).any(): out[j, i] = np.nanmean(w)
    return out


def derive_cases(rows, land_ji, shape=(10, 10)):
    """Provisional case A / case B columns (approved design, phase 2 recomputes them).

    case A (Jia susceptibility, SSP5 only): F_k = mf*(x-0.5)/0.5 per Jia period cohort,
      x = cell-mean s (as coded in rts_postproc jia2026.py --calibrated) -> case_a_F_mean_k,
      x = area_fraction_gt050 (the wording of the approved design) -> case_a_F_gt050_k;
      final disturbed fraction = sum_k F_k (cap 1 never binds at mf<=2e-3).
    case B (Nesterova inventory): A0 = count*1.44 ha/cell_area; r_k = gt050_k/gt050_2030
      (SSP5), smoothed 3x3 over land cells -> case_b_A_k = A0*r_k_smooth.
    """
    ny, nx = shape; land = np.zeros(shape, bool)
    for j, i in land_ji: land[j, i] = True
    pos = {(r['j'], r['i']): r for r in rows}
    for r in rows:
        for p in PERIODS:
            r[f'case_a_F_mean_{p}'] = float(linear_F(r[f'jia_mean_ssp5_{p}'], MF_CENTRAL))
            r[f'case_a_F_gt050_{p}'] = float(linear_F(r[f'jia_gt050_ssp5_{p}'], MF_CENTRAL))
        r['case_b_A0'] = r['nesterova_count'] * ANCHOR_HA * 0.01 / r['cell_area_km2']
    for p in PERIODS:
        ratio = np.full(shape, np.nan)
        for (j, i), r in pos.items():
            g0, g = r['jia_gt050_ssp5_2030'], r[f'jia_gt050_ssp5_{p}']
            if np.isfinite(g0) and np.isfinite(g) and g0 > 0: ratio[j, i] = g / g0
        sm = smooth3x3(ratio, land)
        for (j, i), r in pos.items():
            r[f'case_b_ratio_raw_{p}'] = float(ratio[j, i]); r[f'case_b_ratio_s3_{p}'] = float(sm[j, i])
            r[f'case_b_A_{p}'] = float(r['case_b_A0'] * sm[j, i]) if np.isfinite(sm[j, i]) else (0.0 if r['case_b_A0'] == 0 else float('nan'))
    return rows


def extract(grid_json, sus_dir):
    import xarray as xr
    g = json.load(open(grid_json)); cells = g['cells']
    ds = xr.open_dataset(os.path.join(sus_dir, 'russia_rts_susceptibility_0p5deg.nc')).load()
    inv = xr.open_dataset(os.path.join(sus_dir, 'russia_rts_inventory_count_0p5deg.nc')).rts_count.values
    ws = np.load(os.path.join(sus_dir, 'ws_count.npy')); fobs = np.load(os.path.join(sus_dir, 'darts_fobs.npy'))
    top = float(ds.lat_bnds.values[0, 1]); west = float(ds.lon_bnds.values[0, 0])
    J = match_index(ds.lat.values, [c['lat'] for c in cells]); I = match_index(ds.lon.values, [c['lon'] for c in cells])
    # polygon inventories -> per-cell counts on the aggregation grid
    import geopandas as gpd
    lat0, lat1 = min(c['lat'] for c in cells) - 0.3, max(c['lat'] for c in cells) + 0.3
    lon0, lon1 = min(c['lon'] for c in cells) - 0.3, max(c['lon'] for c in cells) + 0.3
    def grid_points(path, area_col=None, parquet=False):
        df = gpd.read_parquet(path) if parquet else gpd.read_file(path)
        df = df.to_crs(4326); df['geometry'] = df.make_valid(); rp = df.representative_point()
        m = (rp.y > lat0) & (rp.y < lat1) & (rp.x > lon0) & (rp.x < lon1); df, rp = df[m], rp[m]
        rr, cc = point_cell(rp.y.values, rp.x.values, top, west)
        cnt = np.zeros(ds.cell_area.shape); np.add.at(cnt, (rr, cc), 1)
        area = None
        if area_col:
            area = np.zeros(ds.cell_area.shape); np.add.at(area, (rr, cc), df[area_col].values)
        return cnt, area, int(m.sum())
    arts_c, arts_a, n_arts_box = grid_points(os.path.join(sus_dir, 'polygons', 'arts_russia_positive.gpkg'), 'area_ha')
    darts_c, _, n_darts_box = grid_points(os.path.join(sus_dir, 'polygons', 'darts', 'darts_russia_footprints_union.parquet'), parquet=True)
    rows = []
    s5 = dict(scenario='ssp5')
    for c, j, i in zip(cells, J, I):
        r = dict(cell=c['cell'], j=c['j'], i=c['i'], lat=c['lat'], lon=c['lon'], agg_row=int(j), agg_col=int(i))
        for p in PERIODS:
            sel = dict(s5, period=p, lat=ds.lat.values[j], lon=ds.lon.values[i])
            r[f'jia_mean_ssp5_{p}'] = float(ds.rts_susceptibility_mean.sel(sel))
            r[f'jia_gt050_ssp5_{p}'] = float(ds.area_fraction_gt050.sel(sel))
            r[f'jia_gt075_ssp5_{p}'] = float(ds.area_fraction_gt075.sel(sel))
            r[f'valid_fraction_ssp5_{p}'] = float(ds.valid_fraction.sel(sel))
        r['jia_mean_current'] = float(ds.rts_susceptibility_mean_current.values[j, i])
        r['jia_gt050_current'] = float(ds.area_fraction_gt050_current.values[j, i])
        r['jia_gt075_current'] = float(ds.area_fraction_gt075_current.values[j, i])
        r['valid_fraction_current'] = float(ds.valid_fraction_current.values[j, i])
        r['russia_land_fraction'] = float(ds.russia_land_fraction.values[j, i])
        r['cell_area_km2'] = float(ds.cell_area.values[j, i])
        r['nesterova_count'] = float(ws[j, i]); r['jia_inv_count'] = float(inv[j, i])
        r['arts_count'] = float(arts_c[j, i]); r['arts_area_ha'] = float(arts_a[j, i])
        r['darts_fobs'] = float(fobs[j, i])
        r['darts_n_footprints'] = float(darts_c[j, i]) if np.isfinite(fobs[j, i]) else float('nan')
        r['jia_nan_ssp5'] = int(any(np.isnan(r[f'jia_{v}_ssp5_{p}']) for v in ('mean', 'gt050', 'gt075') for p in PERIODS))
        r['jia_nan_current'] = int(any(np.isnan(r[f'jia_{v}_current']) for v in ('mean', 'gt050', 'gt075')))
        rows.append(r)
    meta = dict(n_cells=len(rows), arts_polygons_in_box=n_arts_box, darts_footprints_in_box=n_darts_box,
                anchor_ha=ANCHOR_HA, mf=dict(central=MF_CENTRAL, low=MF_LOW, high=MF_HIGH),
                jia_source=ds.attrs.get('source', ''), agg_top_edge=top)
    return derive_cases(rows, [(r['j'], r['i']) for r in rows]), meta


def write_csv(rows, path):
    import csv
    cols = list(rows[0].keys())
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, cols); w.writeheader()
        for r in rows: w.writerow({k: ('nan' if isinstance(v, float) and np.isnan(v) else v) for k, v in r.items()})


def summarize(rows):
    A = lambda k: np.array([r[k] for r in rows], float)
    ca = A('cell_area_km2'); out = {}
    for k in ('nesterova_count', 'jia_inv_count', 'arts_count'):
        out[k] = dict(total=float(np.nansum(A(k))), cells_zero=int((A(k) == 0).sum()), cells_nan=int(np.isnan(A(k)).sum()))
    out['darts'] = dict(cells_available=int(np.isfinite(A('darts_fobs')).sum()))
    out['jia_nan_cells_ssp5'] = [r['cell'] for r in rows if r['jia_nan_ssp5']]
    out['jia_nan_cells_current'] = [r['cell'] for r in rows if r['jia_nan_current']]
    for x in ('mean', 'gt050'):
        Fs = sum(np.nan_to_num(A(f'case_a_F_{x}_{p}')) for p in PERIODS)
        out[f'case_a_{x}_final_km2'] = dict(central=float((Fs * ca).sum()),
                                            low=float((Fs * ca).sum() * MF_LOW / MF_CENTRAL),
                                            high=float((Fs * ca).sum() * MF_HIGH / MF_CENTRAL),
                                            cells_zero=int((Fs == 0).sum()))
        out[f'case_a_{x}_km2_by_cohort'] = {p: float(np.nansum(A(f'case_a_F_{x}_{p}') * ca)) for p in PERIODS}
    out['case_b_A0_km2'] = float(np.nansum(A('case_b_A0') * ca))
    out['case_b_km2_by_period'] = {p: float(np.nansum(A(f'case_b_A_{p}') * ca)) for p in PERIODS}
    out['case_b_cells_nan'] = {p: int(np.isnan(A(f'case_b_A_{p}')).sum()) for p in PERIODS}
    out['grid_area_km2'] = float(ca.sum())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--grid', required=True); ap.add_argument('--sus-dir', required=True)
    ap.add_argument('-o', '--out', required=True)
    a = ap.parse_args()
    rows, meta = extract(a.grid, a.sus_dir)
    write_csv(rows, a.out)
    meta['summary'] = summarize(rows)
    json.dump(meta, open(os.path.splitext(a.out)[0] + '_meta.json', 'w'), indent=1, default=float)
    print(json.dumps(meta['summary'], indent=1, default=float))


if __name__ == '__main__':
    main()
