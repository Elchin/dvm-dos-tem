#!/usr/bin/env python3
"""Write projected-explicit-thermokarst-rts.nc for ONE RTS initiation year per cell.

Phase 1 of the RTS susceptibility two-case design (Case A Jia et al. 2026 Sci. Data
doi:10.1038/s41597-026-08284-y; Case B Nesterova et al. 2025 West Siberia RTS v0.1;
ARTS 1.44 ha anchor). The cohort-response design runs 1 baseline + 8 single-start-year
RTS runs per cell (2021, 2031, ..., 2091) and blends them offline, so each RTS run
needs a trigger file with exactly one initiation year (no repeat slumps).

File layout (as written by spatial10x10/build_grid.py and read by Thermokarst.cpp
L49-69 via temutil::get_timeseries): time(85, years since 2016) x Y(1) x X(1) int vars
  exp_thermokarst_mask      1 in the initiation year, else 0
  exp_thermokarst_severity  1 (kept from the template)
  exp_jday_of_thermokarst   213 (0-based doy -> doy2month = 7 = August; kept)
  exp_area_of_thermokarst   1 (kept)
Thermokarst::should_initiate fires in tr/sc runs when mask[yr]==1 and
doy2month(jday)==month, with yr indexed from 2016 (2021 -> index 5).

Approach: copy the cell's existing (2021) file as the template, so lat/lon, the CRS
variable and all attributes are preserved, then rewrite only the mask. If the new
mask equals the template mask (e.g. year 2021) the copy is left untouched, so the
result is byte-identical to the template. Otherwise the global `source` attribute
is extended with the initiation year.

Earliest valid year is 2021: the regrow GPP ring needs 60 months of sc-run history
and is not carried in the restart. year=None writes an all-zero mask (no trigger).
"""
import argparse, shutil
import numpy as np
import netCDF4

FIRST_YEAR = 2016          # time origin of the sc-run trigger files
EARLIEST_TRIGGER = 2021    # regrow GPP ring (60 months) is not in the restart
MASK = 'exp_thermokarst_mask'
TK_VARS = (MASK, 'exp_thermokarst_severity', 'exp_jday_of_thermokarst', 'exp_area_of_thermokarst')


def doy2month(doy):
    """Mirror of temutil::doy2month (0-based doy and month)."""
    first = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
    if not 0 <= doy <= 364:
        raise ValueError('doy must be in [0, 364]')
    return max(m for m, d in enumerate(first) if doy >= d)


def mask_for_year(year, ntime):
    m = np.zeros(ntime, dtype=np.int32)
    if year is None:
        return m
    last = FIRST_YEAR + ntime - 1
    if not EARLIEST_TRIGGER <= int(year) <= last:
        raise ValueError(f'initiation year {year} outside [{EARLIEST_TRIGGER}, {last}]')
    m[int(year) - FIRST_YEAR] = 1
    return m


def read_trigger(path):
    """Read the four series at (Y=0, X=0) the way get_timeseries does (plain ints)."""
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        return {v: np.asarray(ds.variables[v][:, 0, 0]) for v in TK_VARS}


def initiation_years(path):
    return [FIRST_YEAR + int(k) for k in np.nonzero(read_trigger(path)[MASK] == 1)[0]]


def write_trigger(template, out, year):
    """Copy `template` to `out` and set a single-year mask. Returns True if the mask changed."""
    if template == out:
        raise ValueError('out must differ from template')
    shutil.copyfile(template, out)
    with netCDF4.Dataset(out) as ds:
        ntime = len(ds.dimensions['time'])
        cur = np.asarray(ds.variables[MASK][:, 0, 0])
    new = mask_for_year(year, ntime)
    if np.array_equal(cur, new):
        return False
    with netCDF4.Dataset(out, 'r+') as ds:
        v = ds.variables[MASK]
        v[:] = new.astype(v.dtype).reshape(-1, 1, 1)
        tag = f'initiation year {year}' if year is not None else 'no initiation'
        ds.setncattr('source', f"{getattr(ds, 'source', '')}; mask rewritten by "
                     f"scripts/rts_sus2case/rts_trigger_file.py ({tag})")
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('template'); ap.add_argument('out')
    ap.add_argument('--year', type=int, help='initiation year (omit for no trigger)')
    a = ap.parse_args(argv)
    changed = write_trigger(a.template, a.out, a.year)
    print(f'{a.out}: years {initiation_years(a.out)} ({"rewritten" if changed else "identical to template"})')


if __name__ == '__main__':
    main()
