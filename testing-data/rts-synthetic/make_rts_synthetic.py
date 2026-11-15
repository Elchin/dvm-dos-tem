#!/usr/bin/env python
"""Regenerate the 1-cell synthetic RTS test input set.

Derived from demo-data Toolik cell (0,0) with controlled overrides:
  * climate: 1901-1930 mean seasonal cycle of that cell, repeated every year
    (no trend, no interannual variability) for historic + projected files
  * veg CMT 5 (tussock tundra; organic-rich soil), drainage 0 (well drained),
    fixed soil texture, flat topography
  * fire: none (fri=2000, explicit burn mask all 0)
  * thermokarst: RTS initiation at transient-year index RTS_YEAR, julian day
    RTS_JDAY. `thermokarst-none.nc` is the rts_init_year=-1 (no-op) variant.
    `projected-explicit-thermokarst-rts.nc` initiates at scenario-year index
    RTS_SC_YEAR (two-stage workflow, scripts/run_rts_two_stage.py).
Ground-ice content is not prescribed directly: it comes from the CMT
initialization plus the cold equilibrium climate (open question).
Usage: python make_rts_synthetic.py [--rts-year N] [--rts-jday D]
"""
import argparse, os, shutil
import numpy as np
import netCDF4 as nc

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'demo-data',
                   'cru-ts40_ar5_rcp85_ncar-ccsm4_toolik_field_station_10x10')

def copy_1x1(fname, out, overrides=None):
    overrides = overrides or {}
    s = nc.Dataset(os.path.join(SRC, fname))
    d = nc.Dataset(os.path.join(HERE, out), 'w', format='NETCDF4')
    for dn, dim in s.dimensions.items():
        size = 1 if dn in ('X', 'Y') else len(dim)
        d.createDimension(dn, None if dim.isunlimited() else size)
    for vn, v in s.variables.items():
        fill = v.getncattr('_FillValue') if '_FillValue' in v.ncattrs() else None
        nv = d.createVariable(vn, v.dtype, v.dimensions, fill_value=fill)
        nv.setncatts({a: v.getncattr(a) for a in v.ncattrs() if a != '_FillValue'})
        sl = tuple(slice(0, 1) if dn in ('X', 'Y') else slice(None) for dn in v.dimensions)
        data = v[sl] if v.dimensions else v[...]
        if vn in overrides:
            data = overrides[vn](np.array(data))
        nv[...] = data
    d.close(); s.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rts-year', type=int, default=10, help='transient year index of RTS initiation (-1: none)')
    ap.add_argument('--rts-jday', type=int, default=213)
    ap.add_argument('--rts-sc-year', type=int, default=5,
                    help='scenario-year index of RTS initiation for the two-stage workflow file')
    a = ap.parse_args()

    # climate: repeated mean annual cycle
    h = nc.Dataset(os.path.join(SRC, 'historic-climate.nc'))
    cyc = {k: h[k][:360, 0, 0].reshape(30, 12).mean(axis=0) for k in ('tair', 'precip', 'nirr', 'vapor_press')}
    h.close()
    clim = {k: (lambda arr, c=c: np.tile(c, arr.shape[0] // 12).reshape(arr.shape)) for k, c in cyc.items()}
    copy_1x1('historic-climate.nc', 'historic-climate.nc', clim)
    copy_1x1('projected-climate.nc', 'projected-climate.nc', clim)
    for f in ('co2.nc', 'projected-co2.nc'):
        shutil.copy(os.path.join(SRC, f), os.path.join(HERE, f))
    one = lambda val: (lambda arr: np.full(arr.shape, val, dtype=arr.dtype))
    copy_1x1('run-mask.nc', 'run-mask.nc', {'run': one(1)})
    copy_1x1('vegetation.nc', 'vegetation.nc', {'veg_class': one(5)})
    copy_1x1('drainage.nc', 'drainage.nc', {'drainage_class': one(0)})
    copy_1x1('soil-texture.nc', 'soil-texture.nc',
             {'pct_sand': one(30.), 'pct_silt': one(50.), 'pct_clay': one(20.)})
    copy_1x1('topo.nc', 'topo.nc', {'slope': one(1.), 'aspect': one(0.), 'elevation': one(500.)})
    copy_1x1('fri-fire.nc', 'fri-fire.nc', {'fri': one(2000), 'fri_severity': one(1),
             'fri_jday_of_burn': one(200), 'fri_area_of_burn': one(1)})
    nofire = {'exp_burn_mask': one(0), 'exp_jday_of_burn': one(200),
              'exp_fire_severity': one(1), 'exp_area_of_burn': one(1)}
    copy_1x1('historic-explicit-fire.nc', 'historic-explicit-fire.nc', nofire)
    copy_1x1('projected-explicit-fire.nc', 'projected-explicit-fire.nc', nofire)

    def tk(year):
        def mask(arr):
            m = np.zeros(arr.shape, dtype=arr.dtype)
            if year >= 0: m[year, ...] = 1
            return m
        return {'exp_thermokarst_mask': mask, 'exp_jday_of_thermokarst': one(a.rts_jday),
                'exp_thermokarst_severity': one(1), 'exp_area_of_thermokarst': one(1)}
    copy_1x1('historic-explicit-thermokarst.nc', 'historic-explicit-thermokarst.nc', tk(a.rts_year))
    copy_1x1('historic-explicit-thermokarst.nc', 'thermokarst-none.nc', tk(-1))
    copy_1x1('projected-explicit-thermokarst.nc', 'projected-explicit-thermokarst.nc', tk(-1))
    # two-stage workflow: RTS in the scenario (sc) stage after a short spin-up
    copy_1x1('projected-explicit-thermokarst.nc', 'projected-explicit-thermokarst-rts.nc', tk(a.rts_sc_year))
    print('wrote synthetic 1-cell RTS inputs to', HERE, 'rts_year =', a.rts_year)

if __name__ == '__main__':
    main()
