"""Tests for scripts/rts_sus2case/rts_trigger_file.py (single-year RTS trigger writer).

The fixture tests/rts/data/c55_projected-explicit-thermokarst-rts.nc is the Gydan
10x10 cell c55 file from spatial10x10/build_grid.py (2021, jday 213).
Set RTS_TRIGGER_INPUTS=<spatial10x10/inputs> to also check all 90 real cells.
"""
import filecmp, glob, os, sys
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts'))
from rts_sus2case import rts_trigger_file as rtf  # noqa: E402

FIX = os.path.join(HERE, 'data', 'c55_projected-explicit-thermokarst-rts.nc')
COHORTS = list(range(2021, 2100, 10))


def test_fixture_is_2021_august():
    d = rtf.read_trigger(FIX)
    assert rtf.initiation_years(FIX) == [2021]
    assert len(d[rtf.MASK]) == 85
    assert set(d['exp_jday_of_thermokarst']) == {213}
    assert rtf.doy2month(213) == 7  # August, 0-based


@pytest.mark.parametrize('doy,month', [(0, 0), (30, 0), (31, 1), (211, 6), (212, 7), (242, 7), (243, 8), (333, 10), (334, 11), (364, 11)])
def test_doy2month_matches_cpp(doy, month):
    assert rtf.doy2month(doy) == month


def test_2021_byte_identical(tmp_path):
    out = tmp_path / 'tk.nc'
    assert rtf.write_trigger(FIX, str(out), 2021) is False
    assert filecmp.cmp(FIX, out, shallow=False)


@pytest.mark.parametrize('year', COHORTS)
def test_round_trip(tmp_path, year):
    out = tmp_path / f'tk{year}.nc'
    rtf.write_trigger(FIX, str(out), year)
    assert rtf.initiation_years(str(out)) == [year]
    a, b = rtf.read_trigger(FIX), rtf.read_trigger(str(out))
    for v in rtf.TK_VARS[1:]:
        np.testing.assert_array_equal(a[v], b[v])          # severity / jday / area untouched
    assert b[rtf.MASK].sum() == 1 and b[rtf.MASK][year - 2016] == 1


def test_other_vars_and_coords_preserved(tmp_path):
    import netCDF4
    out = tmp_path / 'tk2051.nc'
    rtf.write_trigger(FIX, str(out), 2051)
    with netCDF4.Dataset(FIX) as a, netCDF4.Dataset(out) as b:
        assert set(a.variables) == set(b.variables)
        for v in a.variables:
            if v == rtf.MASK:
                continue
            assert a[v].dtype == b[v].dtype and a[v].dimensions == b[v].dimensions
            np.testing.assert_array_equal(np.asarray(a[v][:]), np.asarray(b[v][:]))
            assert a[v].ncattrs() == b[v].ncattrs()
        assert 'initiation year 2051' in b.source


def test_rewrite_back_to_2021_value_identical(tmp_path):
    a, b = tmp_path / 'a.nc', tmp_path / 'b.nc'
    rtf.write_trigger(FIX, str(a), 2061)
    assert rtf.write_trigger(str(a), str(b), 2021) is True     # forces an r+ rewrite
    x, y = rtf.read_trigger(FIX), rtf.read_trigger(str(b))
    for v in rtf.TK_VARS:
        np.testing.assert_array_equal(x[v], y[v])


def test_no_trigger(tmp_path):
    out = tmp_path / 'none.nc'
    assert rtf.write_trigger(FIX, str(out), None) is True
    assert rtf.initiation_years(str(out)) == []


@pytest.mark.parametrize('year', [2015, 2016, 2020, 2101])
def test_invalid_years(tmp_path, year):
    with pytest.raises(ValueError):
        rtf.write_trigger(FIX, str(tmp_path / 'x.nc'), year)


def test_template_not_overwritten():
    with pytest.raises(ValueError):
        rtf.write_trigger(FIX, FIX, 2031)


@pytest.mark.skipif(not os.environ.get('RTS_TRIGGER_INPUTS'), reason='RTS_TRIGGER_INPUTS not set')
def test_all_real_cells(tmp_path):
    files = sorted(glob.glob(os.path.join(os.environ['RTS_TRIGGER_INPUTS'], 'c*', 'projected-explicit-thermokarst-rts.nc')))
    assert len(files) == 90
    for f in files:
        out = tmp_path / 'tk.nc'
        assert rtf.write_trigger(f, str(out), 2021) is False and filecmp.cmp(f, out, shallow=False)
        for y in (2031, 2091):
            rtf.write_trigger(f, str(out), y)
            assert rtf.initiation_years(str(out)) == [y]
