"""Tests for the pure helpers of scripts/rts_sus2case/cell_table.py (no data files needed).

Refs: Jia et al. 2026 Sci. Data doi:10.1038/s41597-026-08284-y (susceptibility);
Nesterova et al. 2025 West Siberia RTS v0.1 (counts); ARTS v1.0 1.44 ha anchor.
"""
import math, os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
from rts_sus2case import cell_table as ct  # noqa: E402


def test_match_index():
    lat = np.arange(81.75, 41.0, -0.5)
    assert list(ct.match_index(lat, [67.75, 72.25])) == [28, 19]
    with pytest.raises(ValueError):
        ct.match_index(lat, [67.5])


def test_point_cell_consistent_with_match_index():
    top = 82.0 - 0.0016
    lat = np.arange(81.75, 41.0, -0.5) - 0.0016; lon = np.arange(-179.75, 180, 0.5)
    r, c = ct.point_cell([67.9, 67.4], [75.9, 75.4], top)
    assert list(r) == list(ct.match_index(lat, [67.75 - 0.0016, 67.25 - 0.0016], tol=0.01))
    assert list(c) == list(ct.match_index(lon, [75.75, 75.25]))


def test_linear_F():
    F = ct.linear_F([0.4, 0.5, 0.75, 1.0, np.nan], 7.5e-4)
    assert F[0] == 0 and F[1] == 0 and math.isclose(F[2], 3.75e-4) and math.isclose(F[3], 7.5e-4)
    assert np.isnan(F[4])


def test_smooth3x3():
    land = np.ones((3, 3), bool); land[0, 0] = False
    f = np.arange(9.0).reshape(3, 3); f[1, 1] = np.nan
    s = ct.smooth3x3(f, land)
    assert np.isnan(s[0, 0])                          # non-land
    assert math.isclose(s[1, 1], np.mean([1, 2, 3, 5, 6, 7, 8]))  # ignores NaN + non-land
    assert math.isclose(s[2, 2], np.mean([5, 7, 8]))


def _rows():
    rows = []
    for j in range(2):
        for i in range(2):
            r = dict(cell=f'c{j}{i}', j=j, i=i, cell_area_km2=1000.0, nesterova_count=10.0 * (j + i),
                     jia_inv_count=0.0, arts_count=1.0, darts_fobs=np.nan,
                     jia_nan_ssp5=0, jia_nan_current=0)
            for k, p in enumerate(ct.PERIODS):
                r[f'jia_mean_ssp5_{p}'] = 0.5 + 0.1 * k
                r[f'jia_gt050_ssp5_{p}'] = 0.4 * (1 + k)
            rows.append(r)
    return rows


def test_derive_cases():
    rows = ct.derive_cases(_rows(), [(0, 0), (0, 1), (1, 0), (1, 1)], shape=(2, 2))
    r = rows[3]
    assert math.isclose(r['case_b_A0'], 20 * 1.44 * 0.01 / 1000)   # count * 1.44 ha / cell area
    assert math.isclose(r['case_b_ratio_s3_2090'], 4.0) and math.isclose(r['case_b_A_2090'], 4 * r['case_b_A0'])
    assert r['case_a_F_mean_2030'] == 0 and math.isclose(r['case_a_F_mean_2090'], 7.5e-4 * 0.3 / 0.5)
    assert r['case_a_F_gt050_2030'] == 0 and math.isclose(r['case_a_F_gt050_2050'], 7.5e-4 * 0.3 / 0.5)
    s = ct.summarize(rows)
    assert s['nesterova_count']['cells_zero'] == 1 and s['case_b_A0_km2'] == pytest.approx(40 * 0.0144)
