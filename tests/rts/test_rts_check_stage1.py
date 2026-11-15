"""Tests for scripts/rts_sus2case/check_stage1.py nbal parsing / N residual."""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'scripts'))
from rts_sus2case import check_stage1 as cs  # noqa: E402


def test_parse_and_residual(tmp_path):
    p = tmp_path / 'nbal.txt'
    lines = []
    for stage in range(3):
        n = 100.0
        for yr in range(2):
            for m in range(12):
                # N1 = N0 + avlnin - avlnlost (closed) + 0.001 error in stage 2
                n1 = n + 0.5 - 0.2 + (0.001 if stage == 2 else 0)
                lines.append(f'NB {yr} {m} {n} {n1} 0.5 0 0.2 0.3 0 0 0 0')
                n = n1
    lines.insert(5, 'some other stderr line')
    p.write_text('\n'.join(lines) + '\n')
    a, seg = cs.parse_nbal(str(p))
    assert a.shape == (72, 12) and list(np.bincount(seg)) == [24, 24, 24]
    res, ol = cs.n_residual(a)
    assert np.allclose(res[seg < 2], 0) and np.allclose(res[seg == 2], 0.001)
    assert np.allclose(ol, 0.3)


def test_residual_includes_regrow_n_input():
    # yr m N0 N1 ai oi al ol tl tv vi rn: N1-N0 = ai - tl + rn
    a = np.array([[0, 7, 100.0, 100.0 + 0.5 - 2.0 + 0.09, 0.5, 0, 0, 0, 2.0, 0, 0, 0.09]])
    res, _ = cs.n_residual(a)
    assert abs(res[0]) < 1e-12
