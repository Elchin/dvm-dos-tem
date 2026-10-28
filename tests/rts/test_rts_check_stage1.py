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
