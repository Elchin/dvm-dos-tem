#!/usr/bin/env python
"""Compare two dvmdostem output directories variable-by-variable.

  python tests/rts/compare_outputs.py DIR_A DIR_B [--rtol 0] [--atol 0]
Exit code 0 if all files/variables match within tolerance.
"""
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rtslib as r

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('a'); ap.add_argument('b')
    ap.add_argument('--rtol', type=float, default=0.0)
    ap.add_argument('--atol', type=float, default=0.0)
    a = ap.parse_args()
    bad = r.compare_dirs(a.a, a.b, a.rtol, a.atol)
    print('\n'.join(bad) if bad else 'MATCH (rtol=%g, atol=%g)' % (a.rtol, a.atol))
    sys.exit(1 if bad else 0)
