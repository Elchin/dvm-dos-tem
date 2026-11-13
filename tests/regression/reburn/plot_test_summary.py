"""One-page summary of the fire-balance test results, from the summary JSONs in FIGDIR:
identity_summary.json (ident.py), n_leak_summary.json, fire_balance_summary.json, overshoot.json,
split_summary.json. Writes FIGDIR/test_summary.png. Options: see common.py (only --figdir is used)."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent)); import common
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
a = common.finish(common.parser(__doc__).parse_args()); F = a.figdir
J = lambda n: json.load(open(F / n))
ident, leak, fb, ov, sp = J("identity_summary.json"), J("n_leak_summary.json"), J("fire_balance_summary.json"), J("overshoot.json"), J("split_summary.json")
fig, ax = plt.subplots(1, 4, figsize=(17, 5.0), dpi=100, gridspec_kw=dict(width_ratios=[1.1, 1.25, 0.7, 1.0]), constrained_layout=True)
# (1) bit identity
a0 = ax[0]; labs = []; ok = []; tot = []
for k, v in ident.items():
    labs.append(k.replace(": ", ":\n")); ok.append(v["identical"]); tot.append(v["total"])
y = np.arange(len(labs))
a0.barh(y, tot, color="0.85", label="stage folders compared"); a0.barh(y, ok, color="tab:green", height=0.5, label="bit-identical")
for i, (o, t) in enumerate(zip(ok, tot)): a0.text(t + 0.5, i, f"{o}/{t}", va="center", fontsize=8)
a0.set_yticks(y); a0.set_yticklabels(labs, fontsize=7.5); a0.invert_yaxis(); a0.set_xlim(0, max(tot) * 1.25)
a0.set_xlabel("stage folders (spin000, main)", fontsize=8)
a0.set_title("Bit identity with the new switches 0\n(all NetCDF outputs + restart files)", fontsize=9); a0.legend(fontsize=7, loc="upper right")
# (2) N balance gap sums
a1 = ax[1]; rows = [("classic, fire fixes off", fb["F1c"]["N_resid_sum"]), ("classic, fire fixes", fb["F1c_on"]["N_resid_sum"])]
for k in ("classic FRI 15", "explicit classic (E1)", "process (Reburn)"):
    rows += [(f"{k}: leak fixes off", leak[k]["N_gap_sum_before"]), (f"{k}: leak fixes on", leak[k]["N_gap_sum_after"])]
y = np.arange(len(rows)); v = [r[1] for r in rows]
a1.barh(y, v, color=["tab:orange" if "off" in r[0] else "k" for r in rows])
for i, x in enumerate(v): a1.text(1.5, i, f"{x:+.2f}", va="center", ha="left", fontsize=7.5)
a1.set_yticks(y); a1.set_yticklabels([r[0] for r in rows], fontsize=7); a1.invert_yaxis(); a1.set_xscale("symlog", linthresh=1)
a1.set_xlim(min(v) * 1.5, 60); a1.axvline(0, color="0.5", lw=.6)
a1.set_title("N balance gap, sum over the 75-yr main stage\n(g N m$^{-2}$; orange = fix off, black = on; symlog)", fontsize=9)
# (3) overshoot
a2 = ax[2]; ks = [("F1c_on", "old restart\n(19af750f)"), ("F1c_on_fb3", "fire state in\nrestart (06c41688)")]
a2.bar([l for _, l in ks], [ov[k]["ratio"] for k, _ in ks], color=["tab:orange", "k"])
for i, (k, _) in enumerate(ks): a2.text(i, ov[k]["ratio"] + 0.02, f"{ov[k]['ratio']:.2f}\n({ov[k]['months_depositing']} months)", ha="center", fontsize=8)
a2.axhline(1, color="0.5", ls="--", lw=.8); a2.set_ylim(0, 1.55); a2.set_ylabel("re-deposited / emitted N")
a2.set_title("Spin-up fire N re-deposited across the\nspin-up → main restart (classic FRI 15)", fontsize=9); a2.tick_params(axis="x", labelsize=7.5)
# (4) restart continuity
a3 = ax[3]; vs = ["DEADC", "DWDC", "DEADN", "BURNAIR2SOILN"]; w = 0.27
for i, (k, lab, c) in enumerate((("E1on_fb3", "fire state in restart", "k"), ("E1on_fb2", "old restart, fixes on", "tab:orange"), ("E1_fb3", "fixes off", "0.6"))):
    a3.bar(np.arange(len(vs)) + (i - 1) * w, [max(sp[k][v]["maxabs"], 1e-9) for v in vs], w, color=c, label=lab)
a3.set_yscale("log"); a3.set_xticks(range(len(vs))); a3.set_xticklabels(vs, fontsize=7.5, rotation=20); a3.legend(fontsize=7, loc="lower left"); a3.set_ylim(1e-9, 300)
a3.set_ylabel("max |restarted − continuous| (g m$^{-2}$)"); a3.set_title("Restart continuity: E1 split in 1985,\n5 yr after a fire (both fire fixes on)", fontsize=9)
fig.suptitle("Fire-balance tests (TVC smoke cases and the 16 regression cases, 20-yr spin-up + 75-yr main); numbers from the committed summary JSONs", fontsize=10)
fig.savefig(F / "test_summary.png"); print("wrote", F / "test_summary.png")
