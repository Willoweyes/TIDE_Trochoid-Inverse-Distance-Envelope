"""Fig. 4 -- (a) one cutting condition per case run with 64 numerical settings
of the forward solution method (cells, step ratio r, edge points) and with TIDE
over sampling only: runtime against signed S_a error; (b) runtime against mean
|error| along the doubling ladder of cells (time step tied to the grid),
averaged over A1-A3; the error axis is reversed so that better is to the right.
Reference: TIDE at 1600 cells.  LNCS text width."""
import sys, csv
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import style
plt = style.apply()
from matplotlib.lines import Line2D
ROOT = Path(__file__).resolve().parents[2]
CT, CF = style.C_EXACT, style.C_FORWARD
CASES = ["A1", "A2", "A3"]; MK = dict(A1="o", A2="s", A3="^")
fig, (a, b) = plt.subplots(1, 2, figsize=(style.W_FULL, 1.82))

# ---------------- (a) settings cloud ----------------
R = list(csv.DictReader(open(ROOT/"data/raw/fsm_settings_cloud.csv")))
shade = {50: "#f3b395", 100: "#e07b4f", 200: "#c1440e", 400: "#7a2a08"}
for c in CASES:
    ref = float([r for r in R if r["case"] == c and r["method"] == "tide" and r["grid"] == "1600"][0]["Sa"])
    for g in shade:
        F = [r for r in R if r["case"] == c and r["method"] == "fsm" and int(r["grid"]) == g]
        a.scatter([100*(float(r["Sa"])/ref - 1) for r in F], [float(r["t_s"]) for r in F], s=7,
                  marker=MK[c], color=shade[g], lw=0, alpha=0.9)
    T_ = sorted([r for r in R if r["case"] == c and r["method"] == "tide" and r["grid"] != "1600"],
                key=lambda r: int(r["grid"]))
    a.plot([100*(float(r["Sa"])/ref - 1) for r in T_], [float(r["t_s"]) for r in T_], "-", marker=MK[c],
           ms=3.2, color=CT, mfc="white", lw=0.9, zorder=5)
a.axvline(0, color="0.4", lw=0.6, ls=":")
a.set_yscale("log"); a.set_xlim(-40, 8); a.set_ylim(1e-3, 1e2)
a.set_xlabel(r"$S_a$ error [%]"); a.set_ylabel("runtime per surface [s]")
a.set_title("(a) one condition, many settings", loc="left")
hand = [Line2D([], [], lw=0, marker="o", ms=3.4, color=shade[g]) for g in shade] + \
       [Line2D([], [], color=CT, lw=0.9, marker="o", ms=3.2, mfc="white")]
a.legend(hand, [f"FSM, {g} cells" for g in shade] + ["TIDE, 50–400 cells"], loc="upper left",
         fontsize=6.0, handlelength=1.4, borderpad=0.3, labelspacing=0.25)

# ---------------- (b) doubling ladder ----------------
D = list(csv.DictReader(open(ROOT/"data/raw/cost_accuracy_dense.csv")))
G = [25, 50, 100, 200, 400, 800, 1600]
ref = {c: {k: float([r for r in D if r["case"] == c and r["grid_base"] == "1600"][0]["tide_" + k])
           for k in ("Sa", "Sz")} for c in CASES}
def path(meth, key):
    e, t = [], []
    for g in G:
        rs = [r for r in D if int(r["grid_base"]) == g]
        e.append(np.mean([abs(float(r[f"{meth}_{key}"])/ref[r["case"]][key] - 1)*100 for r in rs]))
        t.append(np.mean([float(r[f"{meth}_t"]) for r in rs]))
    return np.array(e), np.array(t)
reach = {}
for meth, col, mfc in (("tide", CT, "white"), ("fsm", CF, None)):
    for key, ls in (("Sa", "-"), ("Sz", "--")):
        e, t = path(meth, key); k = e > 0
        b.loglog(e[k], t[k], ls=ls, marker="o", ms=3.0, color=col,
                 mfc=mfc if mfc else col, lw=1.0)
        i = next(i for i in range(len(G)) if np.all(e[i:] <= 1.0))
        reach[(meth, key)] = (G[i], t[i])
        gl = [G[j] for j in range(len(G)) if k[j]]
        if (meth, key) in (("fsm", "Sz"), ("tide", "Sa")):
            off, va, ha = ((0, 5), "bottom", "center") if meth == "fsm" else ((5, -1), "center", "left")
            b.annotate(f"{gl[-1]} cells", (e[k][-1], t[k][-1]), xytext=off, textcoords="offset points",
                       fontsize=6, color=col, ha=ha, va=va)
tT = reach[("tide", "Sa")][1]; tF = reach[("fsm", "Sa")][1]
assert reach[("tide", "Sz")][0] == reach[("tide", "Sa")][0] and reach[("fsm", "Sz")][0] == reach[("fsm", "Sa")][0]
b.axvline(1.0, color="0.4", lw=0.6, ls=":")
def fmt_t(t):     # same rounding as the numbers quoted in the text (paper/compute_numbers.py)
    return f"{t:.3f}" if t < 0.1 else (f"{t:.2f}" if t < 1 else (f"{t:.1f}" if t < 10 else f"{t:.0f}"))
for tt, col, lab, xl, ha, fy, va in ((tT, CT, f"TIDE at 1%: {fmt_t(tT)} s", 0.0058, "right", 1/1.3, "top"),
                                     (tF, CF, f"FSM at 1%: {fmt_t(tF)} s", 90, "left", 1.25, "bottom")):
    b.axhline(tt, color=col, lw=0.6, ls=(0, (4, 2)))
    b.text(xl, tt*fy, lab, fontsize=6.0, color=col, ha=ha, va=va)
b.set_xlim(100, 0.0015); b.set_ylim(1e-3, 2e2)
b.set_xlabel("mean relative error, A1–A3 [%]"); b.set_ylabel("mean runtime [s]")
b.set_title("(b) doubling the cells", loc="left")
hand = [Line2D([], [], color="0.3", lw=1.0), Line2D([], [], color="0.3", lw=1.0, ls="--"),
        Line2D([], [], color=CT, lw=1.0), Line2D([], [], color=CF, lw=1.0)]
b.legend(hand, ["$S_a$", "$S_z$", "TIDE", "FSM"], loc="lower right", fontsize=6.0, ncol=2, handlelength=1.8,
         columnspacing=0.8, borderpad=0.3)
for ax in (a, b): ax.minorticks_off()
fig.tight_layout(w_pad=1.0)
style.save(fig, ROOT, "fig_cost")
print("reach", reach, "ratio %.0f" % (tF/tT))
