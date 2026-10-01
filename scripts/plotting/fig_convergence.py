"""Fig. 2 -- grid refinement of the forward solution method (time step tied to
the grid, r = 0.5) towards TIDE on A1-A3, three panels side by side:
(a) S_a, (b) S_z, (c) relative S_a error of FSM against the reference (TIDE at
1600 cells) on log-log axes with a slope -1 guide.  Colour = method, marker =
case; top axis of (a) and (b) gives the cell size.  LNCS text width."""
import sys, csv
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import style
plt = style.apply()
from matplotlib.lines import Line2D
ROOT = Path(__file__).resolve().parents[2]
rows = list(csv.DictReader(open(ROOT/"data/raw/ladder_dt_regimes.csv")))
CASES = ["A1", "A2", "A3"]; MK = dict(A1="o", A2="s", A3="^")
CT, CF = style.C_EXACT, style.C_FORWARD
LY_UM = 5000.0
fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(style.W_FULL, 1.75))

def sel(c):
    return sorted([r for r in rows if r["case"] == c and r["regime"] == "matched"],
                  key=lambda r: int(r["grid_base"]))

for c in CASES:
    m = sel(c); g = [int(r["grid_base"]) for r in m]
    for ax, key in ((a0, "Sa"), (a1, "Sz")):
        ax.semilogx(g, [float(r["fsm_" + key]) for r in m], marker=MK[c], color=CF, lw=0.9, ms=2.8)
        ax.semilogx(g, [float(r["tide_" + key]) for r in m], ls="--", lw=0.8, marker=MK[c],
                    mfc="white", ms=3.0, color=CT)
for ax, ttl, yl in ((a0, "(a) $S_a$", "$S_a$ [µm]"), (a1, "(b) $S_z$", "$S_z$ [µm]")):
    ax.set_title(ttl, loc="left", pad=3)
    ax.set_ylabel(yl, labelpad=1); ax.set_xlabel("cells $n$", labelpad=1)
    ax.set_xticks([25, 100, 400, 1600]); ax.set_xticklabels(["25", "100", "400", "1600"])
    ax.minorticks_off(); ax.tick_params(labelsize=6.2, pad=1.5)
    top = ax.secondary_xaxis("top", functions=(lambda n: LY_UM / n, lambda d: LY_UM / d))
    top.set_xticks([200, 50, 12.5, 3.125]); top.set_xticklabels(["200", "50", "12.5", "3.1"])
    top.tick_params(labelsize=5.8, pad=1); top.minorticks_off()
    top.set_xlabel("$\\Delta y$ [µm]", fontsize=6.2, labelpad=1.5)
a0.set_ylim(0.65, 2.2); a1.set_ylim(4.2, 10.6)
a0.set_yticks([1.0, 1.5, 2.0]); a1.set_yticks([5, 6, 7, 8, 9, 10])

# (c) S_a error relative to the reference (TIDE at 1600 cells), log-log
orders = []
for c in CASES:
    m = sel(c); g = np.array([int(r["grid_base"]) for r in m])
    ref = float([r for r in m if r["grid_base"] == "1600"][0]["tide_Sa"])
    k = g < 1600
    e = np.array([abs(float(r["fsm_Sa"]) / ref - 1) * 100 for r in m])
    a2.loglog(g[k], e[k], marker=MK[c], ms=2.8, lw=0.9, color=CF)
    kk = (g >= 100) & (g <= 800)
    orders.append(-np.polyfit(np.log(g[kk]), np.log(e[kk]), 1)[0])
a2.loglog([50, 800], [40, 2.5], "k--", lw=0.7)
a2.text(230, 12, "slope $-1$", fontsize=6.0, rotation=-28, ha="center", va="bottom")
a2.set_title("(c) error of FSM in $S_a$", loc="left", pad=3)
a2.set_ylabel("error [%]", labelpad=1); a2.set_xlabel("cells $n$", labelpad=1)
a2.set_xticks([25, 100, 400]); a2.set_xticklabels(["25", "100", "400"])
a2.minorticks_off(); a2.tick_params(labelsize=6.2, pad=1.5)
a2.set_xlim(18, 1100)
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
a2.yaxis.set_major_locator(FixedLocator([1, 3, 10, 30])); a2.yaxis.set_major_formatter(FixedFormatter(["1", "3", "10", "30"]))
a2.yaxis.set_minor_locator(NullLocator())
top = a2.secondary_xaxis("top", functions=(lambda n: LY_UM / n, lambda d: LY_UM / d))
top.set_xticks([200, 50, 12.5]); top.set_xticklabels(["200", "50", "12.5"])
top.tick_params(labelsize=5.8, pad=1); top.minorticks_off()
top.set_xlabel("$\\Delta y$ [µm]", fontsize=6.2, labelpad=1.5)

hand = [Line2D([], [], color=CF, lw=0.9), Line2D([], [], color=CT, lw=0.8, ls="--")] + \
       [Line2D([], [], color="0.3", lw=0, marker=MK[c], ms=3.0, mfc="white") for c in CASES]
a1.legend(hand, ["FSM", "TIDE"] + CASES, loc="lower right", fontsize=5.6, ncol=2,
          columnspacing=0.6, handlelength=1.4, borderpad=0.3, labelspacing=0.2)
print("fitted orders 100-800 cells:", np.round(orders, 3))
fig.tight_layout(w_pad=0.6)
style.save(fig, ROOT, "fig_convergence")
print("ok")
