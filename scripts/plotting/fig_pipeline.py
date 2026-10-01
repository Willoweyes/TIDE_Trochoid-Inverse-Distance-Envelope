"""Fig. 5 -- hybrid roughness model on 36 conditions around A1: RMS error of
(a) the prediction M + delta_hat and (b) the learned discrepancy delta_hat,
against the wall-clock time of the whole model (36 simulations + GP fit).
Mean and standard deviation over 20 splits.  LNCS text width."""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import style
plt = style.apply()
from matplotlib.lines import Line2D
ROOT = Path(__file__).resolve().parents[2]
CT, CF = style.C_EXACT, style.C_FORWARD
P = json.loads((ROOT/"data/processed/pipeline_study.json").read_text())
NOISE = 0.02

def short(name):
    s = name.replace("TIDE-", "").replace("FSM-", "")
    return (s.replace(" released", ", released").replace(" 8n r0.25", ", $8n$,\n$r$ = 0.25").replace(" r0.9", ", $r$ = 0.9")
             .replace(" 1n", ", $1n$"))

# label offsets (points) for panel (b); tuned by eye to avoid overlaps
OFF = {"TIDE-25": (5, 0, "left", "center"), "TIDE-50": (5, -2, "left", "top"),
       "TIDE-100": (5, 0, "left", "center"), "TIDE-200": (5, 0, "left", "center"),
       "FSM-50": (5, 0, "left", "center"), "FSM-100": (5, 0, "left", "center"),
       "FSM-200": (5, 0, "left", "center"), "FSM-400": (-5, 0, "right", "center"),
       "FSM-200 released": (5, 0, "left", "center"), "FSM-100 1n": (-5, 0, "right", "center"),
       "FSM-200 r0.9": (-5, 4.5, "right", "center"), "FSM-400 8n r0.25": (5, 0, "left", "center")}
if (ROOT/"scripts/plotting/pipeline_offsets.json").exists():
    OFF.update({k: tuple(v) for k, v in json.loads((ROOT/"scripts/plotting/pipeline_offsets.json").read_text()).items()})

fig, (a, b) = plt.subplots(1, 2, figsize=(style.W_FULL, 1.72), sharey=True)
for ax, key, sd, ttl in ((a, "pred_err", None, r"(a) prediction $M + \hat\delta$"),
                         (b, "delta_err", "delta_err_sd", r"(b) learned discrepancy $\hat\delta$")):
    for name, v in P.items():
        tide = name.startswith("TIDE")
        col = CT if tide else CF
        yerr = v[sd] if sd else None
        ax.errorbar(v["t_pipeline"], v[key], yerr=yerr, fmt="o", ms=3.4, color=col,
                    mfc="white" if tide else col, mew=0.8, lw=0.7, capsize=1.4, zorder=5)
        if ax is b or name == "FSM-100 1n":
            dx, dy, ha, va = OFF.get(name, (5, 0, "left", "center"))
            ax.annotate(short(name), (v["t_pipeline"], v[key]), xytext=(dx, dy), textcoords="offset points",
                        fontsize=6.0, color=col, ha=ha, va=va, zorder=6)
    ax.axhline(NOISE, color="0.35", lw=0.7, ls=":", zorder=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("wall-clock time of the hybrid model [s]")
    ax.set_title(ttl, loc="left")
    ax.minorticks_off()
a.set_ylabel("RMS error [µm]")
b.tick_params(labelleft=True)
a.text(0.97, NOISE*1.12, "noise level", transform=a.get_yaxis_transform(), fontsize=6.0, color="0.3",
       ha="right", va="bottom")
hand = [Line2D([], [], lw=0, marker="o", ms=3.4, color=CT, mfc="white", mew=0.8),
        Line2D([], [], lw=0, marker="o", ms=3.4, color=CF)]
a.legend(hand, ["TIDE", "FSM"], loc="upper left", fontsize=6.2, handletextpad=0.3, borderpad=0.3)
ts = [v["t_pipeline"] for v in P.values()]
for ax in (a, b):
    ax.set_xlim(min(ts) / 2.2, max(ts) * 7.5)
ys = [v["delta_err"] + v["delta_err_sd"] for v in P.values()] + [v["pred_err"] for v in P.values()]
yl = [v["delta_err"] - v["delta_err_sd"] for v in P.values()] + [v["pred_err"] for v in P.values()]
a.set_ylim(min(0.008, min(yl) / 1.3), max(ys) * 1.35)
fig.tight_layout(w_pad=0.8)
style.save(fig, ROOT, "fig_pipeline")
print("ok")
