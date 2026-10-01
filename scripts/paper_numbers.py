# -*- coding: utf-8 -*-
"""Collect the numbers quoted in the paper from data/ -> results/paper_numbers.json and .txt"""
import csv, json, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW, PROC, RES = ROOT / "data/raw", ROOT / "data/processed", ROOT / "results"
RES.mkdir(exist_ok=True)
N, L = {}, []
def P(s): L.append(s); print(s)
def fmt_t(t):
    return f"{t:.3f}" if t < 0.1 else (f"{t:.2f}" if t < 1 else (f"{t:.1f}" if t < 10 else f"{t:.0f}"))

mi = RES / "machine_info.json"
if mi.exists():
    m = json.loads(mi.read_text())
    P(f"Machine: {m.get('cpu')} | {m.get('physical_cores')} cores/{m.get('logical_cores')} threads | "
      f"{m.get('ram_GB')} GB RAM | {m.get('os')} | compiler {m.get('compiler')} | Python {m.get('python')}")

# ---- Fig. 4b: 1% reach on doubled grids, errors relative to TIDE at 1600 cells (mean over A1-A3)
D = list(csv.DictReader(open(RAW / "cost_accuracy_dense.csv")))
C = ["A1", "A2", "A3"]; G = [25, 50, 100, 200, 400, 800, 1600]
ref = {c: {k: float([r for r in D if r["case"] == c and r["grid_base"] == "1600"][0]["tide_" + k])
           for k in ("Sa", "Sz")} for c in C}
reach = {}
for meth in ("tide", "fsm"):
    e = {k: [] for k in ("Sa", "Sz")}; t = []
    for g in G:
        rs = [r for r in D if int(r["grid_base"]) == g]
        for k in ("Sa", "Sz"):
            e[k].append(np.mean([abs(float(r[f"{meth}_{k}"]) / ref[r["case"]][k] - 1) * 100 for r in rs]))
        t.append(np.mean([float(r[f"{meth}_t"]) for r in rs]))
    both = [max(a, b) for a, b in zip(e["Sa"], e["Sz"])]
    i = next(i for i in range(len(G)) if all(x <= 1 for x in both[i:]))
    reach[meth] = (G[i], t[i])
    P(f"{meth}: mean |err| Sa % {dict(zip(G, [round(float(x), 3) for x in e['Sa']]))}")
    P(f"{meth}: mean |err| Sz % {dict(zip(G, [round(float(x), 3) for x in e['Sz']]))}")
    P(f"{meth}: mean runtime s {dict(zip(G, [round(float(x), 4) for x in t]))}")
ratio = reach["fsm"][1] / reach["tide"][1]
N.update(tT1=fmt_t(reach["tide"][1]), tF1=fmt_t(reach["fsm"][1]), grid_T1=reach["tide"][0],
         grid_F1=reach["fsm"][0], ratio=round(ratio, 1), ratio_floor=int(ratio // 50 * 50))
P(f"1% in Sa and Sz: TIDE {reach['tide'][0]} cells {N['tT1']} s, FSM {reach['fsm'][0]} cells {N['tF1']} s, "
  f"ratio {ratio:.1f} (> {N['ratio_floor']})")

# ---- Sect. 5.3: FSM error at 100 cells and convergence order over 100-800 cells
e100 = [100 * (float(r["fsm_Sa"]) / ref[r["case"]]["Sa"] - 1) for r in D if r["grid_base"] == "100"]
P(f"FSM Sa at 100 cells vs TIDE at 1600 cells: {np.round(e100, 2)} %")
o = {"Sa": [], "Sz": []}
for c in C:
    for k in o:
        gg = [100, 200, 400, 800]
        err = [abs(float([r for r in D if r["case"] == c and int(r["grid_base"]) == g][0][f"fsm_{k}"]) / ref[c][k] - 1)
               for g in gg]
        o[k].append(-np.polyfit(np.log(gg), np.log(err), 1)[0])
N["ord_sa"] = f"{min(o['Sa']):.2f}-{max(o['Sa']):.2f}"; N["ord_sz"] = f"{min(o['Sz']):.2f}-{max(o['Sz']):.2f}"
P(f"convergence order Sa {N['ord_sa']}, Sz {N['ord_sz']}")
lad = RAW / "ladder_dt_regimes.csv"
if lad.exists():
    for r in csv.DictReader(open(lad)):
        if r["grid_base"] == "200" and "released" in r.get("regime", ""):
            P(f"released setting {r['case']}: FSM Sa {100 * (float(r['fsm_Sa']) / ref[r['case']]['Sa'] - 1):.2f} % vs TIDE 1600")
ts = PROC / "timestep_check.json"
if ts.exists():
    P("16-fold smaller dt at 100 cells: " + ", ".join(f"{c} {v['change_pct']:.2f} %" for c, v in json.loads(ts.read_text()).items()))

# ---- Sect. 5.1: runtime spread over v_c and f_z (200 cells)
S = list(csv.DictReader(open(RAW / "scaling.csv")))
for sw, key in (("speed", "spread_vc"), ("feed", "spread_fz")):
    t = np.array([float(r["t_s"]) for r in S if r["sweep"] == sw])
    N[key] = f"{100 * (t.max() / t.min() - 1):.1f}"
P(f"TIDE runtime spread: {N['spread_vc']} % over v_c, {N['spread_fz']} % over f_z")

# ---- Fig. 4a: FSM over 64 settings
R = list(csv.DictReader(open(RAW / "fsm_settings_cloud.csv")))
for c in C:
    rf = float([r for r in R if r["case"] == c and r["method"] == "tide" and r["grid"] == "1600"][0]["Sa"])
    ee = [100 * (float(r["Sa"]) / rf - 1) for r in R if r["case"] == c and r["method"] == "fsm"]
    P(f"settings cloud {c}: Sa error {min(ee):.1f} .. {max(ee):.1f} %")
szA1 = [float(r["Sz"]) for r in R if r["case"] == "A1" and r["method"] == "fsm"]
P(f"settings cloud A1: Sz {min(szA1):.2f} .. {max(szA1):.2f} um")

# ---- Sect. 6: hybrid model
Pj = json.loads((PROC / "pipeline_study.json").read_text())
ref36 = {(float(r["fz"]), float(r["vc"]), float(r["ap"])): float(r["Sa_tide1600"])
         for r in csv.DictReader(open(RAW / "pipeline_reference.csv"))}
CR = list(csv.DictReader(open(RAW / "pipeline_cores.csv")))
for name, v in Pj.items():
    d = np.array([float(r["Sa"]) - ref36[(float(r["fz"]), float(r["vc"]), float(r["ap"]))] for r in CR if r["core"] == name])
    P(f"{name:18s} prediction {v['pred_err']:.4f} | discrepancy {v['delta_err']:.4f} +- {v['delta_err_sd']:.4f} um | "
      f"RMS bias {np.sqrt(np.mean(d ** 2)):.4f} um | model time {v['t_pipeline']:.2f} s")

# ---- other checks quoted in the text
for f in (PROC / "edge_verification.txt", PROC / "radial_shortcut.json", RAW / "resolution_condition.csv"):
    if f.exists(): P(f"--- {f.name}\n" + f.read_text().strip())

(RES / "paper_numbers.json").write_text(json.dumps(N, indent=1))
(RES / "paper_numbers.txt").write_text("\n".join(L) + "\n")
