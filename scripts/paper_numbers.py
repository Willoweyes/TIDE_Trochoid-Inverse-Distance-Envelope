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
         grid_F1=reach["fsm"][0], ratio=round(ratio, 1), ratio_floor=int(ratio // (50 if ratio >= 150 else 10) * (50 if ratio >= 150 else 10)))
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

# ---- Sect. 5.2: runtime spread over v_c and f_z (200 cells)
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

# ---- Sect. 5.2: TIDE runtime over the sampling ladder (maxima over A1-A3, as quoted)
TS = list(csv.DictReader(open(RAW / "tide_sampling.csv")))
g_ = lambda n: [r for r in TS if int(r["grid_base"]) == n]
us_ = [1e6 * float(r["t_s"]) / int(r["n_points"]) for r in TS]
P(f"TIDE runtime: max over cases {max(float(r['t_s']) for r in g_(25)):.4f} s at 25 cells ({g_(25)[0]['n_points']} points), "
  f"{max(float(r['t_s']) for r in g_(1600)):.1f} s at 1600 cells; {min(us_):.1f}-{max(us_):.1f} us per point")

# ---- Sect. 5.3: comparison with the released FSM code (compare_bakhshan/results)
CMP = ROOT / "compare_bakhshan" / "results"
if (CMP / "exp3b_16threads_allcases.csv").exists():
    A_ = list(csv.DictReader(open(CMP / "exp1_accuracy.csv")))
    P(f"released code vs this FSM: {len({(r['regime'], r['case'], r['grid']) for r in A_})} combinations, "
      f"largest height difference {max(float(r['max_abs_dz_um']) for r in A_):.1e} um")
    T_ = {(r["regime"], r["case"], r["grid"], r["impl"]): float(r["t_s"]) for r in csv.DictReader(open(CMP / "exp2_timing_1thread.csv"))}
    sl = [T_[("matched", c, str(g), "gcc_native")] / T_[("matched", c, str(g), "ours")] for c in C for g in G]
    P(f"released code (GCC build) / this FSM, one thread: {min(sl):.1f}-{max(sl):.1f}")
    M_ = list(csv.DictReader(open(CMP / "exp3b_16threads_allcases.csv")))
    t100 = np.mean([T_[("matched", c, "100", "tide")] for c in C])
    for impl in ("gcc_native", "pyd"):
        v = [float(r["t_s"]) for r in M_ if r["impl"] == impl and r["grid"] == "1600"]
        P(f"released FSM {impl}, 16 threads, 1600 cells: mean {np.mean(v):.2f} s (per case {[round(x, 2) for x in v]}); "
          f"TIDE one thread 100 cells {t100:.4f} s; ratio of means {np.mean(v) / t100:.1f}, smallest per case {min(v) / t100:.1f}")

# ---- Table 2, third row: Newton steps (data/processed/newton_steps.txt)
import re
if (PROC / "newton_steps.txt").exists():
    mx = {}
    for ln in (PROC / "newton_steps.txt").read_text().splitlines():
        for s, miss, dz in re.findall(r"(\d) steps: miss ([\d.e+-]+) um, dz ([\d.e+-]+) um", ln):
            a = mx.setdefault(int(s), [0.0, 0.0]); a[0] = max(a[0], float(miss)); a[1] = max(a[1], float(dz))
    P("Newton steps, largest miss of P / height change [um]: " + ", ".join(f"{s}: {a[0]:.1e}/{a[1]:.1e}" for s, a in sorted(mx.items())))

# ---- Sect. 6: known discrepancy over the 36 conditions
dtrue = [8 / (18 * np.sqrt(3)) * (0.004 / 2) * (1 + 5.0 * 0.004 / f ** 2) * 1000 for f in (.15, .25, .35, .45, .55, .65)]
P(f"delta_true for h_min = 4 um: {min(dtrue):.2f} .. {max(dtrue):.2f} um")

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
