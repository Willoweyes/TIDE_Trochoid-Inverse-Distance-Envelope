"""Hybrid-pipeline study (design A).
Part 1  cloud : FSM Sa/Sz over its numerical settings (grid, dt ratio r, edge points)
                vs TIDE over sampling only, cases A1-A3.
Part 2  ref   : independent reference = FSM at 800 and 1600 cells (r = 0.25, 8n edge points),
                Richardson-extrapolated to zero cell size,
                on 36 conditions around A1 (cross-checked against TIDE).
Part 3  cores : every core option on the 36 conditions (Sa + wall-clock).
Part 4  learn : y = S_ref + delta_true + eps; GP on y - M; delta error vs pipeline time.
Usage: python pipeline_study.py [cloud] [ref] [cores] [learn]"""
import sys, csv, json, time, itertools, os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import ROOT, case_cfg, edge_speed, run
RAW, PROC = ROOT/"data/raw", ROOT/"data/processed"
steps = sys.argv[1:] or ["cloud", "ref", "cores", "learn"]
def dt_for(c, gb, r): return r*(c["grid"]["length"]/gb)/edge_speed(c)
def save(p, rows):
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

# ---------------- Part 1 ----------------
if "cloud" in steps:
    partc = RAW/"fsm_settings_cloud_partial.csv"
    rows = list(csv.DictReader(open(partc))) if partc.exists() else []
    done_c = {r["case"] for r in rows}
    for case in ["A1", "A2", "A3"]:
        if case in done_c: continue
        f, c = case_cfg(case)
        for gb in [50, 100, 200, 400, 1600]:
            t = run("tide", f, gb); rows.append(dict(case=case, method="tide", grid=gb, r="", edge_mult="", Sa=t["Sa"], Sz=t["Sz"], t_s=t["t_s"]))
        for gb, r, em in itertools.product([50, 100, 200, 400], [0.1, 0.25, 0.5, 0.9], [1, 2, 4, 8]):
            x = run("fsm", f, gb, dt=dt_for(c, gb, r), edge=em*gb)
            rows.append(dict(case=case, method="fsm", grid=gb, r=r, edge_mult=em, Sa=x["Sa"], Sz=x["Sz"], t_s=x["t_s"]))
            print("cloud", case, gb, r, em, round(x["Sa"], 4), round(x["t_s"], 2), flush=True)
        save(partc, rows)
    save(RAW/"fsm_settings_cloud.csv", rows); partc.unlink()

# ---------------- 36 conditions ----------------
FZ = [.15, .25, .35, .45, .55, .65]; VC = [120, 195, 270]; AP = [.3, .7]
COND = list(itertools.product(FZ, VC, AP))

PREV = {}
if (RAW/"pipeline_reference_1600.csv").exists():   # cached 1600-cell runs
    PREV = {(float(r["fz"]), float(r["vc"]), float(r["ap"])): r for r in csv.DictReader(open(RAW/"pipeline_reference_1600.csv"))}
def ref_one(cond):
    fz, vc, ap = cond; f, c = case_cfg("A1", fz=fz, vc=vc, ap=ap)
    p = PREV.get((fz, vc, ap))
    if p: s16, t16, tt = float(p.get("Sa_fsm1600", p.get("Sa_ref"))), float(p["t_ref"]), float(p["Sa_tide1600"])
    else:
        x = run("fsm", f, 1600, dt=dt_for(c, 1600, 0.25), edge=8*1600); s16, t16 = x["Sa"], x["t_s"]
        tt = run("tide", f, 1600)["Sa"]
    h = run("fsm", f, 800, dt=dt_for(c, 800, 0.25), edge=8*800)
    # FSM converges at first order in the cell size (Sec. 5): Richardson extrapolation to zero cell size
    return dict(fz=fz, vc=vc, ap=ap, Sa_fsm1600=s16, Sa_fsm800=h["Sa"], Sa_ref=2*s16 - h["Sa"],
                Sa_tide1600=tt, t_ref=t16 + h["t_s"])
if "ref" in steps and __name__ == "__main__":      # guard: Windows starts workers by spawn and re-imports this file
    with ProcessPoolExecutor(int(os.environ.get("TIDE_REF_WORKERS", "2"))) as ex: rows = list(ex.map(ref_one, COND))   # results are deterministic; t_ref is not quoted
    save(RAW/"pipeline_reference.csv", rows)
    for k in ("Sa_fsm1600", "Sa_ref"):
        d = np.array([r[k]/r["Sa_tide1600"]-1 for r in rows])*100
        print(k, "vs TIDE1600: %.3f .. %.3f %%" % (d.min(), d.max()))

CORES = [("TIDE-25", "tide", 25, None, None), ("TIDE-50", "tide", 50, None, None),
         ("TIDE-100", "tide", 100, None, None), ("TIDE-200", "tide", 200, None, None)]
for gb in [50, 100, 200, 400]:
    CORES.append((f"FSM-{gb}", "fsm", gb, 0.5, 4))
CORES += [("FSM-200 released", "fsm", 200, "rel", 4), ("FSM-100 1n", "fsm", 100, 0.5, 1),
          ("FSM-200 r0.9", "fsm", 200, 0.9, 4), ("FSM-400 8n r0.25", "fsm", 400, 0.25, 8)]
if "cores" in steps:
    partk = RAW/"pipeline_cores_partial.csv"
    rows = list(csv.DictReader(open(partk))) if partk.exists() else []
    done_k = {r["core"] for r in rows}
    for name, meth, gb, r, em in CORES:
        if name in done_k: continue
        for fz, vc, ap in COND:
            f, c = case_cfg("A1", fz=fz, vc=vc, ap=ap)
            if meth == "tide": x = run("tide", f, gb)
            else:
                dt = 8e-7 if r == "rel" else dt_for(c, gb, r)
                x = run("fsm", f, gb, dt=dt, edge=em*gb)
            rows.append(dict(core=name, fz=fz, vc=vc, ap=ap, Sa=x["Sa"], t_s=x["t_s"]))
        print("core", name, "total t %.2f s" % sum(float(r_["t_s"]) for r_ in rows if r_["core"] == name), flush=True)
        save(partk, rows)
    save(RAW/"pipeline_cores.csv", rows); partk.unlink()

if "learn" in steps:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel as C, RBF, WhiteKernel
    # reference = TIDE at 1600 cells, verified against a brute-force search (Sect. 5.1)
    ref = {(float(r["fz"]), float(r["vc"]), float(r["ap"])): float(r["Sa_tide1600"]) for r in csv.DictReader(open(RAW/"pipeline_reference.csv"))}
    CR = list(csv.DictReader(open(RAW/"pipeline_cores.csv")))
    X = np.array(COND, float); S = np.array([ref[tuple(x)] for x in X])
    K_SA, R_INS, H_MIN = 8/(18*np.sqrt(3)), 5.0, 0.004
    D = K_SA*(H_MIN/2)*(1 + R_INS*H_MIN/X[:, 0]**2)*1000
    Xs = (X - X.mean(0))/X.std(0); out = {}
    import warnings; warnings.filterwarnings("ignore")
    for name, *_ in CORES:
        rr = {(float(r["fz"]), float(r["vc"]), float(r["ap"])): r for r in CR if r["core"] == name}
        M = np.array([float(rr[tuple(x)]["Sa"]) for x in X]); tcore = sum(float(rr[tuple(x)]["t_s"]) for x in X)
        rng = np.random.default_rng(0); rec, pred, tgp = [], [], []
        for rep in range(20):
            y = S + D + rng.normal(0, 0.02, len(D)); tr = np.concatenate([rng.permutation(np.where(X[:, 0] == fz)[0])[:3] for fz in FZ])   # 3 of 6 per feed
            te = np.setdiff1d(np.arange(len(D)), tr)
            t0 = time.perf_counter()
            gp = GaussianProcessRegressor(C(0.1)*RBF([1., 1., 1.]) + WhiteKernel(1e-3), normalize_y=True,
                                          n_restarts_optimizer=3, random_state=rep).fit(Xs[tr], (y - M)[tr])
            dh = gp.predict(Xs[te]); tgp.append(time.perf_counter() - t0)
            rec.append(np.sqrt(np.mean((dh - D[te])**2))); pred.append(np.sqrt(np.mean((M[te] + dh - (S + D)[te])**2)))
        bias = (M/S - 1)*100
        out[name] = dict(delta_err=float(np.mean(rec)), delta_err_sd=float(np.std(rec)), pred_err=float(np.mean(pred)),
                         t_core_36=tcore, t_gp=float(np.mean(tgp)), t_pipeline=tcore + float(np.mean(tgp)),
                         bias_min=float(bias.min()), bias_max=float(bias.max()))
        print(name, {k: round(v, 4) for k, v in out[name].items()}, flush=True)
    (PROC/"pipeline_study.json").write_text(json.dumps(out, indent=1))
