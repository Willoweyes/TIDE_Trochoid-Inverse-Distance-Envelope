"""Compare the ORIGINAL released FSM of Bakhshan et al. (HadiBakhshan/surf-topo, commit 1f5d3e0) with the C++ FSM of this bundle.

Implementations
  pyd          released prebuilt surftopo.cp311-win_amd64.pyd (MSVC, /fp:fast /arch:AVX2), driven as in data_set_builder.py
  gcc_native   the UNMODIFIED simulation.cpp built with GCC 13.3 and the GCC flags of its CMakeLists (-O3 -march=native -fopenmp)
  gcc_generic  same source, -O3 -fopenmp (baseline x86-64)
  ours         code/include/fsm_core.hpp (the C++ FSM used in the paper), single thread

Experiments (results in compare_bakhshan/results/)
  1  accuracy   node-by-node comparison with `ours` on the doubling ladder, two time-step regimes (matched, released)
  2  timing     single thread, median of three, matched-dt ladder for A1-A3 and the released setting
  3  threads    wall-clock against OMP threads (pyd, gcc_native for FSM; TIDE) for A1 at the 1 % points of the paper
Usage: python run_compare.py [acc time threads]      (resumes: rows already in the CSVs are skipped)
"""
import os, sys, csv, json, subprocess, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parent
CODE = BUNDLE
RES = Path(os.environ.get("COMPARE_RESULTS", HERE / "results")); RES.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(CODE / "scripts"))
from common import case_cfg, edge_speed, TMP          # TMP: private temp dir (set TEMP to a folder with space)

GCC_BIN = os.environ.get("GCC_BIN", "")                 # folder of g++ (needed on PATH for the OpenMP runtime on Windows)
SFX = ".exe" if os.name == "nt" else ""
PYD_PY = os.environ.get("PYD_PYTHON", sys.executable)    # python that can import surftopo (cp311 for the prebuilt pyd)
from common import exe
EXE = {"ours": exe("tide_tools"), "gcc_native": str(HERE / ("bakhshan_gcc_native" + SFX)),
       "gcc_generic": str(HERE / ("bakhshan_gcc_generic" + SFX))}
CASES = ["A1", "A2", "A3"]
GRIDS = [int(g) for g in os.environ.get("COMPARE_GRIDS", "25,50,100,200,400,800,1600").split(",")]   # e.g. COMPARE_GRIDS=25,50 for a quick check
def _have_pyd():
    r = subprocess.run([PYD_PY, "-c", "import sys; sys.path.insert(0, sys.argv[1]); import surftopo", os.environ.get("SURFTOPO_BIN", "surf-topo/bin")], capture_output=True)
    return r.returncode == 0
HAVE_PYD = _have_pyd()
if not HAVE_PYD: print("note: surftopo (prebuilt Python module) not importable; the 'pyd' rows are skipped", flush=True)
steps = sys.argv[1:] or ["acc", "time", "threads", "mt3"]

def cmd(impl, cfgf, gb, dt, dump=None):
    if impl == "ours":
        c = [EXE["ours"], "fsm", "--config", cfgf, "--grid-base", str(gb), "--dt", repr(dt)]
    elif impl == "tide":
        c = [EXE["ours"], "tide", "--config", cfgf, "--grid-base", str(gb)]
    elif impl == "pyd":
        c = [PYD_PY, str(HERE / "run_pyd.py"), "--config", cfgf, "--grid-base", str(gb), "--dt", repr(dt)]
    else:
        c = [EXE[impl], "--config", cfgf, "--grid-base", str(gb), "--dt", repr(dt)]
    if dump: c += ["--dump-bin", str(dump)]
    return c

def run(impl, cfgf, gb, dt, threads=1, reps=1, dump=None):
    """median wall-clock of `reps` runs; the first run writes `dump` (results are deterministic)"""
    env = dict(os.environ, PATH=(GCC_BIN + os.pathsep if GCC_BIN else "") + os.environ["PATH"], OMP_NUM_THREADS=str(threads), PYTHONUTF8="1")
    ts, out = [], None
    for k in range(reps):
        o = json.loads(subprocess.run(cmd(impl, cfgf, gb, dt, dump if k == 0 else None), capture_output=True, text=True, env=env, check=True).stdout)
        ts.append(o["t_s"]); out = out or o
    ts.sort(); out["t_s"] = ts[len(ts) // 2]; out["t_all"] = ts
    return out

def matched_dt(c, gb): return (c["grid"]["length"] / gb) / (2 * edge_speed(c))
def regime_dt(c, gb, regime): return matched_dt(c, gb) if regime == "matched" else 8e-7

def load(path): return list(csv.DictReader(open(path))) if Path(path).exists() else []
def append(path, row):
    new = not Path(path).exists()
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row));
        if new: w.writeheader()
        w.writerow(row)
def done(rows, **kv): return any(all(str(r[k]) == str(v) for k, v in kv.items()) for r in rows)

# ------------------------------------------------------------------ 1  accuracy, node by node
if "acc" in steps:
    P = RES / "exp1_accuracy.csv"; rows = load(P)
    for regime in ("matched", "released"):
        for case in CASES:
            cfgf, c = case_cfg(case)
            for gb in GRIDS:
                dt = regime_dt(c, gb, regime)
                if done(rows, regime=regime, case=case, grid=gb): continue
                ref = TMP / "ref.bin"; ro = run("ours", cfgf, gb, dt, dump=ref); zr = np.fromfile(ref)
                for impl in [i for i in ("pyd", "gcc_native", "gcc_generic") if i != "pyd" or HAVE_PYD]:
                    f = TMP / f"{impl}.bin"; o = run(impl, cfgf, gb, dt, threads=8, dump=f); z = np.fromfile(f)
                    d = z - zr
                    append(P, dict(regime=regime, case=case, grid=gb, dt=dt, impl=impl, n_points=z.size,
                                   Sa_ours=ro["Sa"], Sa_impl=o["Sa"], Sq_ours=ro["Sq"], Sq_impl=o["Sq"], Sz_ours=ro["Sz"], Sz_impl=o["Sz"],
                                   rel_Sa=o["Sa"] / ro["Sa"] - 1, rel_Sz=o["Sz"] / ro["Sz"] - 1,
                                   max_abs_dz_um=float(np.abs(d).max()), n_nodes_dz_gt_1e9=int((np.abs(d) > 1e-9).sum()),
                                   n_nodes_dz_gt_1e3=int((np.abs(d) > 1e-3).sum())))
                    print("acc", regime, case, gb, impl, "max|dz| %.2e um, nodes>1e-9: %d" % (np.abs(d).max(), (np.abs(d) > 1e-9).sum()), flush=True)
                rows = load(P)
    # thread determinism of the released implementations (bitwise)
    P2 = RES / "exp1b_thread_determinism.csv"
    if not Path(P2).exists():
        cfgf, c = case_cfg("A1"); gb = 400; dt = matched_dt(c, gb)
        for impl in [i for i in ("pyd", "gcc_native") if i != "pyd" or HAVE_PYD]:
            zs = {}
            for th in (1, 16):
                f = TMP / f"{impl}_{th}.bin"; run(impl, cfgf, gb, dt, threads=th, dump=f); zs[th] = np.fromfile(f)
            append(P2, dict(impl=impl, case="A1", grid=gb, threads_a=1, threads_b=16, bitwise_identical=bool(np.array_equal(zs[1], zs[16])),
                            max_abs_diff_um=float(np.abs(zs[1] - zs[16]).max())))
            print("determinism", impl, np.array_equal(zs[1], zs[16]), flush=True)

# ------------------------------------------------------------------ 2  single-thread timing, median of three
if "time" in steps:
    P = RES / "exp2_timing_1thread.csv"; rows = load(P)
    jobs = [("matched", case, gb) for case in CASES for gb in GRIDS] + [("released", case, 200) for case in CASES]
    for regime, case, gb in jobs:
        cfgf, c = case_cfg(case); dt = regime_dt(c, gb, regime)
        for impl in [i for i in ("ours", "gcc_native", "pyd") if i != "pyd" or HAVE_PYD]:
            if done(rows, regime=regime, case=case, grid=gb, impl=impl): continue
            o = run(impl, cfgf, gb, dt, threads=1, reps=3)
            append(P, dict(regime=regime, case=case, grid=gb, dt=dt, impl=impl, threads=1, Sa=o["Sa"], Sz=o["Sz"], t_s=o["t_s"],
                           t_runs=";".join(f"{t:.4f}" for t in o["t_all"])))
            print("time", regime, case, gb, impl, "%.4f s" % o["t_s"], flush=True)
        # TIDE at the same grid, single thread, for the ratio
        if regime == "matched" and not done(rows, regime=regime, case=case, grid=gb, impl="tide"):
            o = run("tide", cfgf, gb, dt, threads=1, reps=3)
            append(P, dict(regime=regime, case=case, grid=gb, dt=dt, impl="tide", threads=1, Sa=o["Sa"], Sz=o["Sz"], t_s=o["t_s"],
                           t_runs=";".join(f"{t:.4f}" for t in o["t_all"])))
            print("time", regime, case, gb, "tide", "%.4f s" % o["t_s"], flush=True)

# ------------------------------------------------------------------ 3  wall-clock against thread count (A1, the 1 % points)
if "threads" in steps:
    P = RES / "exp3_threads.csv"; rows = load(P)
    cfgf, c = case_cfg("A1")
    jobs = [j for j in [("pyd", 1600), ("gcc_native", 1600), ("tide", 100), ("tide", 1600)] if j[0] != "pyd" or HAVE_PYD]
    for impl, gb in jobs:
        dt = matched_dt(c, gb)
        for th in (1, 2, 4, 8, 16):
            if done(rows, impl=impl, grid=gb, threads=th): continue
            o = run(impl, cfgf, gb, dt, threads=th, reps=3)
            append(P, dict(case="A1", impl=impl, grid=gb, threads=th, Sa=o["Sa"], Sz=o["Sz"], t_s=o["t_s"], t_runs=";".join(f"{t:.4f}" for t in o["t_all"])))
            print("threads", impl, gb, th, "%.4f s" % o["t_s"], flush=True)
# ------------------------------------------------------------------ 3b  16 threads, all three cases (released FSM at 1600 cells, TIDE at 100 and 1600 cells)
if "mt3" in steps:
    P = RES / "exp3b_16threads_allcases.csv"; rows = load(P)
    for case in CASES:
        cfgf, c = case_cfg(case)
        for impl, gb in [j for j in (("pyd", 1600), ("gcc_native", 1600), ("tide", 100), ("tide", 1600)) if j[0] != "pyd" or HAVE_PYD]:
            for th in (16,):
                if done(rows, case=case, impl=impl, grid=gb, threads=th): continue
                o = run(impl, cfgf, gb, matched_dt(c, gb), threads=th, reps=3)
                append(P, dict(case=case, impl=impl, grid=gb, threads=th, Sa=o["Sa"], Sz=o["Sz"], t_s=o["t_s"], t_runs=";".join(f"{t:.4f}" for t in o["t_all"])))
                print("mt3", case, impl, gb, th, "%.4f s" % o["t_s"], flush=True)
print("DONE", steps, flush=True)
