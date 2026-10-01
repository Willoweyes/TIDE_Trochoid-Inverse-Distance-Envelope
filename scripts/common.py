"""Shared helpers: repository root, executables (Linux/macOS/Windows), temp dir."""
import os, sys, json, subprocess, tempfile, math
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / "benchmark" / "test_cases"
TMP = Path(tempfile.mkdtemp(prefix="tide_"))
ENV = dict(os.environ, OMP_NUM_THREADS="1")          # all timings single-thread
sys.path.insert(0, str(ROOT / "src" / "reference_python"))

def exe(name):
    for p in (ROOT/"build"/name, ROOT/"build"/f"{name}.exe",
              ROOT/"build"/"Release"/f"{name}.exe", ROOT/"build"/"Release"/name):
        if p.exists(): return str(p)
    sys.exit(f"{name} not found - build first (see README)")

def case_cfg(case, **mod):
    """Load benchmark/test_cases/paper_<case>.json, override keys, write to TMP."""
    c = json.loads((CASES / f"paper_{case}.json").read_text())
    for k, v in mod.items():
        for sec in ("input", "fixed", "grid"):
            if k in c[sec]: c[sec][k] = v
    f = TMP / (f"{case}_" + "_".join(f"{k}{v}" for k, v in sorted(mod.items())) + ".json")
    f.write_text(json.dumps(c)); return str(f), c

def edge_speed(c):
    """v_e = omega D / 2 in mm/s."""
    Td = c["fixed"]["Td"]; n = c["input"]["vc"] / (Td * 1e-3 * math.pi)
    return n * math.pi / 30 * Td / 2

REPS = int(os.environ.get("TIDE_REPS", "1"))       # timing protocol: median of REPS runs

def run(mode, cfgf, gb, dt=None, dump=None, edge=None, reps=None):
    """Run tide_tools once per repetition; results are deterministic, the
    reported wall-clock t_s is the median over the repetitions."""
    cmd = [exe("tide_tools"), mode, "--config", cfgf, "--grid-base", str(gb)]
    if dt is not None: cmd += ["--dt", repr(dt)]
    if dump: cmd += ["--dump-bin", str(dump)]
    if edge: cmd += ["--edge", str(edge)]
    out, ts = None, []
    for _ in range(reps or REPS):
        out = json.loads(subprocess.run(cmd, capture_output=True, text=True, env=ENV, check=True).stdout)
        ts.append(out["t_s"])
    ts.sort(); out["t_s"] = ts[len(ts) // 2]; out["t_runs"] = len(ts)
    return out

def points(cfgf, Px, Py):
    import numpy as np
    fi, fo = TMP/"p.bin", TMP/"o.bin"
    np.column_stack([Px, Py]).astype(np.float64).tofile(fi)
    subprocess.run([exe("tide_tools"), "points", "--config", cfgf, "--in", str(fi), "--out", str(fo)],
                   check=True, env=ENV)
    return np.fromfile(fo)
