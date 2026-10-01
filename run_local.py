"""Full local reproduction (Step B) on this Windows machine, with logging and resume.

Same steps and protocol as code/run_all.bat (single thread, median of three), with the
MinGW-w64 GCC 13.3 toolchain instead of Visual Studio (any C++17 toolchain works; build first, see README).  Phase 1 runs the steps whose
results carry no timing (verify, surfaces, ref, radial_shortcut, timestep_check) in
parallel; phase 2 runs every timed step one at a time so that they do not disturb each
other; phase 3 redraws the figures and collects the numbers.

Usage:  python run_local.py            (resumes: steps already in results/logs/state.json are skipped)
        python run_local.py --fresh    (ignore state.json)
"""
import os, sys, json, time, subprocess, datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
CODE = HERE
LOGS = CODE / "results" / "logs"; LOGS.mkdir(parents=True, exist_ok=True)
STATE = LOGS / "state.json"
GCC_BIN = os.environ.get("GCC_BIN", "")      # optional: folder of g++ (e.g. MinGW-w64 bin) to put on PATH
PY = sys.executable

ENV = dict(os.environ, PATH=(GCC_BIN + os.pathsep if GCC_BIN else "") + os.environ["PATH"], OMP_NUM_THREADS="1",
           OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONUTF8="1", MPLBACKEND="Agg")
state = {} if "--fresh" in sys.argv or not STATE.exists() else json.loads(STATE.read_text())
def now(): return datetime.datetime.now().isoformat(timespec="seconds")
def save(): STATE.write_text(json.dumps(state, indent=1))

def launch(name, args, reps, extra=None):
    env = dict(ENV, TIDE_REPS=str(reps), **(extra or {}))
    log = open(LOGS / f"{name}.log", "w", encoding="utf-8")
    p = subprocess.Popen([PY, *args], cwd=CODE, env=env, stdout=log, stderr=subprocess.STDOUT)
    return dict(name=name, p=p, t0=time.time(), start=now(), log=log)

def finish(h):
    rc = h["p"].wait(); h["log"].close()
    state[h["name"]] = dict(rc=rc, start=h["start"], end=now(), seconds=round(time.time() - h["t0"], 1))
    save(); print(f"[{now()}] {h['name']}: rc={rc} {state[h['name']]['seconds']} s", flush=True)
    return rc

def todo(name): return not (name in state and state[name]["rc"] == 0)

RE = "scripts/run_experiments.py"; PS = "scripts/experiments/pipeline_study.py"
print(f"[{now()}] start; python {PY}", flush=True)

if todo("machine_info"):
    finish(launch("machine_info", ["scripts/machine_info.py"], 1))

# ---- phase 1: no timing quoted -> run side by side
p1 = [("verify", [RE, "verify"], 1, None), ("surfaces", [RE, "surfaces"], 1, None),
      ("ref", [PS, "ref"], 1, {"TIDE_REF_WORKERS": "8"}),
      ("radial_shortcut", ["scripts/experiments/radial_shortcut.py"], 1, None),
      ("timestep_check", ["scripts/experiments/timestep_check.py"], 1, None)]
hs = [launch(*a) for a in p1 if todo(a[0])]
bad = [h["name"] for h in hs if finish(h) != 0]
if bad: sys.exit(f"phase 1 failed: {bad}")

# ---- phase 2: timed steps, one at a time, median of three
p2 = [("sampling", [RE, "sampling"]), ("ladder", [RE, "ladder"]), ("cost", [RE, "cost"]),
      ("resolution", [RE, "resolution"]), ("scaling", [RE, "scaling"]),
      ("cloud", [PS, "cloud"]), ("cores", [PS, "cores"]), ("learn", [PS, "learn"])]
for name, args in p2:
    if todo(name) and finish(launch(name, args, 3)) != 0: sys.exit(f"{name} failed")

# ---- phase 3: figures and numbers
for f in ("fig_schematic", "fig_convergence", "fig_maps", "fig_cost", "fig_pipeline"):
    if todo(f) and finish(launch(f, [f"scripts/plotting/{f}.py"], 1)) != 0: sys.exit(f"{f} failed")
if todo("paper_numbers") and finish(launch("paper_numbers", ["scripts/paper_numbers.py"], 1)) != 0:
    sys.exit("paper_numbers failed")
print(f"[{now()}] ALL DONE", flush=True)
