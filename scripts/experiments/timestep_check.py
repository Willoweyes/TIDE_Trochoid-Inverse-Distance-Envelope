"""Effect of the time step at a fixed grid (100 cells): S_a of the forward
method with dt tied to the grid (r = 0.5) and with a sixteen-fold smaller dt.
Writes data/processed/timestep_check.json"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import ROOT, case_cfg, edge_speed, run
out = {}
for case in ["A1", "A2", "A3"]:
    f, c = case_cfg(case)
    dt = 0.5 * (c["grid"]["length"] / 100) / edge_speed(c)
    a = run("fsm", f, 100, dt=dt)["Sa"]; b = run("fsm", f, 100, dt=dt / 16)["Sa"]
    t = run("tide", f, 100)["Sa"]
    out[case] = dict(Sa_r05=a, Sa_r05_over16=b, change_pct=100 * (b / a - 1), Sa_tide100=t)
    print(case, out[case], flush=True)
(ROOT / "data/processed/timestep_check.json").write_text(json.dumps(out, indent=1))
