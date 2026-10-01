"""Newton steps needed on the edge-ray crossing f(t) = 0 (third row of Table 2).

For random points in the window of each case, the seeds (six per insert, after
realignment) are refined with 0..4 Newton steps.  For the engaged candidates
(inside the window, ray pointing at P, |l| <= DeltaL) the script reports
  miss  = largest distance |f(t)| between the edge ray and P after s steps [um],
  dz    = largest difference of the TIDE height from the four-step height [um].
Writes data/processed/newton_steps.txt.  Usage: python newton_steps.py [n_points_per_case]"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import ROOT, case_cfg
import tide as T

N = int(sys.argv[1]) if len(sys.argv) > 1 else 120
rng = np.random.default_rng(11)
lines = []
for case in ["A1", "A2", "A3"]:
    _, cfg = case_cfg(case); k = T.kinematics_from_config(cfg); dL = T.engaged_half_width(k)
    Px, Py = rng.uniform(0, k.Lx, N), rng.uniform(0, k.Ly, N)
    miss = {s: 0.0 for s in range(5)}; ncand = 0
    for K in range(k.z_n):
        rho = k.D / 2 + K * k.eps_r; phi0 = k.phi + 2 * np.pi * K / k.z_n; xc = k.x0
        t0 = T._realign(T._candidate_times(Px, Py, k, rho, phi0, xc), Px, Py, k, rho, phi0, xc, gf=k.gamma_f)
        def f_r(t):
            a = phi0 - k.omega * t - k.gamma_f; ca, sa = np.cos(a), np.sin(a)
            dx, dy = Px - xc, Py - (k.y0 + k.v_f * t)
            return ca * dy + sa * dx, ca * dx - sa * dy, k.omega * sa * dy - ca * k.v_f - k.omega * ca * dx
        ts = [t0]
        for _ in range(4):
            t = ts[-1]; f, _r, fp = f_r(t); ts.append(t - f / fp)
        _f, r4, _ = f_r(ts[4]); l4 = r4 - rho
        ok = (ts[4] >= 0) & (ts[4] <= k.t_total) & (r4 > 0) & (np.abs(l4) <= dL)
        ncand += int(ok.sum())
        for s in range(5):
            miss[s] = max(miss[s], float(np.max(np.abs(f_r(ts[s])[0])[ok])) * 1e3)
    zeta = lambda l: (k.R - np.sqrt(k.R ** 2 - np.minimum(np.abs(l), k.R) ** 2)) * np.cos(k.gamma_p)
    def height(steps):
        z = np.full(Px.shape, np.inf)
        for K in range(k.z_n):
            l = T.contact_offset(Px, Py, k, K, cross_steps=steps)
            z = np.minimum(z, np.where(np.abs(l) <= dL, zeta(l) + K * k.eps_a, np.inf).min(axis=0))
        z[np.isinf(z)] = k.a_p; return z * 1e3
    z4 = height(4)
    dz = {s: float(np.max(np.abs(height(s) - z4))) for s in range(5)}
    lines.append(f"{case} n={N} engaged candidates={ncand} | " + " | ".join(
        f"{s} steps: miss {miss[s]:.1e} um, dz {dz[s]:.1e} um" for s in range(5)))
    print(lines[-1], flush=True)
(ROOT / "data/processed/newton_steps.txt").write_text("\n".join(lines) + "\n")
