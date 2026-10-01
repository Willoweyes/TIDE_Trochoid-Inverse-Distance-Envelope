"""Contact along the edge ray (TIDE) versus the radial shortcut.

The shortcut takes the minimum distance d* between P and the tip trajectory
and converts it into an edge offset as if the contact lay along the tool
radius (lateral_offset).  TIDE solves the edge-ray crossing f(t) = 0 instead.
Both heights are computed on the full grid of A1-A3 (200 cells) and compared.
Writes data/processed/radial_shortcut.json"""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import ROOT, case_cfg
import tide as T
from roughness import compute_roughness

GB = int(sys.argv[1]) if len(sys.argv) > 1 else 200
out = {}
for case in ["A1", "A2", "A3"]:
    _, cfg = case_cfg(case, grid_base=GB)
    k = T.kinematics_from_config(cfg)
    X, Y = T.workpiece_grid(k); Px, Py = X.ravel(), Y.ravel()
    dL = T.engaged_half_width(k); cgp = np.cos(k.gamma_p)
    z_ray = np.full(Px.shape, np.inf); z_rad = np.full(Px.shape, np.inf)
    for K in range(k.z_n):
        l = T.contact_offset(Px, Py, k, K)                        # (6, N), signed
        al = np.abs(l)
        z_ray = np.minimum(z_ray, np.where(al <= dL, T.zeta(np.minimum(al, k.R), k.R, cgp) + K*k.eps_a, np.inf).min(axis=0))
        rho = k.D/2 + K*k.eps_r
        d = T.dstar_insert(Px, Py, k, K, newton_steps=2)           # nearest tip distance
        lr = T.lateral_offset(d, rho, k.R, k.gamma_f, k.gamma_p)  # radial conversion
        z_rad = np.minimum(z_rad, np.where(lr <= dL, T.zeta(np.minimum(lr, k.R), k.R, cgp) + K*k.eps_a, np.inf))
    for z in (z_ray, z_rad): z[np.isinf(z)] = k.a_p
    z_ray = (k.z0 + z_ray)*1000; z_rad = (k.z0 + z_rad)*1000       # um
    dz = np.abs(z_rad - z_ray)
    ra, rr = compute_roughness(z_ray.reshape(X.shape)), compute_roughness(z_rad.reshape(X.shape))
    out[case] = dict(grid=GB, n=int(dz.size), max_um=float(dz.max()), p999_um=float(np.percentile(dz, 99.9)),
                     p99_um=float(np.percentile(dz, 99)), mean_um=float(dz.mean()),
                     Sa_ray=float(ra["Sa"]), Sa_radial=float(rr["Sa"]), dSa_pct=float(100*(rr["Sa"]/ra["Sa"]-1)),
                     Sz_ray=float(ra["Sz"]), Sz_radial=float(rr["Sz"]), dSz_pct=float(100*(rr["Sz"]/ra["Sz"]-1)))
    print(case, {kk: (round(v, 5) if isinstance(v, float) else v) for kk, v in out[case].items()}, flush=True)
(ROOT/"data/processed/radial_shortcut.json").write_text(json.dumps(out, indent=1))
