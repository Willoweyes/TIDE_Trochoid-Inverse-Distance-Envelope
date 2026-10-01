"""Run the RELEASED, prebuilt surftopo.pyd (MSVC build, bin/surftopo.cp311-win_amd64.pyd of HadiBakhshan/surf-topo)
exactly as scripts/data_set_builder.py does, for one case file.
  python run_pyd.py --config F --grid-base N --dt DT [--edge NE] [--dump-bin FILE]
Prints one JSON line with Sa, Sq, Sz (um; estimators of tide_core.hpp) and t_s = wall clock of run_surface_simulation."""
import sys, os, json, math, time, argparse
import numpy as np
sys.path.insert(0, os.environ.get("SURFTOPO_BIN", "surf-topo/bin"))   # folder with surftopo.*.pyd / .so
import surftopo
ap = argparse.ArgumentParser()
ap.add_argument("--config"); ap.add_argument("--grid-base", type=int, default=-1); ap.add_argument("--dt", type=float, default=8e-7)
ap.add_argument("--edge", type=int, default=0); ap.add_argument("--dump-bin", default="")
a = ap.parse_args()
c = json.load(open(a.config)); i, f, g = c["input"], c["fixed"], c["grid"]
gb = a.grid_base if a.grid_base > 0 else g["grid_base"]
Td, vc, fz, z_n = f["Td"], i["vc"], i["fz"], i["z"]
ws = vc / (Td * 0.001 * math.pi); vf = fz * z_n * ws / 60; omega = ws * math.pi / 30
m, n = gb * 2, gb; Lx, Ly = g["length"] * 2, g["length"]; grid_t = max(m, n) * 2; total = (1.5 * Td + Ly) / vf
tool = surftopo.Tool(gama_f=f["gama_f"], gama_p=f["gama_p"], D=Td, eps_r=i["eps_r"], eps_a=i["eps_a"], phi=f["phi"],
                     omega=omega, v_f=vf, x_0=f["x_0"], y_0=f["y_0"], z_0=f["z_0"], fz=fz, ri=i["ri"], z_n=z_n)
t0 = time.perf_counter()
s = surftopo.run_surface_simulation(tool, d_x=Lx, m=m, d_y=Ly, n=n, z_level=f["ap"], edge_points=a.edge or grid_t, t_total=total, delta_t=a.dt)
t = time.perf_counter() - t0
z = np.asarray(s)[:, 2] * 1000.0
zc = z - z.mean()
Sq = float(np.sqrt(np.mean(zc ** 2)))
if a.dump_bin: z.astype(np.float64).tofile(a.dump_bin)
print(json.dumps(dict(method="bakhshan_pyd", Sa=float(np.mean(np.abs(zc))), Sq=Sq, Sz=float(z.max() - z.min()), n_points=int(z.size),
                      t_s=t, n_time_steps=int(total / a.dt) + 1, omp_threads=int(os.environ.get("OMP_NUM_THREADS", 0)))))
