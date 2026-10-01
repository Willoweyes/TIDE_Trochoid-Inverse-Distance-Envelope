"""Every experiment of the paper.  Usage:  python scripts/run_experiments.py [step ...]
Steps: verify sampling ladder cost resolution surfaces scaling   (default: all)
All timings are single-thread wall-clock (std::chrono inside tide_tools)."""
import sys, csv, json, subprocess
import numpy as np
from common import ROOT, TMP, ENV, exe, case_cfg, edge_speed, run, points
import tide as T
RAW, PROC = ROOT/"data/raw", ROOT/"data/processed"
RAW.mkdir(parents=True, exist_ok=True); PROC.mkdir(parents=True, exist_ok=True); (RAW/"surfaces").mkdir(exist_ok=True)
ALL = ["verify", "sampling", "ladder", "cost", "resolution", "surfaces", "scaling"]
steps = sys.argv[1:] or ALL
def save_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def matched_dt(c, gb): return (c["grid"]["length"]/gb)/(2*edge_speed(c))     # r = v_c dt / dy = 0.5

# ---- Table 2: TIDE (Python + C++) vs brute-force crossing search, random + edge bands
if "verify" in steps:
    rng = np.random.default_rng(7); lines = []
    for case in ["A1", "A2", "A3"]:
        cfgf, cfg = case_cfg(case); k = T.kinematics_from_config(cfg); dL = T.engaged_half_width(k)
        gf, gp = k.gamma_f, k.gamma_p
        zeta = lambda l: (k.R-np.sqrt(k.R**2-np.minimum(np.abs(l), k.R)**2))*np.cos(gp)
        def f_fp(t, K, px, py):
            a = k.phi+2*np.pi*K/k.z_n-k.omega*t-gf; ux, uy = np.cos(a), -np.sin(a)
            dx, dy = px-k.x0, py-(k.y0+k.v_f*t); return ux*dy-uy*dx, ux*dx+uy*dy
        def oracle(px, py, per_rev=4000):          # dense sign-change scan + 80 bisections
            out = np.empty(px.size)
            for i in range(px.size):
                best = np.inf
                for K in range(k.z_n):
                    rho = k.D/2+K*k.eps_r; tc = (py[i]-k.y0)/k.v_f; span = (rho+dL+1)/k.v_f
                    t0, t1 = max(0, tc-span), min(k.t_total, tc+span)
                    ts = np.linspace(t0, t1, int((t1-t0)/k.T_rev*per_rev)+2); f, dot = f_fp(ts, K, px[i], py[i])
                    for j in np.where((np.sign(f[:-1]) != np.sign(f[1:])) & (dot[:-1] > 0))[0]:
                        a, b, fa = ts[j], ts[j+1], f[j]
                        for _ in range(80):
                            m = .5*(a+b); fm = f_fp(m, K, px[i], py[i])[0]
                            if np.sign(fm) == np.sign(fa): a, fa = m, fm
                            else: b = m
                        l = f_fp(.5*(a+b), K, px[i], py[i])[1]-rho
                        if abs(l) <= dL: best = min(best, zeta(l)+K*k.eps_a)
                out[i] = best if np.isfinite(best) else k.a_p
            return out
        def tide_py(px, py):
            z = np.full(px.shape, np.inf)
            for K in range(k.z_n):
                l = T.contact_offset(px, py, k, K, cross_steps=4)
                z = np.minimum(z, np.where(np.abs(l) <= dL, zeta(l)+K*k.eps_a, np.inf).min(axis=0))
            z[np.isinf(z)] = k.a_p; return z
        N = 250; rho0 = k.D/2
        sets = {"random": (rng.uniform(0, k.Lx, N), rng.uniform(0, k.Ly, N))}
        for nm, (lo, hi) in (("tip-sweep edge", (rho0-.01, rho0+.01)), ("outer band", (rho0-.3, rho0+dL+.05)),
                             ("cut edge", (rho0+dL-.01, rho0+dL+.01))):
            u = rng.uniform(lo, hi, N); side = rng.choice([-1, 1], N)
            sets[nm] = (k.x0+side*u, rng.uniform(.5, k.Ly-.5, N))
        for nm, (px, py) in sets.items():
            zo = oracle(px, py)*1e3; ncut = int((zo < k.a_p*1e3-1e-9).sum())
            ep = np.max(np.abs(tide_py(px, py)*1e3-zo))
            ec = np.max(np.abs(points(cfgf, px, py)-(k.z0*1e3+zo)))
            lines.append(f"{case} {nm} n={N} cut={ncut} max|python-oracle|={ep:.1e} um max|C++-oracle|={ec:.1e} um")
            print(lines[-1], flush=True)
    (PROC/"edge_verification.txt").write_text("\n".join(lines)+"\n")

# ---- TIDE sampling study (doubling grids) -> Sect. 5.1
if "sampling" in steps:
    rows = []
    for case in ["A1", "A2", "A3"]:
        cfgf, _ = case_cfg(case); prev = None
        for gb in [25, 50, 100, 200, 400, 800, 1600]:
            f = TMP/"s.bin"; r = run("tide", cfgf, gb, dump=f); z = np.fromfile(f).reshape(2*gb+1, gb+1)
            d = float(np.max(np.abs(z[::2, ::2]-prev))) if prev is not None else float("nan"); prev = z
            rows.append(dict(case=case, grid_base=gb, n_points=z.size, Sa=r["Sa"], Sq=r["Sq"], Sz=r["Sz"],
                             max_diff_shared_nodes_um=d, t_s=r["t_s"])); print("sampling", case, gb, r["Sa"], r["t_s"], flush=True)
    save_csv(RAW/"tide_sampling.csv", rows)

# ---- Fig. 2: FSM refinement ladder, released dt=8e-7 vs dt matched to the grid
if "ladder" in steps:
    out = []
    for case in ["A1", "A2", "A3"]:
        cfgf, c = case_cfg(case); ve = edge_speed(c)
        for gb in [25, 50, 100, 200, 400, 800, 1600]:
            dy = c["grid"]["length"]/gb; tb = TMP/"t.bin"; t = run("tide", cfgf, gb, dump=tb); T_ = np.fromfile(tb)
            for tag, dt in (("released", 8e-7), ("matched", matched_dt(c, gb))):
                fb = TMP/"f.bin"; r = run("fsm", cfgf, gb, dt=dt, dump=fb); d = np.fromfile(fb)-T_
                out.append(dict(case=case, grid_base=gb, dt=dt, regime=tag, ratio=ve*dt/dy, fsm_Sa=r["Sa"], fsm_Sz=r["Sz"],
                                fsm_t=r["t_s"], tide_Sa=t["Sa"], tide_Sz=t["Sz"], tide_t=t["t_s"],
                                frac_above_2um=float((d > 2).mean()), n_above_2um=int((d > 2).sum()),
                                max_above=float(d.max()), n_points=int(d.size)))
                print("ladder", case, gb, tag, "%.4f/%.4f %.2fs" % (r["Sa"], t["Sa"], r["t_s"]), flush=True)
    save_csv(RAW/"ladder_dt_regimes.csv", out)

# ---- Fig. 4b: runtime vs error (matched dt, doubling grids)
if "cost" in steps:
    rows = []
    for case in ["A1", "A2", "A3"]:
        cfgf, c = case_cfg(case)
        for gb in [25, 50, 100, 200, 400, 800, 1600]:
            dt = matched_dt(c, gb); f = run("fsm", cfgf, gb, dt=dt); t = run("tide", cfgf, gb)
            rows.append(dict(case=case, grid_base=gb, n_points=t["n_points"], dt=dt, tide_Sa=t["Sa"], tide_Sz=t["Sz"],
                             tide_t=t["t_s"], fsm_Sa=f["Sa"], fsm_Sz=f["Sz"], fsm_t=f["t_s"]))
            print("cost", case, gb, "%.3fs %.3fs" % (t["t_s"], f["t_s"]), flush=True)
    save_csv(RAW/"cost_accuracy_dense.csv", rows)

# ---- Sect. 5.3 (Time step): sampling condition r < r_max (A3 at 1600 cells, A1 at 3200 cells)
if "resolution" in steps:
    rows = []
    for case, gb, dts in (("A3", 400, [8e-7]), ("A3", 800, [8e-7]), ("A3", 1600, [8e-7, 2e-7]), ("A1", 3200, [8e-7, 2e-7])):
        cfgf, c = case_cfg(case); ve = edge_speed(c); dy = c["grid"]["length"]/gb
        tb = TMP/"t.bin"; t = run("tide", cfgf, gb, dump=tb); T_ = np.fromfile(tb)
        rows.append(dict(case=case, grid_base=gb, delta_t="", ratio=f"{ve*8e-7/dy:.2f}", method="tide",
                         Sa_um=t["Sa"], Sz_um=t["Sz"], t_s=t["t_s"], cells_gt2um_above_tide=""))
        for dt in dts:
            fb = TMP/"f.bin"; r = run("fsm", cfgf, gb, dt=dt, dump=fb); d = np.fromfile(fb)-T_
            rows.append(dict(case=case, grid_base=gb, delta_t=dt, ratio=f"{ve*dt/dy:.2f}", method="fsm",
                             Sa_um=r["Sa"], Sz_um=r["Sz"], t_s=r["t_s"], cells_gt2um_above_tide=int((d > 2).sum())))
            print("resolution", rows[-1], flush=True)
    save_csv(RAW/"resolution_condition.csv", rows)

# ---- Fig. 3: height maps and profile of A1 (FSM with matched dt)
if "surfaces" in steps:
    cfgf, c = case_cfg("A1"); d = {}
    for meth in ("fsm", "tide"):
        for gb in (25, 100, 400):
            f = TMP/"z.bin"; run(meth, cfgf, gb, dt=matched_dt(c, gb) if meth == "fsm" else None, dump=f)
            m, n = 2*gb, gb; Z = np.fromfile(f).reshape(m+1, n+1)
            d[f"{meth}{gb}_X"] = np.repeat((np.arange(m+1)*(2*c["grid"]["length"]/m))[:, None], n+1, 1)
            d[f"{meth}{gb}_Y"] = np.repeat((np.arange(n+1)*(c["grid"]["length"]/n))[None, :], m+1, 0)
            d[f"{meth}{gb}_Z"] = Z
    d["tide_X"], d["tide_Y"], d["tide_Z"] = d["tide400_X"], d["tide400_Y"], d["tide400_Z"]
    np.savez_compressed(RAW/"surfaces/A1_fsm_vs_tide.npz", **d); print("surfaces ok")

# ---- Sect. 5.2: TIDE runtime vs grid, v_c (50-2400 m/min) and f_z (0.05-1.2 mm)
if "scaling" in steps:
    out = subprocess.run([exe("bench_scaling"), "--reps", "3"], capture_output=True, text=True, env=ENV, check=True).stdout
    (RAW/"scaling.csv").write_text(out); print(out)
