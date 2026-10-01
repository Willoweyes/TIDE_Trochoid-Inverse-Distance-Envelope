# -*- coding: utf-8 -*-
"""
efsm_core.py  -  Pure-Python / NumPy implementation of the EFSM
(Efficient Face-milling Surface Model) Z-MAP forward-solution method.

Portable rewrite (no compiled extension, no platform binary) of the Z-MAP
surface-topography simulator. Produces the machined-surface point cloud and
ISO-25178 areal roughness (Sa, Sq, Sz, ...) for face milling.

Acceleration ideas
------------------
1. The deposited *height* of every cutting-edge point is time-independent
   (z_p = cos(gamma_p)*z_local + z_0 + delta_a*cos(gamma_p)). Only the (x,y)
   landing position changes with spindle rotation + feed -> height precomputed.
2. Each edge point sweeps a 1-D trajectory of cells while depositing a single
   CONSTANT height, so the lower-envelope reduction is a scalar np.minimum over
   the visited cells -> no sorting / no np.minimum.at needed.
3. Time processed in memory-bounded chunks so arbitrarily fine runs fit in RAM.

Reproduces the reference Z-MAP semantics: Z initialised to z_level, all edge
points deposited, grid indices clamped to the border, global minimum kept.
"""
from __future__ import annotations
import math
import time
import numpy as np

# ----------------------------------------------------------------------
# Default parameter set  (= literature Case A1)
# ----------------------------------------------------------------------
DEFAULTS = dict(
    vc=170.0, fz=0.6, ap=0.5,            # cutting conditions
    Td=10.0, ri=5.0, z_n=2,              # tool geometry
    gama_f=0.6, gama_p=0.0, phi=90.0,
    eps_r=0.011, eps_a=0.003,            # runout
    x_0=5.0, y_0=-5.0, z_0=0.0,          # workpiece-origin offset
    grid_base=200, length=5.0,           # simulation knobs
    delta_t=0.8e-6, edge_points=None,
    travel_factor_k=None,                # optional travel = k * Td
    insert="circular", nose_r=0.8,       # insert shape
)

# corner included angle (deg) for the named polygon inserts (ISO-style)
INSERT_ANGLE = dict(triangle=60.0, square=90.0, rhombic=80.0,
                    pentagon=108.0, octagon=135.0)


def compute_derived(p):
    Td = p["Td"]; vc = p["vc"]; fz = p["fz"]; z = p["z_n"]
    gb = int(p["grid_base"]); L = p["length"]
    ws = vc / (Td * 0.001 * math.pi)         # spindle speed (rpm)
    vf = fz * z * ws / 60.0                   # feed speed (mm/s)
    omega = ws * math.pi / 30.0               # angular speed (rad/s)
    grid_m = gb * 2; grid_n = gb
    Lx = L * 2; Ly = L
    grid_t = max(grid_m, grid_n) * 2
    k = p.get("travel_factor_k", None)
    travel = k * Td if k else (1.5 * Td + Ly)
    t_total = travel / vf
    edge_pts = p["edge_points"] or grid_t
    return dict(ws=ws, vf=vf, omega=omega, grid_m=grid_m, grid_n=grid_n,
                Lx=Lx, Ly=Ly, grid_t=grid_t, t_total=t_total,
                edge_points=int(edge_pts),
                n_time_steps=int(t_total / p["delta_t"]) + 1,
                dy=Ly / grid_n, dx=Lx / grid_m)


# ---- insert edge profiles -> (l, z_local) in the tool local frame ----
def edge_profile_circular(R, ap, fz, gama_f, edge_points):
    """Round insert: circular arc of radius R (matches the reference C++)."""
    Lap = math.sqrt(max(0.0, R * R - (R - ap) ** 2))
    Lfz = fz / (2.0 * math.cos(math.radians(gama_f)))
    L = max(Lap, Lfz)
    l = np.linspace(-L, L, edge_points)
    z = R - np.sqrt(np.clip(R * R - l * l, 0.0, None))
    return l, z


def edge_profile_polygon(nose_r, included_angle_deg, ap, fz, gama_f, edge_points):
    """Polygon insert: nose radius r_n tangent-blended into straight flanks.

    Corner included angle phi_c; each flank makes alpha = 90 - phi_c/2 with the
    horizontal.  |l|<=l_b: nose arc  z = r_n - sqrt(r_n^2 - l^2);
                 |l|> l_b: flank     z = z_b + (|l|-l_b)*tan(alpha)
    with tangent continuity at l_b = r_n*sin(alpha). Sampled until z reaches ap.
    """
    rn = float(nose_r)
    alpha = math.radians(90.0 - included_angle_deg / 2.0)
    l_b = rn * math.sin(alpha)
    z_b = rn * (1.0 - math.cos(alpha))
    Lfz = fz / (2.0 * math.cos(math.radians(gama_f)))
    if ap <= z_b:
        L_eng = math.sqrt(max(0.0, rn * rn - (rn - ap) ** 2))
    else:
        L_eng = l_b + (ap - z_b) / math.tan(alpha)
    L = max(L_eng, Lfz)
    l = np.linspace(-L, L, edge_points)
    al = np.abs(l)
    z = np.where(al <= l_b,
                 rn - np.sqrt(np.clip(rn * rn - al * al, 0.0, None)),
                 z_b + (al - l_b) * math.tan(alpha))
    return l, z


def build_edge_profile(p, d):
    kind = p.get("insert", "circular").lower()
    if kind in ("circular", "round", "circle"):
        return edge_profile_circular(p["ri"], p["ap"], p["fz"],
                                     p["gama_f"], d["edge_points"])
    if kind in INSERT_ANGLE:
        return edge_profile_polygon(p.get("nose_r", 0.8), INSERT_ANGLE[kind],
                                    p["ap"], p["fz"], p["gama_f"],
                                    d["edge_points"])
    if kind in ("poly", "polygon"):
        return edge_profile_polygon(p.get("nose_r", 0.8), p["included_angle"],
                                    p["ap"], p["fz"], p["gama_f"],
                                    d["edge_points"])
    raise ValueError("unknown insert geometry: %r" % kind)


# ----------------------------------------------------------------------
# Core surface simulation
# ----------------------------------------------------------------------
def run_surface(p, mem_cap=4_000_000, return_grid=False):
    """Run an EFSM surface simulation.
    Returns (surface[N,3], derived, t_sim) with columns x,y,z (mm)."""
    d = compute_derived(p)
    m, n = d["grid_m"], d["grid_n"]
    delta_x = d["Lx"] / m; delta_y = d["Ly"] / n
    inv_dx = 1.0 / delta_x; inv_dy = 1.0 / delta_y
    z_level = p["ap"]; omega = d["omega"]; vf = d["vf"]
    Nt = d["n_time_steps"]; dt = p["delta_t"]; z_n = int(p["z_n"])
    gamma_f = math.radians(p["gama_f"]); gamma_p = math.radians(p["gama_p"])
    phi = math.radians(p["phi"])
    cos_gp = math.cos(gamma_p); sin_gp = math.sin(gamma_p)

    l_loc, z_loc = build_edge_profile(p, d)
    P = l_loc.shape[0]
    w_loc = l_loc + sin_gp * z_loc
    Zflat = np.full((m + 1) * (n + 1), z_level, dtype=np.float64)

    t0 = time.perf_counter()
    chunk = max(1, int(mem_cap))
    for s in range(0, Nt, chunk):
        e = min(Nt, s + chunk)
        t = np.arange(s, e, dtype=np.float64) * dt
        for K in range(1, z_n + 1):
            dr = p["Td"] / 2.0 + (K - 1) * p["eps_r"]
            da = (K - 1) * p["eps_a"]
            var1 = phi + 2.0 * math.pi * (K - 1) / z_n - omega * t
            cos_diff = np.cos(gamma_f - var1)
            sin_diff = np.sin(gamma_f - var1)
            cx = sin_diff * (da * sin_gp) + p["x_0"]
            cy = cos_diff * (da * sin_gp) + p["y_0"] + vf * t
            zp = cos_gp * z_loc + (p["z_0"] + da * cos_gp)
            w2 = w_loc + dr
            for pi in range(P):
                xp = cos_diff * w2[pi]; xp += cx
                yp = sin_diff * w2[pi]; yp += cy
                xi = (xp * inv_dx).astype(np.int64)
                yi = (yp * inv_dy).astype(np.int64)
                np.clip(xi, 0, m, out=xi)
                np.clip(yi, 0, n, out=yi)
                idx = xi; idx *= (n + 1); idx += yi
                sub = Zflat[idx]
                np.minimum(sub, zp[pi], out=sub)
                Zflat[idx] = sub
    t_sim = time.perf_counter() - t0

    xs = np.arange(m + 1) * delta_x
    ys = np.arange(n + 1) * delta_y
    X = np.repeat(xs, n + 1); Y = np.tile(ys, m + 1)
    surface = np.column_stack([X, Y, Zflat])
    if return_grid:
        return surface, d, t_sim, Zflat.reshape(m + 1, n + 1)
    return surface, d, t_sim


# ---- ISO-25178 areal roughness (z in mm -> reported in micron) ----
def roughness(surface):
    z = surface[:, 2] * 1000.0
    z = z[~np.isnan(z)]
    if z.size == 0:
        return {k: float("nan") for k in
                ["Sa", "Sq", "Sz", "Sp", "Sv", "Ssk", "Sku"]}
    zc = z - z.mean()
    Sq = float(np.sqrt(np.mean(zc ** 2)))
    return dict(Sa=float(np.mean(np.abs(zc))), Sq=Sq,
                Sz=float(z.max() - z.min()), Sp=float(zc.max()),
                Sv=float(-zc.min()),
                Ssk=float(np.mean(zc ** 3) / Sq ** 3) if Sq > 0 else float("nan"),
                Sku=float(np.mean(zc ** 4) / Sq ** 4) if Sq > 0 else float("nan"))


def surface_grid(surface, d):
    m, n = d["grid_m"], d["grid_n"]
    X = surface[:, 0].reshape(m + 1, n + 1)
    Y = surface[:, 1].reshape(m + 1, n + 1)
    Z = surface[:, 2].reshape(m + 1, n + 1) * 1000.0
    return X, Y, Z


def nyquist_ratio(p):
    d = compute_derived(p)
    return p["fz"] / d["dy"]


if __name__ == "__main__":
    p = dict(DEFAULTS)
    surf, d, t = run_surface(p)
    r = roughness(surf)
    print("A1 self-test: points=%s  t_sim=%.3fs  fz/dy=%.1f  Nt=%d"
          % (surf.shape, t, nyquist_ratio(p), d["n_time_steps"]))
    print("  Sa=%.4f  Sz=%.4f  Sq=%.4f um" % (r["Sa"], r["Sz"], r["Sq"]))


# ----------------------------------------------------------------------
# Optional Numba parallel engine (native-code speed). Auto-validated
# against the NumPy reference on first use, so it is safe to enable.
# ----------------------------------------------------------------------
_NUMBA = None
_NUMBA_OK = False


def _get_numba_kernel():
    global _NUMBA
    if _NUMBA is not None:
        return _NUMBA
    try:
        import numba
        from numba import njit, prange

        @njit(parallel=True, fastmath=True, cache=True)
        def _kern(Nt, dt, z_n, phi, omega, gamma_f, cos_gp, sin_gp,
                  Td, eps_r, eps_a, x0, y0, z0, vf,
                  w_loc, z_loc, inv_dx, inv_dy, m, n, z_level):
            ncell = (m + 1) * (n + 1)
            nb = numba.get_num_threads()
            local = np.full((nb, ncell), z_level)
            P = w_loc.shape[0]
            per = (Nt + nb - 1) // nb
            two_pi = 2.0 * math.pi
            for b in prange(nb):
                Z = local[b]
                for t_idx in range(b * per, min(Nt, (b + 1) * per)):
                    t = t_idx * dt
                    for K in range(1, z_n + 1):
                        dr = Td / 2.0 + (K - 1) * eps_r
                        da = (K - 1) * eps_a
                        var1 = phi + two_pi * (K - 1) / z_n - omega * t
                        cd = math.cos(gamma_f - var1)
                        sd = math.sin(gamma_f - var1)
                        cx = sd * (da * sin_gp) + x0
                        cy = cd * (da * sin_gp) + y0 + vf * t
                        for pidx in range(P):
                            w2 = w_loc[pidx] + dr
                            xi = int((cd * w2 + cx) * inv_dx)
                            xi = 0 if xi < 0 else (m if xi > m else xi)
                            yi = int((sd * w2 + cy) * inv_dy)
                            yi = 0 if yi < 0 else (n if yi > n else yi)
                            idx = xi * (n + 1) + yi
                            zp = cos_gp * z_loc[pidx] + z0 + da * cos_gp
                            if zp < Z[idx]:
                                Z[idx] = zp
            out = np.full(ncell, z_level)
            for b in range(nb):
                Zb = local[b]
                for c in range(ncell):
                    if Zb[c] < out[c]:
                        out[c] = Zb[c]
            return out

        _NUMBA = _kern
    except Exception:
        _NUMBA = False
    return _NUMBA


def _run_numba_once(p, with_time=False):
    kern = _get_numba_kernel()
    d = compute_derived(p)
    m, n = d["grid_m"], d["grid_n"]
    delta_x = d["Lx"] / m; delta_y = d["Ly"] / n
    l_loc, z_loc = build_edge_profile(p, d)
    w_loc = l_loc + math.sin(math.radians(p["gama_p"])) * z_loc
    t0 = time.perf_counter()
    Zflat = kern(d["n_time_steps"], p["delta_t"], int(p["z_n"]),
                 math.radians(p["phi"]), d["omega"], math.radians(p["gama_f"]),
                 math.cos(math.radians(p["gama_p"])), math.sin(math.radians(p["gama_p"])),
                 p["Td"], p["eps_r"], p["eps_a"], p["x_0"], p["y_0"], p["z_0"], d["vf"],
                 np.ascontiguousarray(w_loc), np.ascontiguousarray(z_loc),
                 1.0 / delta_x, 1.0 / delta_y, m, n, p["ap"])
    t_sim = time.perf_counter() - t0
    xs = np.arange(m + 1) * delta_x; ys = np.arange(n + 1) * delta_y
    surface = np.column_stack([np.repeat(xs, n + 1), np.tile(ys, m + 1), Zflat])
    return (surface, d, t_sim) if with_time else surface


def run_surface_numba(p, validate=True):
    """Native-speed EFSM via Numba. Returns (surface, derived, t_sim).
    On first call, validates against the NumPy core on a tiny case so a
    broken JIT build is never trusted."""
    global _NUMBA_OK
    if not _get_numba_kernel():
        raise RuntimeError("numba not available; use run_surface (NumPy core).")
    if validate and not _NUMBA_OK:
        ps = dict(DEFAULTS); ps.update(grid_base=40, delta_t=3e-6, edge_points=120)
        a = roughness(run_surface(ps)[0]); b = roughness(_run_numba_once(ps))
        if abs(a["Sa"] - b["Sa"]) > 1e-6 or abs(a["Sz"] - b["Sz"]) > 1e-6:
            raise RuntimeError("numba kernel disagrees with NumPy reference")
        _NUMBA_OK = True
    return _run_numba_once(p, with_time=True)
