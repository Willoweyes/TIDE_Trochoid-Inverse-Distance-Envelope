# -*- coding: utf-8 -*-
"""
TIDE -- Trochoid Inverse-Distance Envelope.
===================================================================================
A closed-form analytical method for 3D face-milling surface topography and areal
roughness.  Author: Duong Duc Tri (sole author).  Personal research project.

TIDE inverts the forward FSM/EFSM simulation: rather than marching through time, it
computes -- for each workpiece grid point -- the closed-form nearest approach of the
insert-tip trochoid and reads the cut depth from the round-nose profile.  The result
is the resolution-independent geometric surface that EFSM only reaches with a fine
(slow) grid.  Detailed derivations: docs/thesis/TIDE_Thesis.docx ([REF] below).

What this module IS
-------------------
* Analytical: the cut depth at each workpiece grid point is obtained from the
  *closed-form* nearest approach of the insert-tip trochoid ([REF] Prop. 1),
  optionally polished by 1-2 Newton steps on the exact stationarity equation
  ([REF] Eq. 23, sign-corrected -- see ERRATUM below).  There is NO sampling of
  the trajectory in time and NO discretization of the cutting edge.
* Complexity O(N_grid * N_pass * z_n) with a small constant ([REF] Eq. 24):
  per grid point, insert and pass, a fixed number (<= 6) of candidate loops is
  examined, independent of spindle speed, swept length and grid resolution.
* Pure kinematic geometry.  NO empirical calibration factor is applied
  (an earlier development build used a least-squares eta fitted on the same cases
  it was scored on; that circular factor was removed deliberately).

What this module is NOT
-----------------------
* Not a vibration / tool-wear / elastic-recovery model: it predicts the ideal
  swept-envelope surface, like the FSM/EFSM it replaces ([REF] Sec. 11).
* Not a time-sampled approximation: methods that sample each trochoid at many
  time points and query a k-d tree silently re-introduce the N_t factor the
  theory removes.  TIDE is fully closed-form.

ERRATUM with respect to [REF]
-----------------------------
[REF] Eq. (23) prints the stationarity condition of ||P - tip(t)||^2 with a
wrong sign on its first term.  The correct equation (verified symbolically,
see src/tide/symbolic_checks.py) is

    g(t) = (X - x_tip) rho_K w sin(Theta_K) + (Y - y_tip)(v_f + rho_K w cos(Theta_K)) = 0.

The Newton step below uses this corrected g and its exact derivative g'.

Conventions (identical to the legacy code and the EFSM benchmark)
-----------------------------------------------------------------
* Insert index K is 0-based in code; [REF] is 1-based, so code `K` == [REF] K-1.
* Units: mm and seconds internally; the returned surface is in micrometres.
* Grid: (m+1) x (n+1) nodes, m = 2*grid_base (X, cross-feed), n = grid_base
  (Y, feed); spacing Lx/m, Ly/n; tool centre line X = x0 (+ p*ae per pass).
* Time window: the tool centre travels y0 -> y0 + v_f*t_total with
  t_total = (1.5*D + Ly)/v_f; trochoid loops outside [0, t_total] do not exist
  and are excluded (exactly like the FSM, which only simulates that window).

Requirements: numpy only.  (scipy is used by validate_tide.py, not here.)

Run on one benchmark case:
    python src/tide/tide.py benchmark/test_cases/paper/paper_A1.json
"""
from dataclasses import dataclass
from pathlib import Path
import sys
import os
import time

import numpy as np

METHOD_NAME = "tide_closed_form"

# roughness.py lives next to this module (src/tide/).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ============================================================================
# 1.  Parameters  ([REF] Sec. 2, Eq. 1)
# ============================================================================
@dataclass(frozen=True)
class Kinematics:
    """All quantities derived from a benchmark config (mm, rad, s)."""
    D: float          # cutting diameter
    R: float          # insert corner radius
    z_n: int          # number of inserts
    f_z: float        # feed per tooth                       [REF] Eq. 1
    a_p: float        # axial depth of cut
    eps_r: float      # radial run-out step between inserts
    eps_a: float      # axial  run-out step between inserts
    gamma_f: float    # radial rake angle [rad]
    gamma_p: float    # axial  rake angle [rad]
    phi: float        # initial phase of insert K=0 [rad]
    x0: float         # tool-centre X (pass p adds p*a_e)
    y0: float         # tool-centre Y at t=0
    z0: float         # workpiece datum
    a_e: float        # step-over between passes (multi-pass only)
    n_pass: int       # number of passes
    omega: float      # spindle angular velocity [rad/s]
    v_f: float        # feed velocity [mm/s]
    t_total: float    # machining-time window [s]
    grid_m: int       # X-nodes - 1
    grid_n: int       # Y-nodes - 1
    Lx: float         # grid extent in X
    Ly: float         # grid extent in Y

    @property
    def T_rev(self):
        """Spindle revolution period; centre advances f_z*z_n per revolution."""
        return 2.0 * np.pi / self.omega


def kinematics_from_config(cfg) -> Kinematics:
    """Parse a benchmark JSON config ([REF] Eq. 1 and the FSM grid conventions).

    n_s = v_c / (pi D)   [rev/min],   omega = 2 pi n_s / 60,
    v_f = f_z z_n n_s / 60           => v_f * T_rev = f_z z_n  (loop spacing).
    """
    inp, fix, grid = cfg["input"], cfg["fixed"], cfg["grid"]
    D, R, z_n, f_z = fix["Td"], inp["ri"], int(inp["z"]), inp["fz"]
    assert D > 0 and R > 0 and z_n >= 1 and f_z > 0, "bad geometry in config"
    assert R >= fix["ap"], "corner radius R must be >= depth of cut ap"

    n_s = inp["vc"] / (D * 1e-3 * np.pi)        # spindle speed [rpm]
    omega = n_s * np.pi / 30.0                  # [rad/s]
    v_f = f_z * z_n * n_s / 60.0                # [mm/s]

    gb, length = grid["grid_base"], grid["length"]
    Lx, Ly = 2.0 * length, length
    passes = cfg.get("passes", {})              # optional multi-pass block
    return Kinematics(
        D=D, R=R, z_n=z_n, f_z=f_z, a_p=fix["ap"],
        eps_r=inp["eps_r"], eps_a=inp["eps_a"],
        gamma_f=np.radians(fix["gama_f"]), gamma_p=np.radians(fix["gama_p"]),
        phi=np.radians(fix["phi"]),
        x0=fix["x_0"], y0=fix["y_0"], z0=fix["z_0"],
        a_e=float(passes.get("ae", 0.0)), n_pass=int(passes.get("n_pass", 1)),
        omega=omega, v_f=v_f, t_total=(1.5 * D + Ly) / v_f,
        grid_m=2 * gb, grid_n=gb, Lx=Lx, Ly=Ly)


def workpiece_grid(k: Kinematics):
    """Workpiece nodes (X: cross-feed, Y: feed), identical to the legacy grid."""
    dx, dy = k.Lx / k.grid_m, k.Ly / k.grid_n
    return np.meshgrid(np.arange(k.grid_m + 1) * dx,
                       np.arange(k.grid_n + 1) * dy, indexing="ij")


# ============================================================================
# 2.  The depth profile zeta and the engaged half-width  ([REF] Eqs. 2-3, 11-12)
# ============================================================================
def zeta(l, R, cos_gp):
    """Height of an edge point at lateral offset |l| above the tip ([REF] Eq. 11):

        zeta(l) = (R - sqrt(R^2 - l^2)) * cos(gamma_p),  0 <= l <= R.

    Strictly increasing on [0, R) ([REF] Eq. 12) -- this monotonicity is what
    makes "deepest cut" equivalent to "nearest approach" (Backward Reduction).
    """
    return (R - np.sqrt(np.clip(R * R - l * l, 0.0, None))) * cos_gp


def engaged_half_width(k: Kinematics):
    """Engaged half-width DeltaL of the edge chord ([REF] Eq. 3):

        DeltaL = max( sqrt(R^2-(R-a_p)^2),  f_z/(2 cos gamma_f) ).

    Term 1 solves zeta(DeltaL)=a_p (cut branch joins the uncut plane z0+a_p
    continuously); term 2 guarantees adjacent feed marks, spaced f_z, connect.
    """
    t1 = np.sqrt(max(k.R**2 - (k.R - k.a_p)**2, 0.0))
    t2 = k.f_z / (2.0 * np.cos(k.gamma_f))
    return max(t1, t2)


# ============================================================================
# 3.  Closed-form nearest approach to one insert's trochoid
#     ([REF] Def. 1, Thm. 2, Eqs. 16-19 (Prop. 1), Eq. 23 corrected)
# ============================================================================
def _tip(t, k: Kinematics, rho, phi0, x_c):
    """Insert-tip trochoid ([REF] Eq. 10):
        Theta(t) = phi0 - omega t,
        x(t) = x_c + rho cos Theta,   y(t) = y0 + v_f t - rho sin Theta.
    """
    Th = phi0 - k.omega * t
    return x_c + rho * np.cos(Th), k.y0 + k.v_f * t - rho * np.sin(Th), Th


def _candidate_times(Px, Py, k: Kinematics, rho, phi0, x_c):
    """O(1) candidate contact times per grid point ([REF] Prop. 1, spelt out).

    Geometry: when the tip passes exactly through P, the centre C = (x_c, Y_c)
    satisfies |P - C| = rho, i.e. (with h = X - x_c the signed distance of P
    from the feed line):
        Y_c = Y + s*sqrt(rho^2 - h^2),   s in {+1,-1}.            [REF Eq. 18]
    At such a crossing the tip must POINT AT P:
        (cos Theta, -sin Theta) = (P - C)/rho  =>  Theta* = atan2(s*disc, h).
    The spindle phase fixes the discrete times realising Theta*:
        t = (phi0 - Theta*)/omega + j*T_rev,  j integer,          [REF Eq. 16]
    and the j actually near the crossing is j ~ (t_cross - t_base)/T_rev with
    t_cross = (Y_c - y0)/v_f.  We return j-1, j, j+1 for both signs s
    (6 candidates; the convexity argument of [REF] Prop. 1 guarantees the
    minimiser is among the lattice neighbours of the continuous solutions).
    For |h| >= rho, disc = 0 and both signs coincide with the perpendicular
    foot, which is the convex minimiser in that regime ([REF] Prop. 1, case 2).
    Times are clipped to the machining window [0, t_total] so that loops which
    are never machined are not consulted (the FSM behaves identically).
    """
    h = Px - x_c
    disc = np.sqrt(np.clip(rho * rho - h * h, 0.0, None))
    cands = []
    for s in (+1.0, -1.0):
        Yc = Py + s * disc
        Theta_star = np.arctan2(s * disc, h)
        t_base = (phi0 - Theta_star) / k.omega
        j = np.round(((Yc - k.y0) / k.v_f - t_base) / k.T_rev)
        for dj in (-1.0, 0.0, 1.0):
            cands.append(t_base + (j + dj) * k.T_rev)
    t = np.stack(cands)                          # shape (6, N)
    return np.clip(t, 0.0, k.t_total)


def _realign(t, Px, Py, k: Kinematics, rho, phi0, x_c, passes=1, gf=0.0):
    """Re-solve the pointing condition ON the snapped loop (fixed point).

    The candidate times of `_candidate_times` use the contact angle of the
    *continuous* crossing Y_c; after snapping to a discrete loop the centre
    has moved by up to f_z z_n/2, so the tip-points-at-P angle changes by
    O(f_z z_n / D).  One fixed-point pass

        Theta* <- atan2(y_c(t) - Y, X - x_c),  t <- (phi0 - Theta*)/omega + j T_rev

    (j chosen to stay on the same loop) restores the alignment to second
    order; it is the discrete analogue of choosing the reference centre at
    the contact phase in [REF] Sec. 8.2.
    """
    for _ in range(passes):
        yc = k.y0 + k.v_f * t
        Theta_star = np.arctan2(yc - Py, Px - x_c)
        t_base = (phi0 - gf - Theta_star) / k.omega
        t = t_base + np.round((t - t_base) / k.T_rev) * k.T_rev
    return np.clip(t, 0.0, k.t_total)


def _newton_polish(t, Px, Py, k: Kinematics, rho, phi0, x_c, steps):
    """Newton iteration on the exact stationarity of ||P - tip(t)||^2.

    CORRECTED [REF] Eq. (23)  (first-term sign fixed, verified symbolically):
        g(t)  = (X-x) rho w sinTh + (Y-y)(v_f + rho w cosTh),
        g'(t) = -rho^2 w^2 sin^2 Th - (X-x) rho w^2 cosTh
                - (v_f + rho w cosTh)^2 + (Y-y) rho w^2 sinTh.
    The realigned candidates already sit in the quadratic basin of the
    nearest root ([REF] Sec. 8.3), so 1-2 steps reach round-off in practice.
    """
    w = k.omega
    for _ in range(steps):
        x, y, Th = _tip(t, k, rho, phi0, x_c)
        sTh, cTh = np.sin(Th), np.cos(Th)
        A, B = Px - x, Py - y
        g = A * rho * w * sTh + B * (k.v_f + rho * w * cTh)
        gp = (-(rho * w * sTh) ** 2 - A * rho * w * w * cTh
              - (k.v_f + rho * w * cTh) ** 2 + B * rho * w * w * sTh)
        t = t - g / np.where(np.abs(gp) > 1e-30, gp, np.inf)
        t = np.clip(t, 0.0, k.t_total)           # window is part of the problem
    return t


def dstar_insert(Px, Py, k: Kinematics, K, p=0, newton_steps=2, mode="newton"):
    """Nearest distance d*_K(P) from each grid point to insert K's trochoid.

    mode="circle":  pure [REF] Prop. 1 -- distance to the circular loop,
                    | sqrt(h^2 + (Y - Y_c^(m))^2) - rho_K |   [REF Eq. 17/19];
                    carries the O(f_z z_n / D) circular-arc error [REF Eq. 22].
    mode="aligned": evaluate the EXACT trochoid at the realigned candidate
                    times (no Newton); error is one order smaller.
    mode="newton":  aligned + `newton_steps` Newton polishes => exact to
                    round-off ([REF] Sec. 8.3).  DEFAULT.
    """
    rho = k.D / 2.0 + K * k.eps_r                # [REF] rho_K, code K 0-based
    phi0 = k.phi + 2.0 * np.pi * K / k.z_n
    x_c = k.x0 + p * k.a_e
    t = _candidate_times(Px, Py, k, rho, phi0, x_c)
    if mode == "circle":
        Yc = k.y0 + k.v_f * t                    # loop centres at contact phase
        h = Px - x_c
        d = np.abs(np.hypot(h, Py - Yc) - rho)
    else:
        t = _realign(t, Px, Py, k, rho, phi0, x_c)
        if mode == "newton":
            t = _newton_polish(t, Px, Py, k, rho, phi0, x_c, newton_steps)
        x, y, _ = _tip(t, k, rho, phi0, x_c)
        d = np.hypot(Px - x, Py - y)
    return d.min(axis=0)                         # best of the <= 6 candidates


def contact_offset(Px, Py, k: Kinematics, K, p=0, newton_steps=0, cross_steps=4):
    """Exact signed edge offset l at which insert K cuts P (one value per candidate).

    In the reference geometry every edge point of insert K lies on one ray from
    the tool centre c(t), direction u(t) = (cos a, -sin a), a = Theta(t) - gamma_f,
    at radius rho_K + l.  P is cut when the ray passes through P:

        f(t) = u(t) x (P - c(t)) = 0,   u(t).(P - c(t)) > 0,

    and then l = u.(P - c) - rho_K exactly (gamma_p = 0; for gamma_p != 0 the
    radius is rho_K + l + sin(gp) z0(l), inverted by Newton).  The tip's closest
    approach t* (candidates + realignment + Newton on g) seeds a Newton
    iteration on f; no radial-direction assumption is made.
    Returns l with shape (6, N); non-crossings are +inf.
    """
    rho = k.D / 2.0 + K * k.eps_r
    phi0 = k.phi + 2.0 * np.pi * K / k.z_n
    x_c = k.x0 + p * k.a_e
    t = _candidate_times(Px, Py, k, rho, phi0, x_c)
    t = _realign(t, Px, Py, k, rho, phi0, x_c, gf=k.gamma_f)   # aim the edge ray at P
    if newton_steps:
        t = _newton_polish(t, Px, Py, k, rho, phi0, x_c, newton_steps)
    w, gf = k.omega, k.gamma_f
    for _ in range(cross_steps):
        a = phi0 - w * t - gf
        ca, sa = np.cos(a), np.sin(a)
        dx, dy = Px - x_c, Py - (k.y0 + k.v_f * t)
        f = ca * dy + sa * dx
        fp = w * sa * dy - ca * k.v_f - w * ca * dx
        t = t - f / np.where(np.abs(fp) > 1e-30, fp, np.inf)
    a = phi0 - w * t - gf
    dx, dy = Px - x_c, Py - (k.y0 + k.v_f * t)
    r = np.cos(a) * dx - np.sin(a) * dy
    ok = (t >= 0.0) & (t <= k.t_total) & (r > 0)
    l = r - rho
    sp_ = np.sin(k.gamma_p)
    if abs(sp_) > 1e-14:                        # r = rho + l + sin(gp) z0(|l|)
        for _ in range(4):
            al = np.minimum(np.abs(l), k.R)
            rt = np.sqrt(np.maximum(k.R**2 - al**2, 0.0))
            F = rho + l + sp_ * (k.R - rt) - r
            dF = 1.0 + sp_ * np.sign(l) * al / np.where(rt > 1e-300, rt, 1.0)
            l = l - F / dF
    return np.where(ok, l, np.inf)


def lateral_offset(d, rho, R, gamma_f, gamma_p):
    """Exact tip-distance -> arc-offset conversion (vectorised in d).

    zeta() takes the ARC OFFSET l along the cutting edge, while d* is an
    IN-PLANE distance from the tip trajectory.  Under the transform chain of
    the reference formulation the edge point at arc offset l sits at in-plane
    radius

        r(l) = | ( rho + l cos gf + sin gf sin gp z0(l),
                  -l sin gf      + cos gf sin gp z0(l) ) |,
        z0(l) = R - sqrt(R^2 - l^2),

    and the contacting point satisfies r(l) = rho + d*.  The common relation
    l = d*/cos(gamma_f) is the first-order truncation of that inverse and omits
    a term of size l^2 sin^2(gf)/(2 rho).  For gamma_p = 0 the inverse is a
    quadratic in l and closes exactly; otherwise three Newton steps from that
    root close it.  r(l) is monotone increasing on [0, R].
    """
    cf, sf = np.cos(gamma_f), np.sin(gamma_f)
    cp, sp_ = np.cos(gamma_p), np.sin(gamma_p)
    target = rho + np.asarray(d, dtype=float)
    l = -rho * cf + np.sqrt(np.maximum(rho * rho * cf * cf
                                       - rho * rho + target * target, 0.0))
    if abs(sp_) < 1e-14:
        return np.minimum(l, R)
    for _ in range(3):
        rt = np.sqrt(np.maximum(R * R - l * l, 0.0))
        z0 = R - rt
        dz0 = np.where(rt > 1e-300, l / np.where(rt > 1e-300, rt, 1.0), 0.0)
        ax, dax = rho + cf * l + sf * sp_ * z0, cf + sf * sp_ * dz0
        ay, day = -sf * l + cf * sp_ * z0, -sf + cf * sp_ * dz0
        r = np.hypot(ax, ay)
        dr = (ax * dax + ay * day) / np.where(r > 1e-300, r, 1.0)
        step = np.where(np.abs(dr) > 1e-300, (r - target) / dr, 0.0)
        l = np.clip(l - step, 0.0, R)
    return np.minimum(l, R)


# ============================================================================
# 4.  Surface = lower envelope over inserts and passes  ([REF] Thm. 2, Cor. 1)
# ============================================================================
def surface_zmap(cfg, mode="newton", newton_steps=2):
    """Machined surface in micrometres ([REF] Eqs. 14-15):

        z_surf(P) = min over (p, K) of [ z0 + K eps_a + zeta(l*_{p,K}) ],
        l*_{p,K}  = lateral_offset(d*_{p,K}(P), rho_K, R, gamma_f, gamma_p)

    deposited only where l* <= DeltaL; otherwise the point keeps the stock
    height z0 + a_p.  The min is taken AFTER adding K*eps_a ([REF] Rem. 4:
    with axial run-out the deepest cut need not come from the nearest insert).

    Returns (z_um 2D array, Kinematics, info dict).
    """
    k = kinematics_from_config(cfg)
    X, Y = workpiece_grid(k)
    Px, Py = X.ravel(), Y.ravel()
    dL = engaged_half_width(k)
    cos_gf, cos_gp = np.cos(k.gamma_f), np.cos(k.gamma_p)

    z = np.full(Px.shape, np.inf)
    for p in range(k.n_pass):
        for K in range(k.z_n):
            l = contact_offset(Px, Py, k, K, p, newton_steps)
            al = np.abs(l)
            depth = np.where(al <= dL,
                             zeta(np.minimum(al, k.R), k.R, cos_gp) + K * k.eps_a,
                             np.inf).min(axis=0)
            np.minimum(z, depth, out=z)

    uncut = np.isinf(z)
    z[uncut] = k.a_p                             # never reached -> stock plane
    info = {"uncut_frac": float(uncut.mean()), "DeltaL": dL, "mode": mode}
    if info["uncut_frac"] > 0.02:
        print(f"[bsm warning] {info['uncut_frac']*100:.1f}% of grid points uncut "
              f"(stock height). Check the measurement window.", file=sys.stderr)
    return ((k.z0 + z) * 1000.0).reshape(X.shape), k, info


# ============================================================================
# 5.  Benchmark-harness entry point (same API as the other methods)
# ============================================================================
def simulate(cfg, mode="newton", newton_steps=2):
    from roughness import compute_roughness     # shared ISO 25178-2 metrics
    t0 = time.time()
    z_um, _, info = surface_zmap(cfg, mode, newton_steps)
    res = compute_roughness(z_um.ravel())
    res["t_sim_s"] = round(time.time() - t0, 4)
    res["method"] = METHOD_NAME if mode == "newton" else f"{METHOD_NAME}_{mode}"
    res["uncut_frac"] = round(info["uncut_frac"], 5)
    return res


if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser(description="Closed-form TIDE (no time sampling, no calibration).")
    ap.add_argument("config", help="path to a benchmark test-case JSON")
    ap.add_argument("--mode", choices=["newton", "aligned", "circle"], default="newton",
                    help="newton = exact (default); aligned/circle expose the error ladder")
    ap.add_argument("--newton-steps", type=int, default=2)
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    o = simulate(cfg, mode=args.mode, newton_steps=args.newton_steps)
    print(f"[{cfg['case_name']}]  Sa={o['Sa']:.4f}  Sq={o['Sq']:.4f}  Sz={o['Sz']:.4f} um  "
          f"t={o['t_sim_s']}s  uncut={o['uncut_frac']*100:.2f}%  ({o['n_points']} pts)  "
          f"[{o['method']}]")
