"""Fig. 1 -- (a) cutting part of a round insert, (b) the six candidate
revolutions at a query point P, (c) Newton iteration on the crossing
condition f(t) = 0 with the full (moving-centre) motion.  Drawn to the
Springer LNCS text width (4.8 in)."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import style
plt = style.apply()
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
ROOT = Path(__file__).resolve().parents[2]
CT, CF, CG = style.C_EXACT, style.C_FORWARD, style.C_GREY
FS, FSS = 7.2, 6.4                                   # label / small label size
halo = dict(fc="white", ec="none", pad=0.4, alpha=0.9)

fig, ax = plt.subplots(1, 3, figsize=(style.W_FULL, 1.8),
                       gridspec_kw=dict(width_ratios=[1.05, 0.9, 1.05]))
for a in ax:
    a.set_aspect("equal"); a.axis("off")

# ---------------- (a) cutting part of the round insert -------------------
a = ax[0]; R = 1.0; ap = 0.36
a.add_patch(Rectangle((-0.95, 0), 1.9, ap, fc="#e8edf3", ec="none"))
a.plot([-0.95, 0.95], [ap, ap], color=CG, lw=0.6, ls="--")
a.text(0.96, ap, "stock", fontsize=FSS, color=CG, va="center")
a.plot([-0.95, 0.95], [0, 0], color=CG, lw=0.6)
th = np.linspace(np.deg2rad(-135), np.deg2rad(-45), 200)
a.plot(R*np.cos(th), R + R*np.sin(th), color="k", lw=1.3)
a.plot([0, R*np.cos(th[-1])], [R, R + R*np.sin(th[-1])], color=CG, lw=0.5)
a.text(0.44, 0.74, "$R$", fontsize=FS)
a.plot(0, R, "k+", ms=4)
a.text(-0.07, R + 0.03, "insert centre", fontsize=FSS, ha="right", va="bottom", color=CG)
a.plot(0, 0, "o", color=CT, ms=3.8, zorder=5)
a.text(0, -0.15, "tip", fontsize=FSS, color=CT, ha="center", va="top")
l = 0.70; zl = R - np.sqrt(R**2 - l**2)
a.plot(l, zl, "o", color=CF, ms=3.8, zorder=5)
yb = -0.36
a.plot([0, l], [yb, yb], color=CF, lw=0.9)
for xx in (0, l): a.plot([xx, xx], [yb + 0.04, yb - 0.04], color=CF, lw=0.9)
a.text(l/2, yb - 0.07, r"$\ell$", color=CF, ha="center", va="top", fontsize=FS + 0.5)
a.plot([l, l + 0.14], [zl, zl], color=CF, lw=0.5, ls=":")
a.annotate("", xy=(l + 0.1, 0.0), xytext=(l + 0.1, zl),
           arrowprops=dict(arrowstyle="<->", lw=0.6, color=CF, mutation_scale=5, shrinkA=0, shrinkB=0))
a.text(l + 0.15, zl/2, r"$\zeta(\ell)$", color=CF, fontsize=FS, va="center")
a.annotate("", xy=(-0.86, 0), xytext=(-0.86, ap), arrowprops=dict(arrowstyle="<->", lw=0.6))
a.text(-0.81, ap/2, r"$a_p$", fontsize=FS, va="center")
xL = np.sqrt(R**2 - (R - ap)**2)
a.annotate("", xy=(-xL, ap + 0.07), xytext=(xL, ap + 0.07), arrowprops=dict(arrowstyle="<->", lw=0.6, color=CG))
a.text(-0.33, ap + 0.12, r"engaged edge $2\Delta L$", fontsize=FSS, ha="center", color=CG)
a.set_xlim(-1.05, 1.3); a.set_ylim(-0.62, 1.22)

# ---------------- (b) the six candidate revolutions ----------------------
a = ax[1]; rho = 1.0; step = 0.3; P = np.array([0.72, 0.0])
disc = np.sqrt(rho**2 - P[0]**2)
thh = np.linspace(-np.pi/2, np.pi/2, 200)
a.plot([0, 0], [-1.78, 1.78], color=CG, lw=0.6, ls="-.")
for s_ in (+1, -1):
    yc = P[1] + s_*disc
    for k in (-2, 2):
        a.plot(rho*np.cos(thh), yc + k*step + rho*np.sin(thh), color="0.82", lw=0.45)
        a.plot(0, yc + k*step, "o", color="0.75", ms=1.8)
    for k in (-1, 0, 1):
        a.plot(rho*np.cos(thh), yc + k*step + rho*np.sin(thh), color="k", lw=0.6)
        a.plot(0, yc + k*step, "s", color="k", ms=2.4, zorder=5)
    a.text(-0.1, yc, f"$s={'+' if s_ > 0 else '-'}1$", fontsize=FSS, ha="right", va="center")
a.plot(*P, "*", color=CF, ms=8, mec="k", mew=0.4, zorder=8)
a.text(P[0] + 0.36, P[1], "$P$", fontsize=FS + 0.5, va="center", ha="center", bbox=halo, zorder=9)
y0 = P[1] + disc + 1*step
for yy in (y0, y0 + step): a.plot([-0.26, -0.06], [yy, yy], color="k", lw=0.5)
a.annotate("", xy=(-0.16, y0), xytext=(-0.16, y0 + step),
           arrowprops=dict(arrowstyle="<->", lw=0.5, mutation_scale=4, shrinkA=0, shrinkB=0))
a.text(-0.3, y0 + step/2, r"$f_z z_n$", fontsize=FSS, ha="right", va="center")
a.text(-0.1, -1.62, "centre\npath", fontsize=FSS, color=CG, ha="right", va="center")
a.set_xlim(-0.8, 1.3); a.set_ylim(-2.05, 1.85)

# ---------------- (c) Newton on the crossing condition -------------------
# exaggerated feed so that the moving centre is visible: centre y = kf * tau,
# edge-ray angle a = a0 - tau (clockwise rotation)
a = ax[2]; rho = 1.0; kf = 0.55; a0 = np.deg2rad(95)
cen = lambda tau: np.array([0.0, kf*tau])
uu = lambda tau: np.array([np.cos(a0 - tau), np.sin(a0 - tau)])
tau2 = np.deg2rad(52); lF = 0.32
Pc = cen(tau2) + (rho + lF)*uu(tau2)
taus = np.linspace(np.deg2rad(5), np.deg2rad(80), 200)
tr = np.array([cen(t) + rho*uu(t) for t in taus])
a.plot(tr[:, 0], tr[:, 1], color="0.55", lw=0.7, ls=":")
a.text(tr[-1, 0] + 0.03, tr[-1, 1] - 0.06, "tip path", fontsize=FSS, color="0.4", va="top")
a.plot([0, 0], [-0.12, kf*np.deg2rad(80) + 0.12], color=CG, lw=0.6, ls="-.")
a.text(-0.07, kf*np.deg2rad(80) + 0.06, "centre\npath", fontsize=FSS, color=CG, ha="right", va="top")
iters = [(np.deg2rad(22), "0", "0.55"), (np.deg2rad(40), "1", "0.3"), (tau2, "2", CT)]
for tau, lab, col in iters:
    c = cen(tau); u = uu(tau); fin = lab == "2"
    a.plot(*c, "s", color=col if fin else "k", ms=2.8 if not fin else 3.2, zorder=6)
    a.plot([c[0], c[0] + 1.62*u[0]], [c[1], c[1] + 1.62*u[1]], color=col,
           lw=1.2 if fin else 0.7, ls="-" if fin else "--")
    a.plot(*(c + rho*u), "o", color=col, ms=3.0 if not fin else 3.6, zorder=6)
    a.text(c[0] + 1.74*u[0], c[1] + 1.74*u[1], lab, fontsize=FS, color=col, ha="center", va="center")
    if not fin:
        foot = c + np.dot(Pc - c, u)*u
        a.plot([Pc[0], foot[0]], [Pc[1], foot[1]], color=CF, lw=0.7, ls=":")
cF, uF = cen(tau2), uu(tau2); tipF = cF + rho*uF
a.plot([tipF[0], Pc[0]], [tipF[1], Pc[1]], color=CF, lw=2.0, solid_capstyle="butt", zorder=7)
a.text((tipF[0] + Pc[0])/2 + 0.07, (tipF[1] + Pc[1])/2 - 0.17, r"$\ell$", color=CF, fontsize=FS + 0.5)
a.plot(*Pc, "*", color=CF, ms=8, mec="k", mew=0.4, zorder=8)
a.text(Pc[0] + 0.22, Pc[1] + 0.02, "$P$", fontsize=FS + 0.5, va="center", ha="center", bbox=halo, zorder=9)
a.annotate("", xy=(0.13, kf*np.deg2rad(50)), xytext=(0.13, kf*np.deg2rad(22)),
           arrowprops=dict(arrowstyle="-|>", lw=0.6, mutation_scale=6))
a.text(0.18, kf*np.deg2rad(33), "$v_f$", fontsize=FS)
hand = [Line2D([], [], color="0.4", lw=0.7, ls="--", marker="s", ms=2.6, mfc="k", mec="k"),
        Line2D([], [], color=CT, lw=1.2, marker="s", ms=2.8),
        Line2D([], [], color=CF, lw=0.7, ls=":")]
a.legend(hand, ["realigned seed (0), Newton step (1)", "converged edge ray (2)", "miss of $P$"],
         loc="lower right", bbox_to_anchor=(1.08, -0.06), fontsize=6.0, handlelength=2.0,
         frameon=False, borderaxespad=0.0)
a.set_xlim(-0.4, 2.1); a.set_ylim(-0.75, 1.95)

fig.tight_layout(w_pad=0.1, rect=(0, 0, 1, 0.92))
for a, t in zip(ax, ["(a) round insert, cutting part", "(b) six candidates at $P$",
                     "(c) Newton on $f(t)=0$"]):
    bb = a.get_position(); fig.text(bb.x0 + 0.005, 0.955, t, fontsize=7.5)
style.save(fig, ROOT, "fig_schematic")
print("ok")
