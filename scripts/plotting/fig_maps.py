"""Fig. 3 -- where the forward error arises, case A1 at 100 cells:
(a) TIDE heights over a 4 x 3 mm part of the window, (b) difference FSM - TIDE
at the same nodes (time step tied to the grid), (c) profile across the feed at
y = 2.5 mm.  LNCS text width."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import style
plt = style.apply()
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
ROOT = Path(__file__).resolve().parents[2]
S = np.load(ROOT/"data/raw/surfaces/A1_fsm_vs_tide.npz")
G = 100; xr, yr = (3.0, 7.0), (0.5, 3.5); YP = 2.5
X, Y = S[f"tide{G}_X"], S[f"tide{G}_Y"]
m = (X[:, 0] >= xr[0] - 1e-9) & (X[:, 0] <= xr[1] + 1e-9)
n = (Y[0] >= yr[0] - 1e-9) & (Y[0] <= yr[1] + 1e-9)
Zt = S[f"tide{G}_Z"][m][:, n]; Zf = S[f"fsm{G}_Z"][m][:, n]; D = Zf - Zt
dx = X[1, 0] - X[0, 0]; dy = Y[0, 1] - Y[0, 0]
ext = (xr[0] - dx/2, xr[1] + dx/2, yr[0] - dy/2, yr[1] + dy/2)
print("window nodes", Zt.shape, "diff min/max %.3f %.3f" % (D.min(), D.max()))

# profile: TIDE at 400 cells (smooth reference curve), FSM at 100 cells
Xt4, Zt4 = S["tide400_X"], S["tide400_Z"]; jt = int(np.argmin(np.abs(S["tide400_Y"][0] - YP)))
Xf, Zfp = S["fsm100_X"], S["fsm100_Z"]; jf = int(np.argmin(np.abs(S["fsm100_Y"][0] - YP)))
assert abs(S["tide400_Y"][0, jt] - YP) < 1e-9 and abs(S["fsm100_Y"][0, jf] - YP) < 1e-9
xs, pt = Xt4[:, jt], Zt4[:, jt]; xf, pf = Xf[:, jf], Zfp[:, jf]
st = (xs >= xr[0]) & (xs <= xr[1]); sf = (xf >= xr[0]) & (xf <= xr[1])

fig = plt.figure(figsize=(style.W_FULL, 1.58))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.0, 1.22],
              left=0.065, right=0.99, bottom=0.33, top=0.90, wspace=0.30)
a0 = fig.add_subplot(gs[0, 0]); a1 = fig.add_subplot(gs[0, 1], sharey=a0); a2 = fig.add_subplot(gs[0, 2])
kw = dict(origin="lower", extent=ext, aspect="equal", interpolation="nearest")
im0 = a0.imshow(Zt.T, cmap="viridis", vmin=0, vmax=10, **kw)
lim = 3.5
im1 = a1.imshow(D.T, cmap="RdBu", norm=TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim), **kw)
for a in (a0, a1):
    a.axhline(YP, color="k", lw=0.6, ls=(0, (3, 2)))
    a.set_xlabel("$x$ [mm]", labelpad=1); a.grid(False)
    a.set_xticks([3, 4, 5, 6, 7]); a.tick_params(labelsize=6.5, pad=1.5)
a0.set_ylabel("$y$ [mm]", labelpad=1); a0.set_yticks([1, 2, 3])
plt.setp(a1.get_yticklabels(), visible=False)
a0.set_title("(a) TIDE heights", loc="left", fontsize=7.5, pad=3)
a1.set_title("(b) FSM $-$ TIDE", loc="left", fontsize=7.5, pad=3)
fig.canvas.draw()
q0, q1 = a0.get_position(), a1.get_position()
c0 = fig.add_axes([q0.x0 + 0.01, q0.y0 - 0.215, q0.width - 0.02, 0.038])
c1 = fig.add_axes([q1.x0 + 0.01, q1.y0 - 0.215, q1.width - 0.02, 0.038])
cb = fig.colorbar(im0, cax=c0, orientation="horizontal", ticks=[0, 5, 10])
cb.set_label("height [µm]", fontsize=6.5, labelpad=1); cb.ax.tick_params(labelsize=6.2, pad=1)
cb2 = fig.colorbar(im1, cax=c1, orientation="horizontal", ticks=[-3, 0, 3])
cb2.set_label("difference [µm]", fontsize=6.5, labelpad=1); cb2.ax.tick_params(labelsize=6.2, pad=1)

a2.fill_between(xs[st], 0, pt[st], color=style.C_EXACT, alpha=0.10, lw=0)
a2.plot(xs[st], pt[st], lw=1.3, color=style.C_EXACT, label="TIDE")
a2.step(xf[sf], pf[sf], where="mid", lw=0.9, color=style.C_FORWARD, label="FSM")
a2.set_xlim(*xr); a2.set_ylim(0, 12.5); a2.set_yticks([0, 4, 8, 12]); a2.set_xticks([3, 4, 5, 6, 7])
a2.set_xlabel("$x$ [mm]", labelpad=1); a2.set_ylabel("height [µm]", labelpad=1)
a2.tick_params(labelsize=6.5, pad=1.5, which="both", top=True, right=True)
a2.set_title("(c) profile at $y$ = 2.5 mm", loc="left", fontsize=7.5, pad=3)
a2.legend(loc="upper center", fontsize=6.2, ncol=2, columnspacing=1.0, handlelength=1.5,
          borderpad=0.3, bbox_to_anchor=(0.5, 1.0))
# match the profile panel height to the maps
fig.canvas.draw()
p0 = a0.get_position(); p2 = a2.get_position()
a2.set_position([p2.x0, p0.y0, p2.width, p0.height])
style.save(fig, ROOT, "fig_maps")
print("ok")
