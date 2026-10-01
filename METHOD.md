# How TIDE computes the height at one point

This note spells out the algorithm step by step. The paper (Sect. 3) gives the
short version; the code is `contact_depth_insert` in `include/tide_core.hpp`
and `contact_offset` in `src/reference_python/tide.py`.

## 1. The real motion

The tool centre moves along `y` at feed speed `v_f`; insert `K` sits at radius
`rho_K` from the centre and turns at `omega`, starting from angle `phi_K`:

    c(t)     = ( x_c , y_0 + v_f t )
    Theta(t) = phi_K - omega t
    tip(t)   = c(t) + rho_K ( cos Theta , -sin Theta )

The two motions together make the tip trace a trochoid.

## 2. The equation to solve (transcendental)

Every point of the cutting edge lies on one ray from the tool centre, turned
by the radial rake `gamma_f` from the tip. Point `P = (X, Y)` is cut when that
ray passes through `P`:

    f(t) = cos a(t) (Y - y_0 - v_f t) + sin a(t) (X - x_c) = 0,
    a(t) = Theta(t) - gamma_f

`t` appears inside `cos`/`sin` and outside them (`v_f t`), so `t` cannot be
isolated: like Kepler's equation, `f(t) = 0` has no closed-form solution.

## 3. Freeze the feed: a circle problem with a closed form

Within one revolution the tool advances only `f_z z_n`, small next to its
radius. Freeze the centre at `(x_c, Y_c)` and ask the tip to pass through `P`:

    (X - x_c)^2 + (Y - Y_c)^2 = rho_K^2

With `h = X - x_c` (sideways distance of `P` from the centre line):

    Y_c     = Y + s sqrt(rho_K^2 - h^2),        s = +1 or -1
    Theta*_s = atan2( s sqrt(rho_K^2 - h^2) , h )

For `|h| >= rho_K` the square root is taken as 0: both signs coincide with the
perpendicular foot of `P` on the tip path.

## 4. From angle to time: why the revolution number is whole

The angle repeats every revolution, so the times that realise `Theta*` are

    t_j = (phi_K - Theta*_s) / omega + j T_rev,   T_rev = 2 pi / omega,  j integer

The centre reaches `Y_c` at `t_c = (Y_c - y_0) / v_f`, which in general is not
one of the `t_j`: the centre being at `Y_c` and the tip being at `Theta*` do not
happen at the same instant. We take the nearest whole revolution,

    j* = round( (t_c - t_0) / T_rev ),   t_0 = (phi_K - Theta*_s) / omega

## 5. What rounding costs

`|t_j* - t_c| <= T_rev / 2`, so the real centre at the seed time lies within

    | y_0 + v_f t_j* - Y_c | <= v_f T_rev / 2 = f_z z_n / 2

of `Y_c`. The seeds are `j* - 1, j*, j* + 1` for both `s`: at most six per
insert. The neighbours are kept because the nearest revolution is not always
the one that cuts deepest. Why these three revolutions bracket the deepest cut
is argued in Sect. 3.3 of the paper (Eq. 10); the argument is checked
numerically near |h| = rho_K, where the two branches merge.

## 6. Realignment: aim again from the real centre

Take the real centre at the seed time, `y_c = y_0 + v_f t`, and recompute the
angle that aims the edge ray at `P`, staying on the same revolution:

    Theta' = atan2( y_c - Y , X - x_c )
    t_b    = (phi_K - gamma_f - Theta') / omega
    t      <- t_b + round( (t - t_b) / T_rev ) T_rev

Still a circle calculation (the centre is held where it is while the ray
turns), but now around the right centre.

## 7. Newton on the real equation

    t <- t - f(t) / f'(t)
    f'(t) = omega sin a (Y - y_0 - v_f t) - v_f cos a - omega cos a (X - x_c)

`f'` carries both the rotation (`omega`) and the feed (`v_f`), so this is the
full trochoidal motion. Near the root `|f'| >= omega |P - c| - v_f > 0`, and
`v_f / (omega |P - c|) = f_z z_n / (2 pi |P - c|)` (about 0.04 for case A1):
the ray turns much faster than the centre moves, so each revolution has exactly
one simple root and `f` is nearly linear around it. Two Newton steps bring
the height to within 2.4e-10 um and three bring the crossing to round-off
(`data/processed/newton_steps.txt`); the code uses four as a margin.

## 8. Height

With the root `t`, the offset along the ray beyond the tip is

    l = cos a (X - x_c) - sin a (Y - y_0 - v_f t) - rho_K

(for non-zero axial rake `gamma_p` the radius is `rho_K + l + sin(gamma_p) z0(|l|)`,
inverted by Newton). Candidates with `|l| > DeltaL` (outside the engaged edge)
are discarded, and

    z(P) = min over inserts K and candidates of
           z_0 + K eps_a + ( R - sqrt(R^2 - l^2) ) cos(gamma_p)

Points reached by no candidate keep the stock height `z_0 + a_p`.

## Verification

Step `verify` of `scripts/run_experiments.py` compares these heights with a
brute-force scan of every revolution for sign changes of `f` (4000 samples per
revolution, 80 bisection steps), on random points and on bands at the edges of
the swept region; the largest difference is 9.1e-13 um
(`data/processed/edge_verification.txt`, Table 2 of the paper).
