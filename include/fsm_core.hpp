// ============================================================================
//  fsm_core.hpp -- C++ port of the forward simulation method (Z-map) of
//  Bakhshan et al., following src/reference_efsm/efsm_core.py::run_surface
//  (circular insert) operation for operation: Z initialised to a_p, the tool
//  stepped through t = i*dt, every edge point stamped into the cell it lands
//  in (indices truncated and clamped to the border), minimum kept.
// ============================================================================
#pragma once
#include <cmath>
#include <cstdint>
#include <vector>
#include <algorithm>
#include "tide_core.hpp"

namespace tide {

struct FsmResult { std::vector<double> z_um; int m = 0, n = 0; long long n_time_steps = 0;
                   int edge_points = 0; double t_s = 0.0; };

inline FsmResult fsm_surface(const Config& c, double delta_t, int edge_points = 0) {
    const double ws = c.vc / (c.Td * 0.001 * PI);
    const double vf = c.fz * c.z * ws / 60.0;
    const double omega = ws * PI / 30.0;
    const int m = 2 * c.grid_base, n = c.grid_base;
    const double Lx = 2.0 * c.length, Ly = c.length;
    const int grid_t = edge_points > 0 ? edge_points : std::max(m, n) * 2;   // released code: 4n
    const double t_total = (1.5 * c.Td + Ly) / vf;
    const int P = grid_t;
    const long long Nt = static_cast<long long>(t_total / delta_t) + 1;
    const double dx = Lx / m, dy = Ly / n, inv_dx = 1.0 / dx, inv_dy = 1.0 / dy;
    const double gf = deg2rad(c.gama_f), gp = deg2rad(c.gama_p), phi = deg2rad(c.phi);
    const double cos_gp = std::cos(gp), sin_gp = std::sin(gp);
    const double R = c.ri;
    const double Lap = std::sqrt(std::max(0.0, R * R - (R - c.ap) * (R - c.ap)));
    const double Lfz = c.fz / (2.0 * std::cos(deg2rad(c.gama_f)));
    const double L = std::max(Lap, Lfz);
    std::vector<double> l(P), zl(P), wl(P);
    for (int i = 0; i < P; ++i) {
        l[i] = (P == 1) ? -L : -L + (2.0 * L) * i / (P - 1);
        zl[i] = R - std::sqrt(std::max(R * R - l[i] * l[i], 0.0));
        wl[i] = l[i] + sin_gp * zl[i];
    }
    FsmResult out; out.m = m; out.n = n; out.n_time_steps = Nt; out.edge_points = P;
    std::vector<double> Z(static_cast<std::size_t>(m + 1) * (n + 1), c.ap);
    std::vector<double> zp(P);
    for (long long it = 0; it < Nt; ++it) {
        const double t = static_cast<double>(it) * delta_t;
        for (int K = 1; K <= c.z; ++K) {
            const double dr = c.Td / 2.0 + (K - 1) * c.eps_r;
            const double da = (K - 1) * c.eps_a;
            const double var1 = phi + 2.0 * PI * (K - 1) / c.z - omega * t;
            const double cd = std::cos(gf - var1), sd = std::sin(gf - var1);
            const double cx = sd * (da * sin_gp) + c.x_0;
            const double cy = cd * (da * sin_gp) + c.y_0 + vf * t;
            const double zoff = c.z_0 + da * cos_gp;
            for (int pi = 0; pi < P; ++pi) {
                const double w2 = wl[pi] + dr;
                const double xp = cd * w2 + cx, yp = sd * w2 + cy;
                long long xi = static_cast<long long>(xp * inv_dx);
                long long yi = static_cast<long long>(yp * inv_dy);
                xi = std::min<long long>(std::max<long long>(xi, 0), m);
                yi = std::min<long long>(std::max<long long>(yi, 0), n);
                const double z = cos_gp * zl[pi] + zoff;
                double& cell = Z[static_cast<std::size_t>(xi) * (n + 1) + yi];
                if (z < cell) cell = z;
            }
        }
    }
    out.z_um.resize(Z.size());
    for (std::size_t i = 0; i < Z.size(); ++i) out.z_um[i] = Z[i] * 1000.0;
    return out;
}

} // namespace tide
