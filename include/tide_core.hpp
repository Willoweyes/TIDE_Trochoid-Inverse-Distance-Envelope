// ============================================================================
//  tide_core.hpp -- TIDE (Trochoid Inverse-Distance Envelope), header-only C++.
//
//  C++ implementation of TIDE as described in the paper (Sec. 3) and in
//  METHOD.md; mirrors src/reference_python/tide.py (contact_offset).
//  Per point and insert: six analytic seeds from the circle family,
//  one realignment of the edge ray, four Newton steps on the edge-ray
//  crossing f(t)=0, then the lower envelope over inserts and passes.
//
//  Units: mm, rad, s internally.  Heights are returned in micrometres.
//
//  Determinism.  Every per-point result depends only on that point's
//  coordinates and the kinematic parameters; there is no reduction across
//  points, no atomics, and no random state.  The only floating-point
//  reduction is the `min` over a fixed, statically ordered candidate list.
//  Results are therefore bitwise identical between serial and parallel runs
//  and across repeated runs on the same binary.  -ffast-math is deliberately
//  NOT used, so IEEE-754 semantics are preserved.
//
//  Author: Duong Duc Tri.
// ============================================================================
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>

namespace tide {

constexpr double PI = 3.14159265358979323846;

// ---------------------------------------------------------------------------
//  Kinematics: every quantity derived from a process/geometry configuration.
//  Mirrors the reference dataclass `Kinematics` field for field.
// ---------------------------------------------------------------------------
struct Kinematics {
    double D       = 10.0;   // cutting diameter                        [mm]
    double R       = 5.0;    // insert corner (nose) radius              [mm]
    int    z_n     = 2;      // number of inserts                         [-]
    double f_z     = 0.6;    // feed per tooth                     [mm/tooth]
    double a_p     = 0.5;    // axial depth of cut                       [mm]
    double eps_r   = 0.011;  // radial run-out step between inserts      [mm]
    double eps_a   = 0.003;  // axial  run-out step between inserts      [mm]
    double gamma_f = 0.0;    // radial rake angle                       [rad]
    double gamma_p = 0.0;    // axial  rake angle                       [rad]
    double phi     = PI/2;   // initial phase of insert K=0             [rad]
    double x0      = 5.0;    // tool-centre X                            [mm]
    double y0      = -5.0;   // tool-centre Y at t=0                     [mm]
    double z0      = 0.0;    // workpiece datum                          [mm]
    double a_e     = 0.0;    // step-over between passes                 [mm]
    int    n_pass  = 1;      // number of passes                          [-]
    double omega   = 0.0;    // spindle angular velocity              [rad/s]
    double v_f     = 0.0;    // feed velocity                          [mm/s]
    double t_total = 0.0;    // machining-time window                     [s]
    int    grid_m  = 400;    // X-nodes - 1  (cross-feed)                 [-]
    int    grid_n  = 200;    // Y-nodes - 1  (feed)                       [-]
    double Lx      = 10.0;   // grid extent in X                         [mm]
    double Ly      = 5.0;    // grid extent in Y                         [mm]

    double T_rev() const { return 2.0 * PI / omega; }
};

// Raw configuration in the units used by the benchmark JSON files.
struct Config {
    double vc = 170.0;       // cutting speed                        [m/min]
    double fz = 0.6;         // feed per tooth                    [mm/tooth]
    double ap = 0.5;         // axial depth of cut                      [mm]
    double Td = 10.0;        // cutting diameter                        [mm]
    double ri = 5.0;         // insert corner radius                    [mm]
    int    z  = 2;           // number of inserts                        [-]
    double eps_r = 0.011;    // radial run-out                          [mm]
    double eps_a = 0.003;    // axial  run-out                          [mm]
    double gama_f = 0.6;     // radial rake                            [deg]
    double gama_p = 0.0;     // axial  rake                            [deg]
    double phi    = 90.0;    // initial phase                          [deg]
    double x_0 = 5.0, y_0 = -5.0, z_0 = 0.0;   //                      [mm]
    int    grid_base = 200;  // grid_m = 2*grid_base, grid_n = grid_base [-]
    double length = 5.0;     // Ly; Lx = 2*length                       [mm]
    double ae = 0.0;         // step-over                               [mm]
    int    n_pass = 1;       // passes                                   [-]
    std::string case_name = "case";
};

inline double deg2rad(double d) { return d * PI / 180.0; }

// ---------------------------------------------------------------------------
//  Derive kinematics.  Reference relations:
//      n_s   = vc / (pi D)          [rev/min], D in metres
//      omega = 2 pi n_s / 60        [rad/s]
//      v_f   = f_z z_n n_s / 60     [mm/s]   =>  v_f * T_rev = f_z z_n
//      t_total = (1.5 D + Ly) / v_f
// ---------------------------------------------------------------------------
inline Kinematics kinematics_from_config(const Config& c) {
    Kinematics k;
    k.D = c.Td; k.R = c.ri; k.z_n = c.z; k.f_z = c.fz; k.a_p = c.ap;
    k.eps_r = c.eps_r; k.eps_a = c.eps_a;
    k.gamma_f = deg2rad(c.gama_f); k.gamma_p = deg2rad(c.gama_p);
    k.phi = deg2rad(c.phi);
    k.x0 = c.x_0; k.y0 = c.y_0; k.z0 = c.z_0;
    k.a_e = c.ae; k.n_pass = c.n_pass;

    const double n_s = c.vc / (c.Td * 1e-3 * PI);   // [rpm]
    k.omega = n_s * PI / 30.0;                      // [rad/s]
    k.v_f   = c.fz * static_cast<double>(c.z) * n_s / 60.0;

    k.grid_m = 2 * c.grid_base;
    k.grid_n = c.grid_base;
    k.Lx = 2.0 * c.length;
    k.Ly = c.length;
    k.t_total = (1.5 * c.Td + k.Ly) / k.v_f;
    return k;
}

// ---------------------------------------------------------------------------
//  Depth profile of the round-nose edge, and the engaged half-width.
//      zeta(l)  = (R - sqrt(R^2 - l^2)) cos(gamma_p),   0 <= l <= R
//      DeltaL   = max( sqrt(R^2-(R-a_p)^2),  f_z/(2 cos gamma_f) )
// ---------------------------------------------------------------------------
inline double zeta(double l, double R, double cos_gp) {
    const double r2 = R * R - l * l;
    return (R - (r2 > 0.0 ? std::sqrt(r2) : 0.0)) * cos_gp;
}

inline double engaged_half_width(const Kinematics& k) {
    const double s = k.R * k.R - (k.R - k.a_p) * (k.R - k.a_p);
    const double t1 = std::sqrt(s > 0.0 ? s : 0.0);
    const double t2 = k.f_z / (2.0 * std::cos(k.gamma_f));
    return std::max(t1, t2);
}

// ---------------------------------------------------------------------------
//  Exact contact of the edge.
//
//  In the reference geometry every edge point of insert K lies on one ray
//  from the tool centre c(t), direction u = (cos a, -sin a),
//  a = Theta(t) - gamma_f, at radius rho + l.  P is cut when the ray passes
//  through it:  f(t) = u x (P - c) = 0,  u.(P - c) > 0,  and then
//  l = u.(P - c) - rho exactly.  Each of the six tip candidates (after
//  realignment) seeds `cross_steps` Newton steps on f.
//  Returns the smallest depth zeta(|l|) + K eps_a over the candidates whose
//  |l| <= dL, or +inf.  Mirrors tide.contact_offset in the Python reference.
// ---------------------------------------------------------------------------
inline double contact_depth_insert(double Px, double Py, const Kinematics& k,
                                   int K, double x_c, double dL,
                                   int cross_steps = 4)
{
    const double rho  = k.D / 2.0 + static_cast<double>(K) * k.eps_r;
    const double phi0 = k.phi + 2.0 * PI * static_cast<double>(K)
                                / static_cast<double>(k.z_n);
    const double omega = k.omega, v_f = k.v_f, y0 = k.y0;
    const double t_total = k.t_total, Trev = k.T_rev();
    const double gf = k.gamma_f, sp = std::sin(k.gamma_p);
    const double cos_gp = std::cos(k.gamma_p);
    const double h = Px - x_c;
    const double val = rho * rho - h * h;
    const double disc = (val > 0.0) ? std::sqrt(val) : 0.0;

    double best = std::numeric_limits<double>::infinity();
    for (int si = 0; si < 2; ++si) {
        const double s = (si == 0) ? 1.0 : -1.0;
        const double Yc = Py + s * disc;
        const double Theta_star = std::atan2(s * disc, h);
        const double t_base = (phi0 - Theta_star) / omega;
        const double j = std::nearbyint(((Yc - y0) / v_f - t_base) / Trev);
        for (int dj = -1; dj <= 1; ++dj) {
            double t = t_base + (j + static_cast<double>(dj)) * Trev;
            t = std::min(std::max(t, 0.0), t_total);
            {   // realign
                const double yc = y0 + v_f * t;
                const double Ths = std::atan2(yc - Py, Px - x_c);
                const double tb  = (phi0 - gf - Ths) / omega;   // aim the edge ray at P
                t = tb + std::nearbyint((t - tb) / Trev) * Trev;
                t = std::min(std::max(t, 0.0), t_total);
            }
            for (int it = 0; it < cross_steps; ++it) {    // Newton on f
                const double a = phi0 - omega * t - gf;
                const double ca = std::cos(a), sa = std::sin(a);
                const double dx = Px - x_c, dy = Py - (y0 + v_f * t);
                const double f  = ca * dy + sa * dx;
                const double fp = omega * sa * dy - ca * v_f - omega * ca * dx;
                if (std::fabs(fp) > 1e-30) t -= f / fp;
            }
            if (!(t >= 0.0 && t <= t_total)) continue;
            const double a = phi0 - omega * t - gf;
            const double dx = Px - x_c, dy = Py - (y0 + v_f * t);
            const double r = std::cos(a) * dx - std::sin(a) * dy;
            if (!(r > 0.0)) continue;
            double l = r - rho;
            if (std::fabs(sp) > 1e-14) {                  // r = rho + l + sp z0(|l|)
                for (int it = 0; it < 4; ++it) {
                    const double al = std::min(std::fabs(l), k.R);
                    const double rt = std::sqrt(std::max(k.R * k.R - al * al, 0.0));
                    const double F  = rho + l + sp * (k.R - rt) - r;
                    const double sg = (l > 0.0) - (l < 0.0);
                    const double dF = 1.0 + sp * sg * al / (rt > 1e-300 ? rt : 1.0);
                    l -= F / dF;
                }
            }
            const double al = std::fabs(l);
            if (al <= dL) {
                const double depth = zeta(std::min(al, k.R), k.R, cos_gp)
                                   + static_cast<double>(K) * k.eps_a;
                if (depth < best) best = depth;
            }
        }
    }
    return best;
}

// ---------------------------------------------------------------------------
//  Height at a single query point: lower envelope over passes and inserts,
//      z(P) = min_{p,K} [ z0 + K eps_a + zeta(|l|) ],
//  with l from the exact contact above.  Returned in micrometres.
// ---------------------------------------------------------------------------
inline double height_at(double Px, double Py, const Kinematics& k,
                        double dL, bool* cut = nullptr)
{
    double best = std::numeric_limits<double>::infinity();
    for (int p = 0; p < k.n_pass; ++p) {
        const double x_c = k.x0 + static_cast<double>(p) * k.a_e;
        for (int K = 0; K < k.z_n; ++K) {
            const double d = contact_depth_insert(Px, Py, k, K, x_c, dL);
            if (d < best) best = d;
        }
    }
    const bool uncut = !(best < std::numeric_limits<double>::infinity());
    if (cut) *cut = !uncut;
    if (uncut) best = k.a_p;
    return (k.z0 + best) * 1000.0;
}

// ---------------------------------------------------------------------------
//  Full surface on the (grid_m+1) x (grid_n+1) node grid, row-major in X.
//  Index layout matches the reference: idx = i*(grid_n+1) + jj.
// ---------------------------------------------------------------------------
struct SurfaceResult {
    std::vector<double> z_um;
    double uncut_frac = 0.0;
    std::size_t n_points = 0;
};

inline SurfaceResult surface(const Kinematics& k) {
    const int MX = k.grid_m + 1, NY = k.grid_n + 1;
    const double dx = k.Lx / k.grid_m, dy = k.Ly / k.grid_n;
    const double dL = engaged_half_width(k);

    SurfaceResult out;
    out.n_points = static_cast<std::size_t>(MX) * static_cast<std::size_t>(NY);
    out.z_um.resize(out.n_points);

    std::int64_t n_uncut = 0;
#ifdef TIDE_OPENMP
    #pragma omp parallel for schedule(static) reduction(+:n_uncut)
#endif
    for (int i = 0; i < MX; ++i) {
        const double Px = static_cast<double>(i) * dx;
        for (int jj = 0; jj < NY; ++jj) {
            const double Py = static_cast<double>(jj) * dy;
            bool cut = false;
            out.z_um[static_cast<std::size_t>(i) * NY + jj] =
                height_at(Px, Py, k, dL, &cut);
            if (!cut) ++n_uncut;
        }
    }
    out.uncut_frac = static_cast<double>(n_uncut)
                   / static_cast<double>(out.n_points);
    return out;
}

// ---------------------------------------------------------------------------
//  ISO 25178-2 areal parameters.  Heights in micrometres.
// ---------------------------------------------------------------------------
struct Roughness {
    double Sa = 0, Sq = 0, Sz = 0, Sp = 0, Sv = 0, Ssk = 0, Sku = 0;
    std::size_t n_points = 0;
};

inline Roughness compute_roughness(const std::vector<double>& z_um) {
    Roughness r;
    const std::size_t n = z_um.size();
    r.n_points = n;
    if (n == 0) return r;

    // Two-pass mean (numerically safer than a single accumulating pass).
    double mean = 0.0;
    for (double v : z_um) mean += v;
    mean /= static_cast<double>(n);

    double s_abs = 0.0, s_sq = 0.0, s_c3 = 0.0, s_c4 = 0.0;
    double zmin = z_um[0], zmax = z_um[0];
    for (double v : z_um) {
        const double c = v - mean;
        s_abs += std::fabs(c);
        s_sq  += c * c;
        s_c3  += c * c * c;
        s_c4  += c * c * c * c;
        if (v < zmin) zmin = v;
        if (v > zmax) zmax = v;
    }
    const double N = static_cast<double>(n);
    r.Sa = s_abs / N;
    r.Sq = std::sqrt(s_sq / N);
    r.Sz = zmax - zmin;
    r.Sp = zmax - mean;
    r.Sv = mean - zmin;
    r.Ssk = (r.Sq > 1e-12) ? (s_c3 / N) / (r.Sq * r.Sq * r.Sq) : 0.0;
    r.Sku = (r.Sq > 1e-12) ? (s_c4 / N) / (r.Sq * r.Sq * r.Sq * r.Sq) : 0.0;
    return r;
}

// ---------------------------------------------------------------------------
//  Label generation: roughness for one parameter vector, no surface stored.
//  This is the operation a surrogate-training pipeline actually calls.
//  Heights are accumulated on the fly, so memory is O(1) in the grid size.
// ---------------------------------------------------------------------------
inline Roughness label_only(const Kinematics& k) {
    const int MX = k.grid_m + 1, NY = k.grid_n + 1;
    const double dx = k.Lx / k.grid_m, dy = k.Ly / k.grid_n;
    const double dL = engaged_half_width(k);
    const std::size_t N = static_cast<std::size_t>(MX)
                        * static_cast<std::size_t>(NY);

    // Deterministic reduction.  A plain OpenMP `reduction(+:...)` makes the
    // summation ORDER depend on the thread count, which perturbs Sa and Sq in
    // the last few bits (~1e-15 relative) and breaks bit-for-bit
    // reproducibility across machines and thread settings.  Reproducibility is
    // a stated requirement for a label generator, so instead each column i
    // accumulates into its own slot and the slots are combined afterwards in a
    // fixed ascending order.  The result is therefore independent of the
    // number of threads and of the scheduling, while the parallel speed-up is
    // preserved.  (Sz, being a min/max reduction, is exact in any order.)
    std::vector<double> col_sum(static_cast<std::size_t>(MX), 0.0);
    std::vector<double> col_min(static_cast<std::size_t>(MX),
                                std::numeric_limits<double>::infinity());
    std::vector<double> col_max(static_cast<std::size_t>(MX),
                               -std::numeric_limits<double>::infinity());

#ifdef TIDE_OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < MX; ++i) {
        const double Px = static_cast<double>(i) * dx;
        double s = 0.0;
        double lo =  std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (int jj = 0; jj < NY; ++jj) {
            const double z = height_at(Px, static_cast<double>(jj) * dy,
                                       k, dL, nullptr);
            s += z;
            if (z < lo) lo = z;
            if (z > hi) hi = z;
        }
        col_sum[static_cast<std::size_t>(i)] = s;
        col_min[static_cast<std::size_t>(i)] = lo;
        col_max[static_cast<std::size_t>(i)] = hi;
    }

    double sum = 0.0;
    double zmin =  std::numeric_limits<double>::infinity();
    double zmax = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < col_sum.size(); ++i) {
        sum += col_sum[i];
        if (col_min[i] < zmin) zmin = col_min[i];
        if (col_max[i] > zmax) zmax = col_max[i];
    }
    const double mean = sum / static_cast<double>(N);

    // Pass 2: central moments, same fixed-order combination.
    std::vector<double> c_abs(static_cast<std::size_t>(MX), 0.0);
    std::vector<double> c_sq (static_cast<std::size_t>(MX), 0.0);
    std::vector<double> c_c3 (static_cast<std::size_t>(MX), 0.0);
    std::vector<double> c_c4 (static_cast<std::size_t>(MX), 0.0);

#ifdef TIDE_OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < MX; ++i) {
        const double Px = static_cast<double>(i) * dx;
        double a = 0.0, q = 0.0, t3 = 0.0, t4 = 0.0;
        for (int jj = 0; jj < NY; ++jj) {
            const double c = height_at(Px, static_cast<double>(jj) * dy,
                                       k, dL, nullptr) - mean;
            a  += std::fabs(c);
            q  += c * c;
            t3 += c * c * c;
            t4 += c * c * c * c;
        }
        c_abs[static_cast<std::size_t>(i)] = a;
        c_sq [static_cast<std::size_t>(i)] = q;
        c_c3 [static_cast<std::size_t>(i)] = t3;
        c_c4 [static_cast<std::size_t>(i)] = t4;
    }

    double s_abs = 0.0, s_sq = 0.0, s_c3 = 0.0, s_c4 = 0.0;
    for (std::size_t i = 0; i < c_abs.size(); ++i) {
        s_abs += c_abs[i]; s_sq += c_sq[i];
        s_c3  += c_c3[i];  s_c4 += c_c4[i];
    }

    Roughness r;
    r.n_points = N;
    r.Sa = s_abs / static_cast<double>(N);
    r.Sq = std::sqrt(s_sq / static_cast<double>(N));
    r.Sz = zmax - zmin;
    r.Sp = zmax - mean;
    r.Sv = mean - zmin;
    r.Ssk = (r.Sq > 1e-12) ? (s_c3 / static_cast<double>(N))
                             / (r.Sq * r.Sq * r.Sq) : 0.0;
    r.Sku = (r.Sq > 1e-12) ? (s_c4 / static_cast<double>(N))
                             / (r.Sq * r.Sq * r.Sq * r.Sq) : 0.0;
    return r;
}

} // namespace tide
