// Thin driver around the UNMODIFIED simulation.cpp / simulation.hpp of surf-topo
// (github.com/HadiBakhshan/surf-topo, commit 1f5d3e0).  It reproduces
// scripts/data_set_builder.py::compute_derived + run_surface_simulation for one case file.
//   bakhshan_gcc --config F --grid-base N --dt DT [--edge NE] [--dump-bin FILE]
// Output: one JSON line (Sa, Sq, Sz in um with the estimators of tide_core.hpp; t_s = wall clock of compute_surface).
#include "simulation.hpp"
#include "tide_json.hpp"      // same case-file reader as tide_tools (inputs identical by construction)
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <omp.h>

int main(int argc, char** argv) {
    std::string cfgf, dump; int gb = -1, edge = 0; double dt = 8e-7;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--config") cfgf = argv[++i]; else if (a == "--grid-base") gb = std::atoi(argv[++i]);
        else if (a == "--dt") dt = std::atof(argv[++i]); else if (a == "--edge") edge = std::atoi(argv[++i]);
        else if (a == "--dump-bin") dump = argv[++i];
    }
    tide::Config c = tide::config_from_json_file(cfgf); if (gb > 0) c.grid_base = gb;
    // ---- compute_derived (data_set_builder.py)
    const double PI = 3.14159265358979323846;
    const double ws = c.vc / (c.Td * 0.001 * PI), vf = c.fz * c.z * ws / 60.0, omega = ws * PI / 30.0;
    const int grid_m = c.grid_base * 2, grid_n = c.grid_base;
    const double Lx = c.length * 2, Ly = c.length;
    const int grid_t = std::max(grid_m, grid_n) * 2, edge_pts = edge > 0 ? edge : grid_t;
    const double total_time = (1.5 * c.Td + Ly) / vf;
    Tool tool(c.gama_f, c.gama_p, c.Td, c.eps_r, c.eps_a, c.phi, omega, vf, c.x_0, c.y_0, c.z_0, c.fz, c.ri, c.z);
    Simulation sim(tool, Lx, grid_m, Ly, grid_n, c.ap, edge_pts, total_time, dt);
    std::vector<Eigen::Vector3d> pts;
    const auto t0 = std::chrono::steady_clock::now();
    sim.compute_surface(pts);
    const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::vector<double> z(pts.size());
    for (size_t i = 0; i < pts.size(); ++i) z[i] = pts[i][2] * 1000.0;      // mm -> um
    const tide::Roughness r = tide::compute_roughness(z);
    if (!dump.empty()) { std::ofstream o(dump, std::ios::binary); o.write(reinterpret_cast<const char*>(z.data()), z.size() * sizeof(double)); }
    std::printf("{\"method\":\"bakhshan_gcc\",\"Sa\":%.15g,\"Sq\":%.15g,\"Sz\":%.15g,\"n_points\":%zu,\"t_s\":%.9g,\"n_time_steps\":%lld,\"omp_threads\":%d}\n",
                r.Sa, r.Sq, r.Sz, z.size(), t, (long long)(total_time / dt) + 1, omp_get_max_threads());
    return 0;
}
