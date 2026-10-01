// bench_scaling -- TIDE runtime versus grid, cutting speed and feed (case A1).
// TIDE has no time loop, so the cost must not depend on v_c or f_z.
// Output: CSV on stdout, median of --reps runs, single thread.
#include "tide_core.hpp"
#include <chrono>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace tide;

static Config base_cfg() {
    Config c;                       // literature case A1
    c.vc = 170.0; c.fz = 0.6; c.ap = 0.5;
    c.Td = 10.0;  c.ri = 5.0; c.z = 2;
    c.eps_r = 0.011; c.eps_a = 0.003;
    c.gama_f = 0.6; c.gama_p = 0.0; c.phi = 90.0;
    c.x_0 = 5.0; c.y_0 = -5.0; c.z_0 = 0.0;
    c.grid_base = 200; c.length = 5.0;
    return c;
}

static double timed(const Kinematics& k, int reps, Roughness* out) {
    std::vector<double> ts;
    for (int i = 0; i < reps; ++i) {
        const auto t0 = std::chrono::steady_clock::now();
        Roughness r = label_only(k);
        const auto t1 = std::chrono::steady_clock::now();
        ts.push_back(std::chrono::duration<double>(t1 - t0).count());
        if (out) *out = r;                       // deterministic result
    }
    std::sort(ts.begin(), ts.end());
    return ts[ts.size() / 2];                    // median
}

int main(int argc, char** argv) {
    int reps = 3;
    for (int i = 1; i < argc; ++i)
        if (!std::strcmp(argv[i], "--reps") && i + 1 < argc)
            reps = std::atoi(argv[++i]);

    std::printf("sweep,knob_value,grid_base,n_points,z_n,n_pass,"
                "vc,fz,t_s,points_per_s,Sa,Sq,Sz\n");

    auto emit = [&](const char* sweep, double knob, const Config& c,
                    int) {
        const Kinematics k = kinematics_from_config(c);
        Roughness r;
        const double t = timed(k, reps, &r);
        std::printf("%s,%.10g,%d,%zu,%d,%d,%.10g,%.10g,"
                    "%.9g,%.9g,%.10g,%.10g,%.10g\n",
                    sweep, knob, c.grid_base, r.n_points, c.z, c.n_pass,
                    c.vc, c.fz, t,
                    static_cast<double>(r.n_points) / t, r.Sa, r.Sq, r.Sz);
        std::fflush(stdout);
    };

    // --- grid refinement: N_grid = (2*gb+1)*(gb+1) ~ 2 gb^2 ---
    for (int gb : {25, 35, 50, 71, 100, 141, 200, 283, 400}) {
        Config c = base_cfg(); c.grid_base = gb;
        emit("grid", gb, c, 2);
    }
    // --- cutting speed: cost must be flat (no time loop) ---
    for (double vc : {50.0, 100.0, 170.0, 300.0, 600.0, 1200.0, 2400.0}) {
        Config c = base_cfg(); c.vc = vc; c.grid_base = 200;
        emit("speed", vc, c, 2);
    }
    // --- feed per tooth: cost must be flat ---
    for (double fz : {0.05, 0.1, 0.2, 0.4, 0.6, 0.9, 1.2}) {
        Config c = base_cfg(); c.fz = fz; c.grid_base = 200;
        emit("feed", fz, c, 2);
    }
    return 0;
}
