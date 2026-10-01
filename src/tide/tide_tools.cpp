// tide_tools -- one C++ driver for every experiment in the paper.
//   tide_tools fsm    --config F --grid-base N --dt DT [--edge NE] [--dump-bin FILE]
//   tide_tools tide   --config F --grid-base N          [--dump-bin FILE]
//   tide_tools points --config F --in PTS.bin --out OUT.bin
// fsm/tide print one JSON line (Sa, Sq, Sz, n_points, t_s, ...).
// points reads (x,y) float64 pairs and writes one TIDE height [um] per point.
#include "tide_core.hpp"
#include "tide_json.hpp"
#include "fsm_core.hpp"
#include <chrono>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
using namespace tide;

static void write_bin(const std::string& f, const std::vector<double>& v) {
    std::ofstream o(f, std::ios::binary); o.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(double));
}
static std::vector<double> read_bin(const std::string& f) {
    std::ifstream i(f, std::ios::binary | std::ios::ate); std::size_t n = i.tellg() / sizeof(double);
    i.seekg(0); std::vector<double> v(n); i.read(reinterpret_cast<char*>(v.data()), n * sizeof(double)); return v;
}
static void print_r(const char* what, const Roughness& r, std::size_t n, double t, long long nt) {
    std::printf("{\"method\":\"%s\",\"Sa\":%.15g,\"Sq\":%.15g,\"Sz\":%.15g,\"n_points\":%zu,\"t_s\":%.9g,\"n_time_steps\":%lld}\n",
                what, r.Sa, r.Sq, r.Sz, n, t, nt);
}
int main(int argc, char** argv) {
    if (argc < 2) return 2;
    std::string mode = argv[1], cfgf, dump, fin, fout; int gb = -1, edge = 0; double dt = 8e-7;
    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--config") cfgf = argv[++i]; else if (a == "--grid-base") gb = std::atoi(argv[++i]);
        else if (a == "--dt") dt = std::atof(argv[++i]); else if (a == "--dump-bin") dump = argv[++i];
        else if (a == "--in") fin = argv[++i]; else if (a == "--edge") edge = std::atoi(argv[++i]); else if (a == "--out") fout = argv[++i];
    }
    Config c = config_from_json_file(cfgf); if (gb > 0) c.grid_base = gb;
    const auto t0 = std::chrono::steady_clock::now();
    if (mode == "fsm") {
        FsmResult f = fsm_surface(c, dt, edge);
        const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        Roughness r = compute_roughness(f.z_um); if (!dump.empty()) write_bin(dump, f.z_um);
        print_r("fsm", r, f.z_um.size(), t, f.n_time_steps);
    } else if (mode == "tide") {
        const Kinematics k = kinematics_from_config(c); SurfaceResult s = surface(k);
        const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        Roughness r = compute_roughness(s.z_um); if (!dump.empty()) write_bin(dump, s.z_um);
        print_r("tide", r, s.z_um.size(), t, 0);
    } else if (mode == "points") {
        const Kinematics k = kinematics_from_config(c); const double dL = engaged_half_width(k);
        std::vector<double> p = read_bin(fin), o; const std::size_t N = p.size() / 2; o.reserve(N);
        for (std::size_t i = 0; i < N; ++i) o.push_back(height_at(p[2 * i], p[2 * i + 1], k, dL));
        write_bin(fout, o);
    }
    return 0;
}
