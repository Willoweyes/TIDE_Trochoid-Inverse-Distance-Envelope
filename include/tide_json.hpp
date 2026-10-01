// ============================================================================
//  tide_json.hpp -- minimal, dependency-free reader for the benchmark JSON
//  configuration files.  Deliberately NOT a general JSON parser: it scans for
//  "key" : number / "key" : "string" pairs, which is all the flat benchmark
//  configs require.  Keeping this in-tree removes the only third-party
//  dependency the executable would otherwise need.
// ============================================================================
#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include "tide_core.hpp"

namespace tide {

inline std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open config: " + path);
    std::ostringstream ss; ss << f.rdbuf();
    return ss.str();
}

// Find "key" and return the numeric literal that follows the next ':'.
inline bool json_number(const std::string& s, const std::string& key,
                        double* out)
{
    const std::string pat = "\"" + key + "\"";
    std::size_t pos = 0;
    while ((pos = s.find(pat, pos)) != std::string::npos) {
        std::size_t c = s.find(':', pos + pat.size());
        if (c == std::string::npos) return false;
        std::size_t i = c + 1;
        while (i < s.size() && (s[i]==' '||s[i]=='\t'||s[i]=='\n'||s[i]=='\r'))
            ++i;
        if (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i]))
                             || s[i]=='-' || s[i]=='+' || s[i]=='.')) {
            *out = std::strtod(s.c_str() + i, nullptr);
            return true;
        }
        pos += pat.size();          // key present but not numeric -> keep looking
    }
    return false;
}

inline bool json_string(const std::string& s, const std::string& key,
                        std::string* out)
{
    const std::string pat = "\"" + key + "\"";
    std::size_t pos = s.find(pat);
    if (pos == std::string::npos) return false;
    std::size_t c = s.find(':', pos + pat.size());
    if (c == std::string::npos) return false;
    std::size_t a = s.find('"', c + 1);
    if (a == std::string::npos) return false;
    std::size_t b = s.find('"', a + 1);
    if (b == std::string::npos) return false;
    *out = s.substr(a + 1, b - a - 1);
    return true;
}

inline Config config_from_json_file(const std::string& path) {
    const std::string s = read_file(path);
    Config c;
    double v;
    if (json_number(s, "vc", &v))     c.vc = v;
    if (json_number(s, "fz", &v))     c.fz = v;
    if (json_number(s, "ap", &v))     c.ap = v;
    if (json_number(s, "Td", &v))     c.Td = v;
    if (json_number(s, "ri", &v))     c.ri = v;
    if (json_number(s, "z", &v))      c.z = static_cast<int>(v);
    if (json_number(s, "eps_r", &v))  c.eps_r = v;
    if (json_number(s, "eps_a", &v))  c.eps_a = v;
    if (json_number(s, "gama_f", &v)) c.gama_f = v;
    if (json_number(s, "gama_p", &v)) c.gama_p = v;
    if (json_number(s, "phi", &v))    c.phi = v;
    if (json_number(s, "x_0", &v))    c.x_0 = v;
    if (json_number(s, "y_0", &v))    c.y_0 = v;
    if (json_number(s, "z_0", &v))    c.z_0 = v;
    if (json_number(s, "grid_base", &v)) c.grid_base = static_cast<int>(v);
    if (json_number(s, "length", &v)) c.length = v;
    if (json_number(s, "ae", &v))     c.ae = v;
    if (json_number(s, "n_pass", &v)) c.n_pass = static_cast<int>(v);
    std::string nm;
    if (json_string(s, "case_name", &nm)) c.case_name = nm;
    return c;
}

} // namespace tide
