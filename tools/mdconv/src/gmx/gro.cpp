#include "mdconv/gmx/gro.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace mdconv {

static bool parse_double(const std::string& s, double& out)
{
    std::istringstream in(s);
    in >> out;
    return in && in.eof();
}

GroFrame parse_gro_from_string(const std::string& text)
{
    std::istringstream in(text);
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        lines.push_back(line);
    }
    while (!lines.empty() && lines.back().empty()) {
        lines.pop_back();
    }

    if (lines.size() < 3) {
        throw std::runtime_error("invalid .gro file");
    }

    std::istringstream count_line(lines[1]);
    int n_atoms = 0;
    count_line >> n_atoms;
    if (!count_line || n_atoms < 0) {
        throw std::runtime_error("invalid .gro atom count");
    }

    if (lines.size() != static_cast<size_t>(n_atoms) + 3) {
        throw std::runtime_error("atom count does not match .gro contents");
    }

    GroFrame frame;
    frame.x_nm.reserve(static_cast<size_t>(n_atoms));
    for (int i = 0; i < n_atoms; ++i) {
        const std::string& atom_line = lines[2 + static_cast<size_t>(i)];
        std::istringstream atom_stream(atom_line);
        std::vector<std::string> tokens;
        std::string tok;
        while (atom_stream >> tok) {
            tokens.push_back(tok);
        }
        if (tokens.size() < 3) {
            throw std::runtime_error("invalid .gro atom line");
        }
        const size_t n = tokens.size();
        Vec3 v;
        if (!parse_double(tokens[n - 3], v.x_nm)
            || !parse_double(tokens[n - 2], v.y_nm)
            || !parse_double(tokens[n - 1], v.z_nm)) {
            throw std::runtime_error("invalid .gro atom line");
        }
        frame.x_nm.push_back(v);
    }

    std::istringstream box_line(lines[2 + static_cast<size_t>(n_atoms)]);
    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    if (!(box_line >> lx >> ly >> lz)) {
        throw std::runtime_error("invalid .gro box line");
    }
    double extra = 0.0;
    if (box_line >> extra) {
        throw std::runtime_error("invalid .gro box line");
    }
    frame.box_nm.lx_nm = lx;
    frame.box_nm.ly_nm = ly;
    frame.box_nm.lz_nm = lz;
    return frame;
}

GroFrame parse_gro_from_file(const std::string& path)
{
    std::ifstream in(path.c_str(), std::ios::in);
    if (!in) {
        throw std::runtime_error("failed to open .gro file: " + path);
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return parse_gro_from_string(buffer.str());
}

} // namespace mdconv
