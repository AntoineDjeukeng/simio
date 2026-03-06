// mini_gro_topology.cpp
#include <array>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

static std::string trim(const std::string& s) {
    size_t b = 0, e = s.size();
    while (b < e && std::isspace((unsigned char)s[b])) ++b;
    while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
    return s.substr(b, e - b);
}

static std::string upper(std::string s) {
    for (char& c : s) c = (char)std::toupper((unsigned char)c);
    return s;
}

static std::unordered_set<std::string> parse_names(const std::string& s) {
    std::unordered_set<std::string> out;
    std::string cur;
    auto flush = [&]() {
        if (!cur.empty()) out.insert(upper(trim(cur)));
        cur.clear();
    };
    for (char c : s) {
        if (c == ',' || std::isspace((unsigned char)c)) flush();
        else cur.push_back(c);
    }
    flush();
    return out;
}

static std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c == '\\' || c == '"') {
            out.push_back('\\');
            out.push_back(c);
        } else if (c == '\n') {
            out += "\\n";
        } else if (c == '\r') {
            out += "\\r";
        } else if (c == '\t') {
            out += "\\t";
        } else {
            out.push_back(c);
        }
    }
    return out;
}

static void write_json_output(
    const std::string& path,
    const std::array<double,6>& channel,
    const std::array<double,6>& box_local,
    const std::array<double,3>& gro_box,
    int natoms_filtered,
    const std::vector<std::tuple<std::string, int, int>>& molecules_by_type_ordered) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open json output file: " + path);

    out << std::fixed << std::setprecision(6);
    out << "{\n";
    out << "  \"channel\": {\n";
    out << "    \"min\": [" << channel[0] << ", " << channel[1] << ", " << channel[2] << "],\n";
    out << "    \"max\": [" << channel[3] << ", " << channel[4] << ", " << channel[5] << "]\n";
    out << "  },\n";
    out << "  \"box_local\": {\n";
    out << "    \"min\": [" << box_local[0] << ", " << box_local[1] << ", " << box_local[2] << "],\n";
    out << "    \"max\": [" << box_local[3] << ", " << box_local[4] << ", " << box_local[5] << "]\n";
    out << "  },\n";
    out << "  \"gro_box\": [" << gro_box[0] << ", " << gro_box[1] << ", " << gro_box[2] << "],\n";
    out << "  \"natoms_filtered\": " << natoms_filtered << ",\n";
    out << "  \"types\": [\n";
    for (size_t i = 0; i < molecules_by_type_ordered.size(); ++i) {
        const auto& kv = molecules_by_type_ordered[i];
        out << "    {\"name\": \"" << json_escape(std::get<0>(kv))
            << "\", \"molecules\": " << std::get<1>(kv)
            << ", \"atoms_per_molecule\": " << std::get<2>(kv) << "}";
        out << (i + 1 < molecules_by_type_ordered.size() ? ",\n" : "\n");
    }
    out << "  ]\n";
    out << "}\n";
}

int main(int argc, char** argv) {
    // usage: ./mini_gro_topology system.gro "SOL,NA,CL" 2 1.25 [output.json]
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " file.gro \"SOL,NA,CL\" axis middle [output.json]\n";
        return 1;
    }

    const std::string gro_path = argv[1];
    const auto solutions = parse_names(argv[2]);
    const int axis = std::stoi(argv[3]);       // 0=x,1=y,2=z
    const double middle = std::stod(argv[4]);  // nm
    const std::string json_out_path = (argc >= 6) ? argv[5] : "";
    if (axis < 0 || axis > 2) throw std::runtime_error("axis must be 0..2");

    std::ifstream in(gro_path);
    if (!in) throw std::runtime_error("cannot open gro file");

    std::string line;
    std::getline(in, line); // title
    if (!std::getline(in, line)) throw std::runtime_error("missing natoms line");
    const int natoms = std::stoi(trim(line));

    const double INF = std::numeric_limits<double>::infinity();
    std::array<double,6> box_local{INF, INF, INF, -INF, -INF, -INF};
    std::array<double,6> channel  {INF, INF, INF, -INF, -INF, -INF};
    std::array<double,3> gro_box  {0.0, 0.0, 0.0};
    channel[axis] = -INF;      // lower wall = max below middle
    channel[axis + 3] = +INF;  // upper wall = min above middle

    int natoms_filtered = 0;
    std::vector<std::tuple<std::string, int, int>> molecules_by_type_ordered;
    std::unordered_map<std::string, size_t> type_to_index;

    // residue block tracking (for molecule counting)
    bool have_res = false;
    int cur_resnr = -1;
    std::string cur_resname;
    int cur_res_atoms = 0;
    auto flush_res = [&]() {
        if (!have_res) return;
        if (cur_res_atoms <= 0) return;
        auto it = type_to_index.find(cur_resname);
        if (it == type_to_index.end()) {
            type_to_index[cur_resname] = molecules_by_type_ordered.size();
            molecules_by_type_ordered.push_back({cur_resname, 1, cur_res_atoms});
        } else {
            auto& entry = molecules_by_type_ordered[it->second];
            if (std::get<2>(entry) != cur_res_atoms) {
                throw std::runtime_error(
                    "inconsistent atoms_per_molecule for type " + cur_resname +
                    ": saw " + std::to_string(std::get<2>(entry)) +
                    " and " + std::to_string(cur_res_atoms));
            }
            std::get<1>(entry)++;
        }
    };

    for (int i = 0; i < natoms; ++i) {
        if (!std::getline(in, line) || line.size() < 44) throw std::runtime_error("bad atom line");

        const int resnr = std::stoi(line.substr(0, 5));
        const std::string resname = upper(trim(line.substr(5, 5)));
        const double x = std::stod(line.substr(20, 8));
        const double y = std::stod(line.substr(28, 8));
        const double z = std::stod(line.substr(36, 8));
        const double p[3] = {x, y, z};

        if (!have_res) {
            have_res = true; cur_resnr = resnr; cur_resname = resname; cur_res_atoms = 1;
        } else if (resnr != cur_resnr || resname != cur_resname) {
            flush_res();
            cur_resnr = resnr; cur_resname = resname; cur_res_atoms = 1;
        } else {
            cur_res_atoms++;
        }

        const bool is_solution = solutions.count(resname) > 0;
        if (is_solution) {
            natoms_filtered++;
            for (int k = 0; k < 3; ++k) {
                if (p[k] < box_local[k]) box_local[k] = p[k];
                if (p[k] > box_local[k + 3]) box_local[k + 3] = p[k];
            }
        } else {
            for (int k = 0; k < 3; ++k) {
                if (k == axis) {
                    if (p[k] < middle && p[k] > channel[axis]) channel[axis] = p[k];
                    if (p[k] > middle && p[k] < channel[axis + 3]) channel[axis + 3] = p[k];
                } else {
                    if (p[k] < channel[k]) channel[k] = p[k];
                    if (p[k] > channel[k + 3]) channel[k + 3] = p[k];
                }
            }
        }
    }
    flush_res();

    if (!std::getline(in, line)) throw std::runtime_error("missing box line");
    std::istringstream boxs(line);
    boxs >> gro_box[0] >> gro_box[1] >> gro_box[2];

    std::cout << "channel min: " << channel[0] << " " << channel[1] << " " << channel[2] << "\n";
    std::cout << "channel max: " << channel[3] << " " << channel[4] << " " << channel[5] << "\n";
    std::cout << "box_local min/max: "
              << box_local[0] << " " << box_local[1] << " " << box_local[2] << " / "
              << box_local[3] << " " << box_local[4] << " " << box_local[5] << "\n";
    std::cout << "gro_box: " << gro_box[0] << " " << gro_box[1] << " " << gro_box[2] << "\n";
    std::cout << "natoms_filtered: " << natoms_filtered << "\n";
    for (const auto& kv : molecules_by_type_ordered) {
        const std::string& name = std::get<0>(kv);
        const int nmol = std::get<1>(kv);
        const int apm = std::get<2>(kv);
        const bool selected = solutions.count(name) > 0;
        std::cout << "type " << name << " molecules=" << nmol
                  << " atoms_per_molecule=" << apm
                  << " selected=" << (selected ? 1 : 0) << "\n";
    }

    if (!json_out_path.empty()) {
        write_json_output(
            json_out_path,
            channel,
            box_local,
            gro_box,
            natoms_filtered,
            molecules_by_type_ordered);
        std::cout << "json written: " << json_out_path << "\n";
    }

    return 0;
}
