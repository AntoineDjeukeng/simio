#include "reactor_setup.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace simio::middle_reservoir {
namespace {

struct GroAtom {
    int residue_id = 0;
    std::string residue_name;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct WallComponent {
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    std::vector<const GroAtom*> atoms;
};

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return value;
}

std::string read_text(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Failed to open file: " + path);
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

std::string json_object(const std::string& text, const std::string& key) {
    const std::string marker = "\"" + key + "\"";
    const size_t key_pos = text.find(marker);
    if (key_pos == std::string::npos) {
        throw std::runtime_error("Insertion report is missing object: " + key);
    }
    const size_t begin = text.find('{', key_pos + marker.size());
    if (begin == std::string::npos) {
        throw std::runtime_error("Invalid JSON object for: " + key);
    }
    int depth = 0;
    bool in_string = false;
    bool escaped = false;
    for (size_t i = begin; i < text.size(); ++i) {
        const char c = text[i];
        if (in_string) {
            if (escaped) {
                escaped = false;
            } else if (c == '\\') {
                escaped = true;
            } else if (c == '"') {
                in_string = false;
            }
            continue;
        }
        if (c == '"') {
            in_string = true;
        } else if (c == '{') {
            ++depth;
        } else if (c == '}') {
            --depth;
            if (depth == 0) return text.substr(begin, i - begin + 1);
        }
    }
    throw std::runtime_error("Unterminated JSON object for: " + key);
}

int json_int(const std::string& object, const std::string& key) {
    const std::regex expression("\"" + key + "\"\\s*:\\s*(-?[0-9]+)");
    std::smatch match;
    if (!std::regex_search(object, match, expression)) {
        throw std::runtime_error("Insertion report is missing integer: " + key);
    }
    return std::stoi(match[1].str());
}

double json_double(const std::string& object, const std::string& key) {
    const std::regex expression(
        "\\\"" + key + "\\\"\\s*:\\s*([-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)");
    std::smatch match;
    if (!std::regex_search(object, match, expression)) {
        throw std::runtime_error("Insertion report is missing number: " + key);
    }
    return std::stod(match[1].str());
}

ReservoirReport parse_reservoir_report(const std::string& reservoirs,
                                       const std::string& na_allocation,
                                       const std::string& cl_allocation,
                                       const std::string& name) {
    const std::string reservoir = json_object(reservoirs, name);
    ReservoirReport out;
    out.total_waters = json_int(reservoir, "total_waters");
    out.eligible_waters = json_int(reservoir, "eligible_waters");
    out.allocated_na = json_int(na_allocation, name);
    out.allocated_cl = json_int(cl_allocation, name);
    return out;
}

IonInsertionReportSummary parse_report_summary(const std::string& report) {
    IonInsertionReportSummary out;

    const std::string cleanup = json_object(report, "cleanup");
    out.waters_before_cleanup = json_int(cleanup, "waters_before_cleanup");
    out.waters_removed_as_overlap = json_int(cleanup, "waters_removed_as_overlap");
    out.waters_after_cleanup = json_int(cleanup, "waters_after_cleanup");
    out.overlap_cutoff_nm = json_double(cleanup, "overlap_cutoff_nm");
    out.minimum_water_surface_distance_nm =
        json_double(cleanup, "minimum_water_surface_distance_nm");

    const std::string surface_charge = json_object(report, "surface_charge");
    out.summed_surface_charge_e = json_double(surface_charge, "summed_charge_e");
    out.integer_surface_charge_e = json_int(surface_charge, "integer_charge_e");

    const std::string ions = json_object(report, "ion_counts");
    out.target_molality = json_double(ions, "target_molality");
    out.salt_pairs = json_int(ions, "salt_pairs");
    out.extra_na_for_surface = json_int(ions, "extra_na_for_surface");
    out.extra_cl_for_surface = json_int(ions, "extra_cl_for_surface");
    out.total_na = json_int(ions, "total_na");
    out.total_cl = json_int(ions, "total_cl");

    const std::string reservoirs = json_object(report, "reservoirs");
    const std::string allocation = json_object(report, "allocation");
    const std::string na_allocation = json_object(allocation, "na");
    const std::string cl_allocation = json_object(allocation, "cl");
    out.left = parse_reservoir_report(reservoirs, na_allocation, cl_allocation, "left");
    out.middle = parse_reservoir_report(reservoirs, na_allocation, cl_allocation, "middle");
    out.right = parse_reservoir_report(reservoirs, na_allocation, cl_allocation, "right");

    const std::string expected = json_object(report, "expected_final");
    out.expected_water = json_int(expected, "water_molecules");
    out.expected_na = json_int(expected, "na");
    out.expected_cl = json_int(expected, "cl");
    out.expected_total_charge_e = json_int(expected, "total_charge_e");

    if (out.waters_before_cleanup - out.waters_removed_as_overlap != out.waters_after_cleanup) {
        throw std::runtime_error("Insertion report cleanup water balance is inconsistent");
    }
    if (out.waters_after_cleanup - out.total_na - out.total_cl != out.expected_water) {
        throw std::runtime_error("Insertion report final water balance is inconsistent");
    }
    if (out.total_na != out.expected_na || out.total_cl != out.expected_cl) {
        throw std::runtime_error("Insertion report ion_counts and expected_final disagree");
    }
    const int allocated_na = out.left.allocated_na + out.middle.allocated_na + out.right.allocated_na;
    const int allocated_cl = out.left.allocated_cl + out.middle.allocated_cl + out.right.allocated_cl;
    if (allocated_na != out.total_na || allocated_cl != out.total_cl) {
        throw std::runtime_error("Insertion report reservoir allocations do not sum to ion totals");
    }
    return out;
}

std::map<std::string, int> report_wall_counts(const std::string& report) {
    const std::string details = json_object(report, "details");
    std::map<std::string, int> counts;
    const std::regex entry(
        "\"([^\"]+)\"\\s*:\\s*\\{[^\\{\\}]*\"atom_count\"\\s*:\\s*([0-9]+)",
        std::regex::icase);
    for (std::sregex_iterator it(details.begin(), details.end(), entry), end; it != end; ++it) {
        counts[upper((*it)[1].str())] = std::stoi((*it)[2].str());
    }
    if (counts.empty()) {
        throw std::runtime_error("Insertion report surface_charge.details has no atom counts");
    }
    return counts;
}

simio::MolType classify_residue(const std::string& residue_name) {
    const std::string name = upper(residue_name);
    if (name == "SOL" || name == "WAT" || name == "HOH") return simio::MolType::Water;
    if (name == "NA" || name == "NA+" || name == "SOD") return simio::MolType::Cation;
    if (name == "CL" || name == "CL-" || name == "CLA") return simio::MolType::Anion;
    return simio::MolType::Other;
}

std::vector<GroAtom> read_gro(const std::string& path, std::array<double, 3>& box) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Failed to open GRO: " + path);

    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("GRO is missing title: " + path);
    if (!std::getline(input, line)) throw std::runtime_error("GRO is missing atom count: " + path);
    const int natoms = std::stoi(trim(line));
    if (natoms <= 0) throw std::runtime_error("GRO atom count must be positive: " + path);

    std::vector<GroAtom> atoms;
    atoms.reserve(static_cast<size_t>(natoms));
    for (int i = 0; i < natoms; ++i) {
        if (!std::getline(input, line) || line.size() < 44) {
            throw std::runtime_error("Invalid GRO atom line " + std::to_string(i + 1));
        }
        GroAtom atom;
        atom.residue_id = std::stoi(line.substr(0, 5));
        atom.residue_name = trim(line.substr(5, 5));
        atom.x = std::stod(line.substr(20, 8));
        atom.y = std::stod(line.substr(28, 8));
        atom.z = std::stod(line.substr(36, 8));
        atoms.push_back(std::move(atom));
    }

    if (!std::getline(input, line)) throw std::runtime_error("GRO is missing box line: " + path);
    std::istringstream box_line(line);
    if (!(box_line >> box[0] >> box[1] >> box[2])) {
        throw std::runtime_error("Only orthorhombic GRO box lines are supported: " + path);
    }
    return atoms;
}

simio::Topology build_topology(const std::vector<GroAtom>& atoms) {
    simio::Topology topology;
    size_t first = 0;
    while (first < atoms.size()) {
        size_t end = first + 1;
        while (end < atoms.size() && atoms[end].residue_id == atoms[first].residue_id &&
               atoms[end].residue_name == atoms[first].residue_name) {
            ++end;
        }
        topology.mols.push_back(simio::MolSpan{static_cast<int>(first),
                                               static_cast<int>(end - first),
                                               classify_residue(atoms[first].residue_name)});
        first = end;
    }
    topology.build_type_lists();
    return topology;
}

simio::Topology mobile_only_topology(const simio::Topology& gro_topology) {
    simio::Topology topology;
    int first = 0;
    for (const simio::MolSpan& molecule : gro_topology.mols) {
        if (molecule.type == simio::MolType::Other) continue;
        topology.mols.push_back(simio::MolSpan{first, molecule.natoms, molecule.type});
        first += molecule.natoms;
    }
    topology.build_type_lists();
    return topology;
}

std::vector<WallComponent> find_wall_components(const std::vector<GroAtom>& atoms) {
    std::vector<double> x_values;
    for (const GroAtom& atom : atoms) {
        if (classify_residue(atom.residue_name) == simio::MolType::Other) x_values.push_back(atom.x);
    }
    if (x_values.empty()) throw std::runtime_error("GRO contains no wall atoms");
    std::sort(x_values.begin(), x_values.end());
    x_values.erase(std::unique(x_values.begin(), x_values.end(), [](double a, double b) {
                       return std::abs(a - b) < 1e-6;
                   }),
                   x_values.end());

    std::vector<double> spacings;
    for (size_t i = 1; i < x_values.size(); ++i) {
        const double gap = x_values[i] - x_values[i - 1];
        if (gap > 1e-6) spacings.push_back(gap);
    }
    if (spacings.empty()) throw std::runtime_error("Wall geometry has no x extent");
    std::sort(spacings.begin(), spacings.end());
    const double lattice_spacing = spacings[spacings.size() / 2];
    const double split_gap = std::max(0.5, 4.0 * lattice_spacing);

    std::vector<std::pair<double, double>> ranges;
    double begin = x_values.front();
    double previous = begin;
    for (size_t i = 1; i < x_values.size(); ++i) {
        if (x_values[i] - previous > split_gap) {
            ranges.emplace_back(begin, previous);
            begin = x_values[i];
        }
        previous = x_values[i];
    }
    ranges.emplace_back(begin, previous);
    if (ranges.size() != 2) {
        throw std::runtime_error("Expected exactly two separated wall components, found " +
                                 std::to_string(ranges.size()));
    }

    std::vector<WallComponent> components(2);
    for (size_t i = 0; i < ranges.size(); ++i) {
        components[i].xmin = ranges[i].first;
        components[i].xmax = ranges[i].second;
    }
    for (const GroAtom& atom : atoms) {
        if (classify_residue(atom.residue_name) != simio::MolType::Other) continue;
        for (WallComponent& component : components) {
            if (atom.x >= component.xmin - 1e-6 && atom.x <= component.xmax + 1e-6) {
                component.atoms.push_back(&atom);
                break;
            }
        }
    }
    return components;
}

GateGeometry gate_from_component(const WallComponent& component, double gate_x) {
    std::map<std::string, std::vector<double>> z_by_type;
    for (const GroAtom* atom : component.atoms) {
        z_by_type[upper(atom->residue_name)].push_back(atom->z);
    }
    if (z_by_type.size() < 2) {
        throw std::runtime_error("Wall component has no distinct channel-surface atom type");
    }

    // CAG is the neutral scaffold in generated reactors. Of the remaining
    // types, the most populated one defines this channel's two facing surfaces.
    auto selected = z_by_type.end();
    for (auto it = z_by_type.begin(); it != z_by_type.end(); ++it) {
        if (it->first == "CAG") continue;
        if (selected == z_by_type.end() || it->second.size() > selected->second.size()) selected = it;
    }
    if (selected == z_by_type.end()) {
        throw std::runtime_error("Wall component only contains CAG; cannot infer channel aperture");
    }

    const auto z_bounds = std::minmax_element(selected->second.begin(), selected->second.end());
    if (*z_bounds.second <= *z_bounds.first) {
        throw std::runtime_error("Channel-surface type " + selected->first +
                                 " does not define two z surfaces");
    }
    return GateGeometry{gate_x, *z_bounds.first, *z_bounds.second, selected->first};
}

int molecule_count(const simio::Topology& topology, simio::MolType type) {
    return static_cast<int>(topology.mol_ids_by_type[static_cast<int>(type)].size());
}

}  // namespace

ReactorSetup load_reactor_setup(const std::string& gro_path,
                                const std::string& ion_insertion_report_path) {
    const std::string report = read_text(ion_insertion_report_path);
    ReactorSetup setup;
    setup.report = parse_report_summary(report);
    setup.nwater = setup.report.expected_water;
    setup.nna = setup.report.expected_na;
    setup.ncl = setup.report.expected_cl;
    const std::map<std::string, int> expected_wall_counts = report_wall_counts(report);

    std::vector<GroAtom> atoms = read_gro(gro_path, setup.gro_box_nm);
    std::map<std::string, int> gro_wall_counts;
    for (const GroAtom& atom : atoms) {
        if (classify_residue(atom.residue_name) == simio::MolType::Other) {
            ++gro_wall_counts[upper(atom.residue_name)];
        }
    }
    for (const auto& [name, expected_count] : expected_wall_counts) {
        const int actual_count = gro_wall_counts[name];
        if (actual_count != expected_count) {
            throw std::runtime_error("GRO/report wall count mismatch for " + name + ": GRO=" +
                                     std::to_string(actual_count) + " report=" +
                                     std::to_string(expected_count));
        }
    }
    setup.gro_natoms = static_cast<int>(atoms.size());
    setup.topology = mobile_only_topology(build_topology(atoms));
    for (const simio::MolSpan& molecule : setup.topology.mols) setup.natoms += molecule.natoms;

    const int gro_water = molecule_count(setup.topology, simio::MolType::Water);
    const int gro_na = molecule_count(setup.topology, simio::MolType::Cation);
    const int gro_cl = molecule_count(setup.topology, simio::MolType::Anion);
    if (gro_water != setup.nwater || gro_na != setup.nna || gro_cl != setup.ncl) {
        throw std::runtime_error("GRO/report composition mismatch: GRO water/Na/Cl=" +
                                 std::to_string(gro_water) + "/" + std::to_string(gro_na) + "/" +
                                 std::to_string(gro_cl) + ", report=" +
                                 std::to_string(setup.nwater) + "/" + std::to_string(setup.nna) +
                                 "/" + std::to_string(setup.ncl));
    }

    const std::vector<WallComponent> walls = find_wall_components(atoms);
    setup.left_gate = gate_from_component(walls[0], walls[0].xmax);
    setup.right_gate = gate_from_component(walls[1], walls[1].xmin);
    if (setup.right_gate.x_nm <= setup.left_gate.x_nm) {
        throw std::runtime_error("Derived middle-reservoir interval is empty");
    }
    return setup;
}

}  // namespace simio::middle_reservoir
