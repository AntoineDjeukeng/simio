#include "middle_reservoir_monitor.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace simio::middle_reservoir {
namespace {

constexpr double kCrossEps = 1e-12;

enum class CrossingDirection : uint8_t { None, LeftToRight, RightToLeft };

struct Crossing {
    CrossingDirection direction = CrossingDirection::None;
    double fraction = 0.0;
};

Crossing detect_crossing(double previous_x, double dx_minimum_image, double plane_x, double box_x) {
    const double image = plane_x + std::round((previous_x - plane_x) / box_x) * box_x;
    const double d0 = previous_x - image;
    const double d1 = previous_x + dx_minimum_image - image;
    if (!((d0 < 0.0 && d1 >= 0.0) || (d0 >= 0.0 && d1 < 0.0))) return {};

    const double denominator = d1 - d0;
    if (std::abs(denominator) <= kCrossEps) return {};
    double fraction = -d0 / denominator;
    if (fraction < -kCrossEps || fraction > 1.0 + kCrossEps) return {};
    fraction = std::max(0.0, std::min(1.0, fraction));
    return Crossing{d0 < 0.0 ? CrossingDirection::LeftToRight
                             : CrossingDirection::RightToLeft,
                    fraction};
}

int species_index(simio::MolType type) {
    if (type == simio::MolType::Water) return 0;
    if (type == simio::MolType::Cation) return 1;
    if (type == simio::MolType::Anion) return 2;
    return -1;
}

bool in_z_aperture(double z, const GateGeometry& gate, const simio::Pbc3D& pbc) {
    const double zw = pbc.wrap_pos(2, z);
    const double zmin = pbc.wrap_pos(2, gate.zmin_nm);
    const double zmax = pbc.wrap_pos(2, gate.zmax_nm);
    if (zmin <= zmax) return zw >= zmin && zw < zmax;
    return zw >= zmin || zw < zmax;
}

std::string json_escape(const std::string& value) {
    std::string out;
    for (const char c : value) {
        if (c == '\\' || c == '"') out.push_back('\\');
        out.push_back(c);
    }
    return out;
}

}  // namespace

MiddleReservoirMonitor::MiddleReservoirMonitor(const ReactorSetup& setup) : setup_(setup) {
    previous_key_.resize(setup_.topology.mols.size());
    has_previous_.assign(setup_.topology.mols.size(), 0);
}

simio::Vec3d MiddleReservoirMonitor::molecule_key(const simio::Frame& frame,
                                                   const simio::MolSpan& molecule) const {
    if (molecule.type == simio::MolType::Water) {
        if (molecule.natoms != 3) {
            throw std::runtime_error("MiddleReservoirMonitor requires three-site water");
        }
        return simio::water_com_pbc_unwrapped(frame.pbc,
                                              frame.atoms.pos(molecule.first),
                                              frame.atoms.pos(molecule.first + 1),
                                              frame.atoms.pos(molecule.first + 2));
    }
    return frame.atoms.pos(molecule.first);
}

void MiddleReservoirMonitor::process_frame(const simio::Frame& frame, int frame_index) {
    if (frame.pbc.L[0] <= 0.0 || frame.pbc.L[2] <= 0.0) {
        throw std::runtime_error("MiddleReservoirMonitor received invalid box lengths");
    }

    FrameRow row;
    row.frame_index = frame_index;
    row.step = frame.step;
    row.time_ps = frame.time_ps;

    for (size_t molecule_id = 0; molecule_id < setup_.topology.mols.size(); ++molecule_id) {
        const simio::MolSpan& molecule = setup_.topology.mols[molecule_id];
        const int species = species_index(molecule.type);
        if (species < 0) continue;

        simio::Vec3d current = molecule_key(frame, molecule);
        frame.pbc.wrap_pos3(current);
        if (current.v[0] >= setup_.left_gate.x_nm && current.v[0] < setup_.right_gate.x_nm) {
            ++row.count[static_cast<size_t>(species)];
        }

        if (has_previous_[molecule_id]) {
            const simio::Vec3d delta = frame.pbc.min_image(current - previous_key_[molecule_id]);

            const Crossing left = detect_crossing(previous_key_[molecule_id].v[0],
                                                  delta.v[0],
                                                  setup_.left_gate.x_nm,
                                                  frame.pbc.L[0]);
            if (left.direction != CrossingDirection::None) {
                const double crossing_z = previous_key_[molecule_id].v[2] + left.fraction * delta.v[2];
                if (in_z_aperture(crossing_z, setup_.left_gate, frame.pbc)) {
                    if (left.direction == CrossingDirection::LeftToRight) {
                        ++row.left.entered[static_cast<size_t>(species)];
                        ++left_cumulative_net_[static_cast<size_t>(species)];
                    } else {
                        ++row.left.exited[static_cast<size_t>(species)];
                        --left_cumulative_net_[static_cast<size_t>(species)];
                    }
                }
            }

            const Crossing right = detect_crossing(previous_key_[molecule_id].v[0],
                                                   delta.v[0],
                                                   setup_.right_gate.x_nm,
                                                   frame.pbc.L[0]);
            if (right.direction != CrossingDirection::None) {
                const double crossing_z = previous_key_[molecule_id].v[2] + right.fraction * delta.v[2];
                if (in_z_aperture(crossing_z, setup_.right_gate, frame.pbc)) {
                    if (right.direction == CrossingDirection::RightToLeft) {
                        ++row.right.entered[static_cast<size_t>(species)];
                        ++right_cumulative_net_[static_cast<size_t>(species)];
                    } else {
                        ++row.right.exited[static_cast<size_t>(species)];
                        --right_cumulative_net_[static_cast<size_t>(species)];
                    }
                }
            }
        }

        previous_key_[molecule_id] = current;
        has_previous_[molecule_id] = 1;
    }

    row.left_cumulative_net = left_cumulative_net_;
    row.right_cumulative_net = right_cumulative_net_;
    rows_.push_back(row);
}

void MiddleReservoirMonitor::write_csv(const std::string& path) const {
    if (rows_.empty()) throw std::runtime_error("No middle-reservoir frames were processed");
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Failed to open output CSV: " + path);

    output << "frame_idx,step,time_ps,water_count,na_count,cl_count";
    constexpr const char* species_names[] = {"water", "na", "cl"};
    for (const char* side : {"left", "right"}) {
        for (const char* species : species_names) {
            output << ',' << side << '_' << species << "_in"
                   << ',' << side << '_' << species << "_out"
                   << ',' << side << '_' << species << "_cum_net_in";
        }
    }
    output << '\n';

    output << std::setprecision(12);
    for (const FrameRow& row : rows_) {
        output << row.frame_index << ',' << row.step << ',' << row.time_ps << ',' << row.count[0]
               << ',' << row.count[1] << ',' << row.count[2];
        for (size_t species = 0; species < 3; ++species) {
            output << ',' << row.left.entered[species] << ',' << row.left.exited[species] << ','
                   << row.left_cumulative_net[species];
        }
        for (size_t species = 0; species < 3; ++species) {
            output << ',' << row.right.entered[species] << ',' << row.right.exited[species] << ','
                   << row.right_cumulative_net[species];
        }
        output << '\n';
    }
}

void write_setup_json(const ReactorSetup& setup,
                      const std::string& gro_path,
                      const std::string& report_path,
                      const std::string& output_path) {
    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("Failed to open setup output: " + output_path);
    const IonInsertionReportSummary& report = setup.report;
    const double middle_length = setup.middle_length_nm();

    output << std::setprecision(12);
    output << "{\n"
           << "  \"source_gro\": \"" << json_escape(gro_path) << "\",\n"
           << "  \"ion_insertion_report\": \"" << json_escape(report_path) << "\",\n"
           << "  \"composition\": {\"water\": " << setup.nwater << ", \"na\": " << setup.nna
           << ", \"cl\": " << setup.ncl << ", \"trajectory_natoms\": " << setup.natoms
           << ", \"gro_natoms\": " << setup.gro_natoms << "},\n"
           << "  \"box_nm\": [" << setup.gro_box_nm[0] << ", " << setup.gro_box_nm[1] << ", "
           << setup.gro_box_nm[2] << "],\n"
           << "  \"report_summary\": {\n"
           << "    \"cleanup\": {\"waters_before\": " << report.waters_before_cleanup
           << ", \"waters_removed\": " << report.waters_removed_as_overlap
           << ", \"waters_after\": " << report.waters_after_cleanup
           << ", \"overlap_cutoff_nm\": " << report.overlap_cutoff_nm
           << ", \"minimum_surface_distance_nm\": "
           << report.minimum_water_surface_distance_nm << "},\n"
           << "    \"surface_charge\": {\"summed_charge_e\": " << report.summed_surface_charge_e
           << ", \"integer_charge_e\": " << report.integer_surface_charge_e << "},\n"
           << "    \"ion_counts\": {\"target_molality\": " << report.target_molality
           << ", \"salt_pairs\": " << report.salt_pairs
           << ", \"extra_na_for_surface\": " << report.extra_na_for_surface
           << ", \"extra_cl_for_surface\": " << report.extra_cl_for_surface
           << ", \"total_na\": " << report.total_na << ", \"total_cl\": " << report.total_cl
           << "},\n"
           << "    \"reservoirs\": {\n"
           << "      \"left\": {\"total_waters\": " << report.left.total_waters
           << ", \"eligible_waters\": " << report.left.eligible_waters
           << ", \"allocated_na\": " << report.left.allocated_na
           << ", \"allocated_cl\": " << report.left.allocated_cl
           << ", \"expected_water_after_insertion\": "
           << report.left.expected_water_after_insertion() << "},\n"
           << "      \"middle\": {\"total_waters\": " << report.middle.total_waters
           << ", \"eligible_waters\": " << report.middle.eligible_waters
           << ", \"allocated_na\": " << report.middle.allocated_na
           << ", \"allocated_cl\": " << report.middle.allocated_cl
           << ", \"expected_water_after_insertion\": "
           << report.middle.expected_water_after_insertion() << "},\n"
           << "      \"right\": {\"total_waters\": " << report.right.total_waters
           << ", \"eligible_waters\": " << report.right.eligible_waters
           << ", \"allocated_na\": " << report.right.allocated_na
           << ", \"allocated_cl\": " << report.right.allocated_cl
           << ", \"expected_water_after_insertion\": "
           << report.right.expected_water_after_insertion() << "}\n"
           << "    },\n"
           << "    \"expected_final\": {\"water\": " << report.expected_water
           << ", \"na\": " << report.expected_na << ", \"cl\": " << report.expected_cl
           << ", \"total_charge_e\": " << report.expected_total_charge_e << "}\n"
           << "  },\n"
           << "  \"middle_reservoir\": {\n"
           << "    \"origin_x_absolute_nm\": " << setup.left_gate.x_nm << ",\n"
           << "    \"length_nm\": " << middle_length << ",\n"
           << "    \"absolute_interval_nm\": [" << setup.left_gate.x_nm << ", "
           << setup.right_gate.x_nm << "],\n"
           << "    \"relative_interval_nm\": [0, " << middle_length << "],\n"
           << "    \"left_gate\": {\"x_absolute_nm\": " << setup.left_gate.x_nm
           << ", \"x_relative_nm\": 0, \"zmin_nm\": " << setup.left_gate.zmin_nm
           << ", \"zmax_nm\": " << setup.left_gate.zmax_nm << ", \"wall_type\": \""
           << setup.left_gate.wall_type << "\"},\n"
           << "    \"right_gate\": {\"x_absolute_nm\": " << setup.right_gate.x_nm
           << ", \"x_relative_nm\": " << middle_length << ", \"zmin_nm\": "
           << setup.right_gate.zmin_nm << ", \"zmax_nm\": " << setup.right_gate.zmax_nm
           << ", \"wall_type\": \"" << setup.right_gate.wall_type << "\"}\n"
           << "  }\n"
           << "}\n";
}

}  // namespace simio::middle_reservoir
