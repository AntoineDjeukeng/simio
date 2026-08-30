#pragma once

#include "simio/simio.hpp"

#include <array>
#include <string>

namespace simio::middle_reservoir {

struct GateGeometry {
    double x_nm = 0.0;
    double zmin_nm = 0.0;
    double zmax_nm = 0.0;
    std::string wall_type;
};

struct ReservoirReport {
    int total_waters = 0;
    int eligible_waters = 0;
    int allocated_na = 0;
    int allocated_cl = 0;

    int expected_water_after_insertion() const {
        return total_waters - allocated_na - allocated_cl;
    }
};

struct IonInsertionReportSummary {
    int waters_before_cleanup = 0;
    int waters_removed_as_overlap = 0;
    int waters_after_cleanup = 0;
    double overlap_cutoff_nm = 0.0;
    double minimum_water_surface_distance_nm = 0.0;

    double summed_surface_charge_e = 0.0;
    int integer_surface_charge_e = 0;

    double target_molality = 0.0;
    int salt_pairs = 0;
    int extra_na_for_surface = 0;
    int extra_cl_for_surface = 0;
    int total_na = 0;
    int total_cl = 0;

    ReservoirReport left;
    ReservoirReport middle;
    ReservoirReport right;

    int expected_water = 0;
    int expected_na = 0;
    int expected_cl = 0;
    int expected_total_charge_e = 0;
};

struct ReactorSetup {
    simio::Topology topology;
    int natoms = 0;
    int gro_natoms = 0;
    int nwater = 0;
    int nna = 0;
    int ncl = 0;
    std::array<double, 3> gro_box_nm{0.0, 0.0, 0.0};
    GateGeometry left_gate;
    GateGeometry right_gate;
    IonInsertionReportSummary report;

    double middle_length_nm() const { return right_gate.x_nm - left_gate.x_nm; }
};

// The insertion report is authoritative for composition and preparation counts.
// The final GRO is authoritative for atom ordering and wall geometry.
ReactorSetup load_reactor_setup(const std::string& gro_path,
                                const std::string& ion_insertion_report_path);

}  // namespace simio::middle_reservoir
