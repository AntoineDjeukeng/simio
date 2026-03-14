#pragma once

#include <string>

#include "coord_x.hpp"
#include "density_x.hpp"
#include "dipole_x.hpp"
#include "water_atom_density_x.hpp"

namespace simio::analysis {

struct AxialProfileXExportConfig {
    double xmin = 0.0;
    double xmax = 0.0;
    double box_length_x_nm = 0.0;
};

void write_axial_profile_x_csv(const std::string& path,
                               const AxialProfileXExportConfig& cfg,
                               const DensityXAnalyzer& density,
                               const WaterAtomDensityXAnalyzer& water_atom_density,
                               const DipoleXAnalyzer& dipole,
                               const CoordXAnalyzer& coord);

}  // namespace simio::analysis
