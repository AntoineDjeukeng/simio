#include <cassert>
#include <cmath>
#include <numeric>

#include "simio/layered_pipeline.hpp"
#include "simio/properties/DensityXZ.hpp"

namespace {
constexpr double kTol = 1e-12;
}

int main() {
    simio::Topology topo;
    topo.mols.push_back(simio::MolSpan{0, 1, simio::MolType::Cation});
    topo.build_type_lists();

    simio::Frame fr;
    fr.pbc = simio::Pbc3D(10.0, 5.0, 3.0);
    fr.atoms.x.assign(1, 7.5);
    fr.atoms.y.assign(1, 2.0);
    fr.atoms.z.assign(1, 1.2);

    simio::layered::Layer0Config l0;
    l0.nx = 10;
    l0.nz = 6;
    l0.channel_length_x = 10.0;
    l0.channel_height_z = 3.0;
    simio::layered::Pipeline pipe(l0, {});

    simio::properties::DensityXZConfig dcfg;
    dcfg.roi.x_min = 7.0;
    dcfg.roi.x_max = 9.0;
    dcfg.roi.z_min = 1.0;
    dcfg.roi.z_max = 2.0;
    dcfg.roi.x_mode = simio::properties::DensityXMode::ChannelWindow;
    dcfg.nx = 4;
    dcfg.nz = 2;
    dcfg.normalize_number_density = true;

    auto density = std::make_unique<simio::properties::DensityXZProperty>(dcfg);
    simio::properties::DensityXZProperty* density_ptr = density.get();
    pipe.add_property(std::move(density));

    fr.step = 0;
    fr.time_ps = 0.0;
    pipe.process_frame(topo, fr);

    const auto& f0 = density_ptr->frames().back();
    assert(f0.selected_total == 1);
    assert(f0.binned_total == 1);
    assert(f0.bin_oob_count == 0);
    assert(f0.selected_cation == 1);
    assert(f0.selected_water == 0);
    assert(f0.selected_anion == 0);
    assert(f0.selected_other == 0);
    assert(std::accumulate(f0.counts.begin(), f0.counts.end(), int64_t{0}) == 1);
    assert(std::fabs(f0.dx - 0.5) < kTol);
    assert(std::fabs(f0.dz - 0.5) < kTol);
    assert(std::fabs(f0.bin_volume - 1.25) < kTol);
    assert(std::fabs(f0.rho[1] - 0.8) < kTol);

    fr.step = 1;
    fr.time_ps = 1.0;
    fr.atoms.x[0] = 9.5;  // outside [7,9)
    fr.atoms.z[0] = 1.2;
    pipe.process_frame(topo, fr);

    const auto& f1 = density_ptr->frames().back();
    assert(f1.selected_total == 0);
    assert(f1.binned_total == 0);
    assert(std::accumulate(f1.counts.begin(), f1.counts.end(), int64_t{0}) == 0);

    return 0;
}

