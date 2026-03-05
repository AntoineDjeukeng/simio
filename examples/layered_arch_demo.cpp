#include <iostream>
#include <memory>
#include <vector>

#include "simio/layered_pipeline.hpp"
#include "simio/properties/DensityXZ.hpp"
#include "simio/properties/DipoleXZ.hpp"
#include "simio/properties/GatingCenterPlane.hpp"
#include "simio/properties/GatingPbcWrap.hpp"

namespace {

simio::Topology make_topology() {
    simio::Topology topo;
    topo.mols.push_back(simio::MolSpan{0, 3, simio::MolType::Water});
    topo.mols.push_back(simio::MolSpan{3, 3, simio::MolType::Water});
    topo.mols.push_back(simio::MolSpan{6, 1, simio::MolType::Cation});
    topo.mols.push_back(simio::MolSpan{7, 1, simio::MolType::Anion});
    topo.build_type_lists();
    return topo;
}

void fill_frame(simio::Frame& fr, int frame_id) {
    fr.step = frame_id;
    fr.time_ps = frame_id * 1.0;

    const double w0x = (frame_id == 0) ? 2.45 : ((frame_id == 1) ? 2.60 : 2.40);
    const double clx = (frame_id == 0) ? 4.90 : ((frame_id == 1) ? 0.10 : 4.90);

    fr.atoms.x = {w0x, w0x + 0.10, w0x - 0.08,
                  3.00 + 0.03 * frame_id, 3.10 + 0.03 * frame_id, 2.95 + 0.03 * frame_id,
                  2.00 + 0.01 * frame_id, clx};

    fr.atoms.y = {1.0, 1.02, 0.95,
                  2.0, 2.02, 1.95,
                  4.0, 4.2};

    fr.atoms.z = {5.0 + 0.01 * frame_id, 5.05 + 0.01 * frame_id, 4.95 + 0.01 * frame_id,
                  5.8 + 0.01 * frame_id, 5.85 + 0.01 * frame_id, 5.75 + 0.01 * frame_id,
                  5.4 + 0.01 * frame_id, 5.2 + 0.01 * frame_id};

    for (size_t i = 0; i < fr.atoms.x.size(); ++i) {
        simio::Vec3d r{{fr.atoms.x[i], fr.atoms.y[i], fr.atoms.z[i]}};
        fr.pbc.wrap_pos3(r);
        fr.atoms.x[i] = r.v[0];
        fr.atoms.y[i] = r.v[1];
        fr.atoms.z[i] = r.v[2];
    }
}

}  // namespace

int main() {
    simio::Topology topo = make_topology();

    simio::Frame fr;
    fr.pbc = simio::Pbc3D(5.0, 5.0, 10.0);
    fr.atoms.q = {
        -0.834, 0.417, 0.417,
        -0.834, 0.417, 0.417,
        1.0,
        -1.0,
    };

    simio::layered::Layer0Config l0_cfg;
    l0_cfg.nx = 10;
    l0_cfg.nz = 20;
    l0_cfg.x_center_wrapped = 2.5;
    l0_cfg.z_center_wrapped = 5.0;
    l0_cfg.channel_length_x = 5.0;
    l0_cfg.channel_height_z = 10.0;
    l0_cfg.gate_left_rel_z = -1.0;
    l0_cfg.gate_right_rel_z = 1.0;

    simio::layered::Layer1Config l1_cfg;
    l1_cfg.time_key_uses_ref_for_water = false;

    simio::layered::Pipeline pipe(l0_cfg, l1_cfg);

    simio::properties::DensityXZConfig density_cfg;
    density_cfg.selection = simio::properties::MoleculeSelection::All;
    density_cfg.roi.x_min = 1.0;
    density_cfg.roi.x_max = 4.0;
    density_cfg.roi.z_min = 4.6;
    density_cfg.roi.z_max = 6.6;
    density_cfg.roi.x_mode = simio::properties::DensityXMode::ChannelWindow;
    density_cfg.nx = 12;
    density_cfg.nz = 8;
    auto density = std::make_unique<simio::properties::DensityXZProperty>(density_cfg);
    simio::properties::DensityXZProperty* density_ptr = density.get();
    pipe.add_property(std::move(density));

    auto dipole = std::make_unique<simio::properties::DipoleXZProperty>();
    simio::properties::DipoleXZProperty* dipole_ptr = dipole.get();
    pipe.add_property(std::move(dipole));

    auto gate_center = std::make_unique<simio::properties::GatingCenterPlaneProperty>();
    simio::properties::GatingCenterPlaneProperty* gate_center_ptr = gate_center.get();
    pipe.add_property(std::move(gate_center));

    auto gate_wrap = std::make_unique<simio::properties::GatingPbcWrapProperty>();
    simio::properties::GatingPbcWrapProperty* gate_wrap_ptr = gate_wrap.get();
    pipe.add_property(std::move(gate_wrap));

    for (int f = 0; f < 3; ++f) {
        fill_frame(fr, f);
        pipe.process_frame(topo, fr);

        const auto& l0 = pipe.layer0_data();
        const auto& l1 = pipe.layer1_data();
        const auto& density_f = density_ptr->frames().back();
        const auto& dipole_f = dipole_ptr->frames().back();
        const auto& center_f = gate_center_ptr->frames().back();
        const auto& wrap_f = gate_wrap_ptr->frames().back();

        std::cout << "[frame " << f << "]"
                  << " mol0 bin(ix,iz)=(" << l0.ix[0] << "," << l0.iz[0] << ")"
                  << " region=" << static_cast<int>(l0.region[0])
                  << " vel0=(" << l1.velocity[0].v[0] << "," << l1.velocity[0].v[1] << ","
                  << l1.velocity[0].v[2] << ") density(selected/binned/oob)="
                  << density_f.selected_total << "/" << density_f.binned_total << "/"
                  << density_f.bin_oob_count
                  << " dipole_mean_mu=" << dipole_f.mean_mu
                  << " gate_center(NL,NR,dN)=(" << center_f.n_left << "," << center_f.n_right
                  << "," << center_f.dn << ") gate_wrap(NL,NR,dN)=(" << wrap_f.n_left << ","
                  << wrap_f.n_right << "," << wrap_f.dn << ")\n";
    }

    pipe.finalize();

    const auto& center_total = gate_center_ptr->cumulative();
    const auto& wrap_total = gate_wrap_ptr->cumulative();
    std::cout << "[gating cumulative] center(NL,NR,dN)=(" << center_total.n_left << ","
              << center_total.n_right << "," << center_total.dn << ")"
              << " wrap(NL,NR,dN)=(" << wrap_total.n_left << "," << wrap_total.n_right << ","
              << wrap_total.dn << ")\n";

    return 0;
}
