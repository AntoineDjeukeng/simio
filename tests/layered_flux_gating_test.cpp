#include <cassert>
#include <vector>

#include "simio/layered_pipeline.hpp"
#include "simio/properties/GatingCenterPlane.hpp"
#include "simio/properties/GatingPbcWrap.hpp"

int main() {
    struct RunResult {
        std::vector<simio::properties::GatingFrameCount> center_frames;
        std::vector<simio::properties::GatingFrameCount> wrap_frames;
        simio::properties::GatingFrameCount center_total;
        simio::properties::GatingFrameCount wrap_total;
    };

    auto run_case = [](const std::vector<double>& xs) -> RunResult {
        simio::Topology topo;
        topo.mols.push_back(simio::MolSpan{0, 1, simio::MolType::Cation});
        topo.build_type_lists();

        simio::Frame fr;
        fr.pbc = simio::Pbc3D(5.0, 5.0, 5.0);
        fr.atoms.x.assign(1, 0.0);
        fr.atoms.y.assign(1, 1.0);
        fr.atoms.z.assign(1, 1.0);

        simio::layered::Layer0Config l0;
        l0.nx = 10;
        l0.nz = 10;
        l0.x_center_wrapped = 2.5;
        l0.z_center_wrapped = 2.5;
        l0.channel_length_x = 5.0;
        l0.channel_height_z = 5.0;

        simio::layered::Layer1Config l1;
        simio::layered::Pipeline pipe(l0, l1);

        auto center = std::make_unique<simio::properties::GatingCenterPlaneProperty>();
        simio::properties::GatingCenterPlaneProperty* center_ptr = center.get();
        pipe.add_property(std::move(center));

        auto wrap = std::make_unique<simio::properties::GatingPbcWrapProperty>();
        simio::properties::GatingPbcWrapProperty* wrap_ptr = wrap.get();
        pipe.add_property(std::move(wrap));

        for (size_t i = 0; i < xs.size(); ++i) {
            fr.step = (int64_t)i;
            fr.time_ps = (double)i;
            fr.atoms.x[0] = xs[i];
            pipe.process_frame(topo, fr);
        }

        return RunResult{
            center_ptr->frames(),
            wrap_ptr->frames(),
            center_ptr->cumulative(),
            wrap_ptr->cumulative(),
        };
    };

    // Case 1: seam crossing only (4.9 -> 0.1 -> 4.9).
    const RunResult seam_case = run_case({4.9, 0.1, 4.9});
    assert(seam_case.center_frames.size() == 3);
    assert(seam_case.wrap_frames.size() == 3);

    assert(seam_case.center_frames[0].n_left == 0);
    assert(seam_case.center_frames[0].n_right == 0);
    assert(seam_case.wrap_frames[0].n_left == 0);
    assert(seam_case.wrap_frames[0].n_right == 0);

    assert(seam_case.center_frames[1].n_left == 0);
    assert(seam_case.center_frames[1].n_right == 0);
    assert(seam_case.center_frames[1].dn == 0);
    assert(seam_case.wrap_frames[1].n_left == 0);
    assert(seam_case.wrap_frames[1].n_right == 1);
    assert(seam_case.wrap_frames[1].dn == 1);

    assert(seam_case.center_frames[2].n_left == 0);
    assert(seam_case.center_frames[2].n_right == 0);
    assert(seam_case.center_frames[2].dn == 0);
    assert(seam_case.wrap_frames[2].n_left == 1);
    assert(seam_case.wrap_frames[2].n_right == 0);
    assert(seam_case.wrap_frames[2].dn == -1);

    assert(seam_case.center_total.n_left == 0);
    assert(seam_case.center_total.n_right == 0);
    assert(seam_case.center_total.dn == 0);
    assert(seam_case.wrap_total.n_left == 1);
    assert(seam_case.wrap_total.n_right == 1);
    assert(seam_case.wrap_total.dn == 0);

    // Case 2: center crossing only (2.4 -> 2.6 -> 2.4).
    const RunResult center_case = run_case({2.4, 2.6, 2.4});
    assert(center_case.center_frames.size() == 3);
    assert(center_case.wrap_frames.size() == 3);

    assert(center_case.wrap_frames[1].n_left == 0);
    assert(center_case.wrap_frames[1].n_right == 0);
    assert(center_case.wrap_frames[2].n_left == 0);
    assert(center_case.wrap_frames[2].n_right == 0);
    assert(center_case.wrap_total.n_left == 0);
    assert(center_case.wrap_total.n_right == 0);
    assert(center_case.wrap_total.dn == 0);

    assert(center_case.center_frames[1].n_left == 0);
    assert(center_case.center_frames[1].n_right == 1);
    assert(center_case.center_frames[1].dn == 1);
    assert(center_case.center_frames[2].n_left == 1);
    assert(center_case.center_frames[2].n_right == 0);
    assert(center_case.center_frames[2].dn == -1);
    assert(center_case.center_total.n_left == 1);
    assert(center_case.center_total.n_right == 1);
    assert(center_case.center_total.dn == 0);

    return 0;
}
