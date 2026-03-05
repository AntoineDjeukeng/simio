#include <cassert>
#include <cmath>
#include <vector>

#include "simio/simio.hpp"

namespace {

constexpr double kTol = 1e-9;

}  // namespace

int main() {
    simio::Topology topo;
    topo.mols.push_back(simio::MolSpan{0, 1, simio::MolType::Cation});
    topo.build_type_lists();

    simio::Frame fr;
    fr.pbc = simio::Pbc3D(5.0, 5.0, 5.0);
    fr.atoms.x.assign(1, 0.0);
    fr.atoms.y.assign(1, 1.0);
    fr.atoms.z.assign(1, 1.0);

    std::vector<simio::MolState> ms(1);
    auto pipe = simio::make_default_pipeline(1, 0.5);

    fr.step = 0;
    fr.time_ps = 0.0;
    fr.atoms.x[0] = 4.9;
    pipe.process_frame(topo, fr, ms);

    fr.step = 1;
    fr.time_ps = 0.5;
    fr.atoms.x[0] = 0.1;
    pipe.process_frame(topo, fr, ms);

    const auto& st = ms[0];
    assert(st.cache.flags & simio::MolCache::HAS_KEY);
    assert(simio::is_wrapped3_box(st.cache.key_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2]));

    assert(st.time.pbc_correction_axes >= 1);
    assert(!st.time.step_anomaly);

    assert(std::abs(st.time.last_raw_dr.v[0] + 4.8) < kTol);
    assert(std::abs(st.time.last_dr.v[0] - 0.2) < kTol);
    assert(std::abs(st.time.last_step_norm - 0.2) < kTol);
    assert(std::abs(st.time.key_cont.v[0] - 5.1) < kTol);

    return 0;
}
