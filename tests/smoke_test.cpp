#include <cassert>
#include <vector>

#include "simio/simio.hpp"

int main() {
    simio::Topology topo;
    topo.mols.push_back(simio::MolSpan{0, 3, simio::MolType::Water});
    topo.build_type_lists();

    simio::Frame fr;
    fr.pbc = simio::Pbc3D(5.0, 5.0, 5.0);
    fr.atoms.x = {1.0, 1.1, 0.9};
    fr.atoms.y = {1.0, 1.0, 1.1};
    fr.atoms.z = {1.0, 1.1, 1.0};

    std::vector<simio::MolState> ms(topo.mols.size());
    auto pipe = simio::make_default_pipeline(1, 0.5);
    pipe.process_frame(topo, fr, ms);

    assert(ms[0].cache.flags & simio::MolCache::HAS_KEY);
    assert(!fr.grid.cells.empty());

    return 0;
}
