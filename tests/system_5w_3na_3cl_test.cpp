#include <cassert>
#include <cmath>
#include <vector>

#include "simio/simio.hpp"

namespace {

constexpr int kWaterCount = 5;
constexpr int kNaCount = 3;
constexpr int kClCount = 3;
constexpr int kThreads = 4;
constexpr int kFrames = 4;

int total_molecules() { return kWaterCount + kNaCount + kClCount; }
int total_atoms() { return kWaterCount * 3 + kNaCount + kClCount; }

simio::Topology build_topology() {
    simio::Topology topo;
    int atom_cursor = 0;

    for (int i = 0; i < kWaterCount; ++i) {
        topo.mols.push_back(simio::MolSpan{atom_cursor, 3, simio::MolType::Water});
        atom_cursor += 3;
    }
    for (int i = 0; i < kNaCount; ++i) {
        topo.mols.push_back(simio::MolSpan{atom_cursor, 1, simio::MolType::Cation});
        atom_cursor += 1;
    }
    for (int i = 0; i < kClCount; ++i) {
        topo.mols.push_back(simio::MolSpan{atom_cursor, 1, simio::MolType::Anion});
        atom_cursor += 1;
    }

    assert(atom_cursor == total_atoms());
    topo.build_type_lists();
    return topo;
}

void assign_charges(simio::Frame& fr) {
    fr.atoms.q.assign((size_t)total_atoms(), 0.0);

    for (int w = 0; w < kWaterCount; ++w) {
        const int O = w * 3 + 0;
        const int H1 = w * 3 + 1;
        const int H2 = w * 3 + 2;
        fr.atoms.q[(size_t)O] = -0.834;
        fr.atoms.q[(size_t)H1] = 0.417;
        fr.atoms.q[(size_t)H2] = 0.417;
    }

    const int na_start = kWaterCount * 3;
    const int cl_start = na_start + kNaCount;
    for (int i = 0; i < kNaCount; ++i) fr.atoms.q[(size_t)(na_start + i)] = 1.0;
    for (int i = 0; i < kClCount; ++i) fr.atoms.q[(size_t)(cl_start + i)] = -1.0;
}

void fill_frame_positions(simio::Frame& fr, int frame_id) {
    fr.atoms.x.assign((size_t)total_atoms(), 0.0);
    fr.atoms.y.assign((size_t)total_atoms(), 0.0);
    fr.atoms.z.assign((size_t)total_atoms(), 0.0);

    auto small_drift = [frame_id](int atom_id, double scale) {
        const double f = static_cast<double>(frame_id);
        const double a = static_cast<double>(atom_id);
        return scale * std::sin(0.30 * f + 0.17 * a);
    };

    for (int w = 0; w < kWaterCount; ++w) {
        const double bx = 0.8 + 0.75 * w;
        const double by = 0.9 + 0.35 * (w % 3);
        const double bz = 1.0 + 0.9 * (w / 2);

        const int O = w * 3 + 0;
        const int H1 = w * 3 + 1;
        const int H2 = w * 3 + 2;

        const double dx = 0.03 * frame_id + small_drift(O, 0.015);
        const double dy = 0.02 * frame_id + small_drift(O, 0.012);
        const double dz = 0.01 * frame_id + small_drift(O, 0.010);

        fr.atoms.x[(size_t)O] = bx + dx;
        fr.atoms.y[(size_t)O] = by + dy;
        fr.atoms.z[(size_t)O] = bz + dz;

        fr.atoms.x[(size_t)H1] = bx + 0.095 + dx;
        fr.atoms.y[(size_t)H1] = by + dy;
        fr.atoms.z[(size_t)H1] = bz + 0.028 + dz;

        fr.atoms.x[(size_t)H2] = bx - 0.031 + dx;
        fr.atoms.y[(size_t)H2] = by + 0.089 + dy;
        fr.atoms.z[(size_t)H2] = bz - 0.020 + dz;
    }

    const int na_start = kWaterCount * 3;
    for (int i = 0; i < kNaCount; ++i) {
        const int a = na_start + i;
        fr.atoms.x[(size_t)a] = 1.2 + 1.4 * i + 0.04 * frame_id + small_drift(a, 0.02);
        fr.atoms.y[(size_t)a] = 3.4 + 0.3 * i + 0.02 * frame_id + small_drift(a, 0.01);
        fr.atoms.z[(size_t)a] = 2.2 + 0.2 * i + 0.01 * frame_id + small_drift(a, 0.01);
    }

    const int cl_start = na_start + kNaCount;
    for (int i = 0; i < kClCount; ++i) {
        const int a = cl_start + i;
        fr.atoms.x[(size_t)a] = 0.7 + 1.5 * i + 0.03 * frame_id + small_drift(a, 0.02);
        fr.atoms.y[(size_t)a] = 4.2 - 0.4 * i + 0.02 * frame_id + small_drift(a, 0.01);
        fr.atoms.z[(size_t)a] = 7.5 - 0.5 * i + 0.01 * frame_id + small_drift(a, 0.01);
    }

    for (int a = 0; a < total_atoms(); ++a) {
        simio::Vec3d r{{fr.atoms.x[(size_t)a], fr.atoms.y[(size_t)a], fr.atoms.z[(size_t)a]}};
        fr.pbc.wrap_pos3(r);
        fr.atoms.x[(size_t)a] = r.v[0];
        fr.atoms.y[(size_t)a] = r.v[1];
        fr.atoms.z[(size_t)a] = r.v[2];
    }
}

int count_molecules_in_grid(const simio::Frame::MolGrid& grid) {
    int total = 0;
    for (const auto& bucket : grid.cells) total += (int)bucket.size();
    return total;
}

void validate_grid_membership(const simio::Frame::MolGrid& grid, int nmol) {
    std::vector<int> seen((size_t)nmol, 0);
    for (const auto& bucket : grid.cells) {
        for (int mid : bucket) {
            assert(mid >= 0 && mid < nmol);
            ++seen[(size_t)mid];
        }
    }
    for (int c : seen) assert(c == 1);
}

}  // namespace

int main() {
    simio::Topology topo = build_topology();

    simio::Frame fr;
    fr.pbc = simio::Pbc3D(5.0, 5.0, 10.0);
    assign_charges(fr);

    std::vector<simio::MolState> ms((size_t)total_molecules());
    auto pipe = simio::make_default_pipeline(kThreads, 0.5);

    for (int f = 0; f < kFrames; ++f) {
        fr.step = f;
        fr.time_ps = 0.5 * f;
        fill_frame_positions(fr, f);

        pipe.process_frame(topo, fr, ms);

        int key_count = 0;
        for (const auto& st : ms) {
            if (st.cache.flags & simio::MolCache::HAS_KEY) ++key_count;
            if (st.cache.flags & simio::MolCache::HAS_KEY) {
                assert(simio::is_wrapped3_box(st.cache.key_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2]));
            }
            if (st.cache.flags & simio::MolCache::HAS_REF) {
                assert(simio::is_wrapped3_box(st.cache.ref_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2]));
            }
            assert(!st.time.step_anomaly);
        }
        const int grid_count = count_molecules_in_grid(fr.grid);

        assert(key_count == total_molecules());
        assert(grid_count == total_molecules());
        validate_grid_membership(fr.grid, total_molecules());
    }

    return 0;
}
