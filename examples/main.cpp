#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <vector>

#include "simio/simio.hpp"

#if !defined(SIMIO_TRACE)
#define SIMIO_TRACE 0
#endif

#if !defined(SIMIO_DEBUG)
#define SIMIO_DEBUG 0
#endif

namespace {

constexpr int kWaterCount = 5;
constexpr int kNaCount = 3;
constexpr int kClCount = 3;
constexpr int kThreads = 4;
constexpr int kFrames = 4;
constexpr double kGridCellSizeNm = 0.5;
constexpr bool kVerboseThreadTrace = (SIMIO_TRACE != 0);
constexpr std::int64_t kTlsCellWarnThreshold = 500000;

int total_molecules() { return kWaterCount + kNaCount + kClCount; }
int total_atoms() { return kWaterCount * 3 + kNaCount + kClCount; }

const char* mol_type_name(simio::MolType t) {
    switch (t) {
        case simio::MolType::Water:
            return "Water";
        case simio::MolType::Cation:
            return "Na";
        case simio::MolType::Anion:
            return "Cl";
        case simio::MolType::Other:
            return "Other";
        case simio::MolType::MOLTYPE_N:
            return "Invalid";
    }
    return "Invalid";
}

struct TypeSummary {
    int count = 0;
    int min_natoms = std::numeric_limits<int>::max();
    int max_natoms = std::numeric_limits<int>::min();

    void update(int natoms) {
        ++count;
        if (natoms < min_natoms) min_natoms = natoms;
        if (natoms > max_natoms) max_natoms = natoms;
    }

    int min_or_zero() const { return count > 0 ? min_natoms : 0; }
    int max_or_zero() const { return count > 0 ? max_natoms : 0; }
};

struct TopologySummary {
    int nmol_total = 0;
    TypeSummary water;
    TypeSummary na;
    TypeSummary cl;
    TypeSummary other;
};

TopologySummary summarize_topology(const simio::Topology& topo) {
    TopologySummary s;
    s.nmol_total = (int)topo.mols.size();

    for (const auto& m : topo.mols) {
        if (m.type == simio::MolType::Water) {
            s.water.update(m.natoms);
        } else if (m.type == simio::MolType::Cation) {
            s.na.update(m.natoms);
        } else if (m.type == simio::MolType::Anion) {
            s.cl.update(m.natoms);
        } else {
            s.other.update(m.natoms);
        }
    }

    return s;
}

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

void resize_atoms(simio::Frame& fr) {
    fr.atoms.x.assign((size_t)total_atoms(), 0.0);
    fr.atoms.y.assign((size_t)total_atoms(), 0.0);
    fr.atoms.z.assign((size_t)total_atoms(), 0.0);
}

void fill_frame_positions(simio::Frame& fr, int frame_id) {
    resize_atoms(fr);

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

struct GridValidation {
    int entries = 0;
    int duplicates = 0;
    int missing = 0;
    int out_of_range = 0;
};

struct GridBucketStats {
    int max_bucket = 0;
    int nonempty_buckets = 0;
    double mean_nonempty_bucket = 0.0;
};

GridValidation validate_grid(const simio::Frame::MolGrid& grid, int nmol) {
    GridValidation out;
    std::vector<int> seen((size_t)nmol, 0);

    for (const auto& bucket : grid.cells) {
        out.entries += (int)bucket.size();
        for (int mid : bucket) {
            if (mid < 0 || mid >= nmol) {
                ++out.out_of_range;
                continue;
            }
            ++seen[(size_t)mid];
        }
    }

    for (int c : seen) {
        if (c == 0) ++out.missing;
        if (c > 1) out.duplicates += (c - 1);
    }

    return out;
}

GridBucketStats compute_grid_bucket_stats(const simio::Frame::MolGrid& grid) {
    GridBucketStats out;
    int sum_nonempty = 0;

    for (const auto& bucket : grid.cells) {
        const int sz = (int)bucket.size();
        if (sz > out.max_bucket) out.max_bucket = sz;
        if (sz > 0) {
            ++out.nonempty_buckets;
            sum_nonempty += sz;
        }
    }

    if (out.nonempty_buckets > 0) {
        out.mean_nonempty_bucket = (double)sum_nonempty / (double)out.nonempty_buckets;
    }
    return out;
}

struct StepStats {
    double min_step = 0.0;
    double max_step = 0.0;
    double mean_step = 0.0;
    int pbc_anom = 0;
    int pbc_cross_count = 0;
    double max_abs_dr_raw_component = 0.0;
    double max_abs_dr_min_component = 0.0;
};

StepStats compute_step_stats(const std::vector<simio::MolState>& ms) {
    StepStats out;
    if (ms.empty()) return out;

    double sum = 0.0;
    out.min_step = std::numeric_limits<double>::max();

    for (const auto& st : ms) {
        const double step = st.time.last_step_norm;
        if (step < out.min_step) out.min_step = step;
        if (step > out.max_step) out.max_step = step;
        sum += step;
        if (st.time.step_anomaly) ++out.pbc_anom;
        if (st.time.pbc_correction_axes > 0) ++out.pbc_cross_count;

        for (int a = 0; a < 3; ++a) {
            const double raw_abs = std::abs(st.time.last_raw_dr.v[a]);
            const double min_abs = std::abs(st.time.last_dr.v[a]);
            if (raw_abs > out.max_abs_dr_raw_component) out.max_abs_dr_raw_component = raw_abs;
            if (min_abs > out.max_abs_dr_min_component) out.max_abs_dr_min_component = min_abs;
        }
    }

    out.mean_step = sum / (double)ms.size();
    return out;
}

struct WrappedViolations {
    int key = 0;
    int ref = 0;
};

WrappedViolations count_wrapped_violations(const simio::Frame& fr, const std::vector<simio::MolState>& ms) {
    WrappedViolations out;

    for (const auto& st : ms) {
        if ((st.cache.flags & simio::MolCache::HAS_KEY) &&
            !simio::is_wrapped3_box(st.cache.key_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2])) {
            ++out.key;
        }
        if ((st.cache.flags & simio::MolCache::HAS_REF) &&
            !simio::is_wrapped3_box(st.cache.ref_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2])) {
            ++out.ref;
        }
    }

    return out;
}

void print_vec3(const simio::Vec3d& v) {
    std::cout << "(" << v.v[0] << ", " << v.v[1] << ", " << v.v[2] << ")";
}

void print_topology_report(const TopologySummary& ts) {
    std::cout << "[topology] nmol=" << ts.nmol_total << " water=" << ts.water.count
              << " na=" << ts.na.count << " cl=" << ts.cl.count << " other=" << ts.other.count
              << "\n";
    std::cout << "[topology] natoms: water=3 (min=" << ts.water.min_or_zero()
              << " max=" << ts.water.max_or_zero() << ")"
              << " na=1 (min=" << ts.na.min_or_zero() << " max=" << ts.na.max_or_zero() << ")"
              << " cl=1 (min=" << ts.cl.min_or_zero() << " max=" << ts.cl.max_or_zero() << ")\n";
    std::cout << "[threads] nthreads=" << kThreads << " loop=stride(i=tid; i<n; i+=nthreads)\n";
    std::cout << "[keys] grid_key=COM  time_key=COM  ref_key=O\n";
}

#if SIMIO_TRACE
std::mutex g_trace_mu;

struct TraceMolPropertyB : simio::Property {
    const char* name() const override { return "B_neighbor_probe"; }
    simio::Scope scope() const override { return simio::Scope::Mol; }
    simio::Stage stage() const override { return simio::Stage::Neighbor; }

    void run_mol(const simio::Topology& topo, simio::Frame& fr, int mol_id,
                 std::vector<simio::MolState>&, int tid) override {
        const auto& span = topo.mols[(size_t)mol_id];
        std::lock_guard<std::mutex> lock(g_trace_mu);
        std::cout << "[trace][frame " << fr.step << "][T" << tid << "] MolProperty '" << name()
                  << "' -> molecule " << mol_id << " (" << mol_type_name(span.type) << ")\n";
    }
};

struct TraceAtomPropertyC : simio::Property {
    const char* name() const override { return "C_atom_probe"; }
    simio::Scope scope() const override { return simio::Scope::Atom; }
    simio::Stage stage() const override { return simio::Stage::Neighbor; }

    void run_atom(const simio::Topology&, simio::Frame& fr, int atom_id,
                  std::vector<simio::MolState>&, int tid) override {
        std::lock_guard<std::mutex> lock(g_trace_mu);
        std::cout << "[trace][frame " << fr.step << "][T" << tid << "] AtomProperty '" << name()
                  << "' -> atom " << atom_id << "\n";
    }
};
#endif

void print_time_probe(int mol_id, const std::vector<simio::MolState>& ms) {
    const auto& st = ms[(size_t)mol_id];
    std::cout << "[dbg][mol " << mol_id << "] key_wrapped=";
    print_vec3(st.cache.key_wrapped);
    std::cout << " prev_before=";
    print_vec3(st.time.prev_key_before_update);
    std::cout << " prev_after=";
    print_vec3(st.time.prev_key_after_update);
    std::cout << " dr_raw=";
    print_vec3(st.time.last_raw_dr);
    std::cout << " dr_min_image=";
    print_vec3(st.time.last_dr);
    std::cout << " key_cont=";
    print_vec3(st.time.key_cont);
    std::cout << "\n";
}

}  // namespace

int main() {
    std::cout << "=== simio pipeline check ===\n";

    simio::Topology topo = build_topology();
    const TopologySummary ts = summarize_topology(topo);
    print_topology_report(ts);

    SIMIO_DEBUG_ASSERT(ts.nmol_total == total_molecules());
    SIMIO_DEBUG_ASSERT(ts.water.count == kWaterCount && ts.water.min_or_zero() == 3 && ts.water.max_or_zero() == 3);
    SIMIO_DEBUG_ASSERT(ts.na.count == kNaCount && ts.na.min_or_zero() == 1 && ts.na.max_or_zero() == 1);
    SIMIO_DEBUG_ASSERT(ts.cl.count == kClCount && ts.cl.min_or_zero() == 1 && ts.cl.max_or_zero() == 1);

    simio::Frame fr;
    fr.pbc = simio::Pbc3D(5.0, 5.0, 10.0);
    assign_charges(fr);
    fr.grid.init(fr.pbc, kGridCellSizeNm);

    const std::int64_t tls_slots = (std::int64_t)kThreads * (std::int64_t)fr.grid.ncells();
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[box] L=(" << fr.pbc.L[0] << ", " << fr.pbc.L[1] << ", " << fr.pbc.L[2]
              << ") cell=(" << fr.grid.cx << ", " << fr.grid.cy << ", " << fr.grid.cz << ")"
              << " n=(" << fr.grid.nx << ", " << fr.grid.ny << ", " << fr.grid.nz << ")"
              << " ncells=" << fr.grid.ncells() << " tls_vectors=" << tls_slots << "\n";
    if (tls_slots > kTlsCellWarnThreshold) {
        std::cout << "[warn] large TLS grid footprint nthreads*ncells=" << tls_slots << "\n";
    }
    std::cout << "[neighbor] counting=central-only\n";

    std::vector<simio::MolState> ms((size_t)total_molecules());
    auto pipe = simio::make_default_pipeline(kThreads, kGridCellSizeNm);

#if SIMIO_TRACE
    if (kVerboseThreadTrace) {
        pipe.add(std::make_unique<TraceMolPropertyB>());
        pipe.add(std::make_unique<TraceAtomPropertyC>());
    }
#endif

    for (int f = 0; f < kFrames; ++f) {
        fr.step = f;
        fr.time_ps = 0.5 * f;

        fill_frame_positions(fr, f);
        pipe.process_frame(topo, fr, ms);

        int keys_ok = 0;
        for (const auto& st : ms) {
            if (st.cache.flags & simio::MolCache::HAS_KEY) ++keys_ok;
        }

        const WrappedViolations wr = count_wrapped_violations(fr, ms);
        const GridValidation gv = validate_grid(fr.grid, (int)ms.size());
        const GridBucketStats gb = compute_grid_bucket_stats(fr.grid);
        const StepStats ss = compute_step_stats(ms);

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "[frame " << f << "] step=" << fr.step << " time_ps=" << fr.time_ps << "\n";
        std::cout << "  keys=" << keys_ok << "/" << total_molecules() << " key_range_viol=" << wr.key
                  << " ref_range_viol=" << wr.ref << "\n";
        std::cout << "  grid: cell=" << kGridCellSizeNm << " nx=" << fr.grid.nx << " ny=" << fr.grid.ny
                  << " nz=" << fr.grid.nz << " ncells=" << fr.grid.ncells()
                  << " entries=" << gv.entries << " dup=" << gv.duplicates
                  << " miss=" << gv.missing << " out_of_range=" << gv.out_of_range
                  << " max_bucket=" << gb.max_bucket
                  << " mean_nonempty_bucket=" << gb.mean_nonempty_bucket << "\n";
        std::cout << "  step(|dr|): min=" << ss.min_step << " max=" << ss.max_step
                  << " mean=" << ss.mean_step << " pbc_anom=" << ss.pbc_anom
                  << " pbc_cross=" << ss.pbc_cross_count
                  << " max|dr_raw_comp|=" << ss.max_abs_dr_raw_component
                  << " max|dr_min_comp|=" << ss.max_abs_dr_min_component << "\n";

#if SIMIO_DEBUG
        print_time_probe(0, ms);
        print_time_probe(kWaterCount, ms);
        print_time_probe(kWaterCount + kNaCount, ms);
#endif

        const bool frame_ok =
            (keys_ok == total_molecules()) && (wr.key == 0) && (wr.ref == 0) &&
            (gv.entries == total_molecules()) && (gv.duplicates == 0) && (gv.missing == 0) &&
            (gv.out_of_range == 0) && (ss.pbc_anom == 0);

        if (!frame_ok) {
            std::cerr << "[error] frame " << f << " failed validation checks\n";
            return 1;
        }
    }

    std::cout << "All frames processed successfully.\n";
    return 0;
}
