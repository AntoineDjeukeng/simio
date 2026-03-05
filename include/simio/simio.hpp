#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <memory>
#include <thread>
#include <utility>
#include <vector>

#if !defined(SIMIO_DEBUG)
#define SIMIO_DEBUG 0
#endif

#if SIMIO_DEBUG
#define SIMIO_DEBUG_ASSERT(expr) assert(expr)
#else
#define SIMIO_DEBUG_ASSERT(expr) ((void)0)
#endif

namespace simio {

struct Vec3d {
    double v[3]{0.0, 0.0, 0.0};
    double& operator[](int a) { return v[a]; }
    double operator[](int a) const { return v[a]; }
};

inline Vec3d operator+(const Vec3d& a, const Vec3d& b) {
    return Vec3d{{a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2]}};
}
inline Vec3d operator-(const Vec3d& a, const Vec3d& b) {
    return Vec3d{{a.v[0] - b.v[0], a.v[1] - b.v[1], a.v[2] - b.v[2]}};
}
inline Vec3d operator*(const Vec3d& a, double s) {
    return Vec3d{{a.v[0] * s, a.v[1] * s, a.v[2] * s}};
}
inline Vec3d operator/(const Vec3d& a, double s) {
    return Vec3d{{a.v[0] / s, a.v[1] / s, a.v[2] / s}};
}
inline double dot3(const Vec3d& a, const Vec3d& b) {
    return a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2];
}
inline double norm3(const Vec3d& a) { return std::sqrt(dot3(a, a)); }

inline bool is_wrapped_axis(double x, double L, double eps = 1e-10) {
    return (x >= -eps) && (x < L + eps);
}

inline bool is_wrapped3_box(const Vec3d& r, double Lx, double Ly, double Lz, double eps = 1e-10) {
    return is_wrapped_axis(r.v[0], Lx, eps) && is_wrapped_axis(r.v[1], Ly, eps) &&
           is_wrapped_axis(r.v[2], Lz, eps);
}

struct Pbc3D {
    double L[3] = {0.0, 0.0, 0.0};

    Pbc3D() = default;
    Pbc3D(double Lx, double Ly, double Lz) {
        L[0] = Lx;
        L[1] = Ly;
        L[2] = Lz;
    }

    bool valid_axis(int a) const { return (0 <= a && a < 3 && L[a] > 0.0); }

    double wrap_pos(int a, double x) const {
        const double La = L[a];
        x -= La * std::floor(x / La);
        if (x >= La) x -= La;
        if (x < 0.0) x += La;
        return x;
    }

    double wrap_delta(int a, double d) const {
        const double La = L[a];
        d -= La * std::round(d / La);
        return d;
    }

    void wrap_pos3(Vec3d& r) const {
        r.v[0] = wrap_pos(0, r.v[0]);
        r.v[1] = wrap_pos(1, r.v[1]);
        r.v[2] = wrap_pos(2, r.v[2]);
    }

    void wrap_delta3(Vec3d& dr) const {
        dr.v[0] = wrap_delta(0, dr.v[0]);
        dr.v[1] = wrap_delta(1, dr.v[1]);
        dr.v[2] = wrap_delta(2, dr.v[2]);
    }

    Vec3d min_image(const Vec3d& dr) const {
        Vec3d out = dr;
        wrap_delta3(out);
        return out;
    }
};

inline Vec3d pbc_unwrap_to_ref(const Pbc3D& pbc, const Vec3d& ref, const Vec3d& x) {
    Vec3d out = x;
    out.v[0] = ref.v[0] + pbc.wrap_delta(0, x.v[0] - ref.v[0]);
    out.v[1] = ref.v[1] + pbc.wrap_delta(1, x.v[1] - ref.v[1]);
    out.v[2] = ref.v[2] + pbc.wrap_delta(2, x.v[2] - ref.v[2]);
    return out;
}

inline Vec3d water_com_pbc_unwrapped(const Pbc3D& pbc, const Vec3d& O, const Vec3d& H1,
                                     const Vec3d& H2) {
    const Vec3d H1u = pbc_unwrap_to_ref(pbc, O, H1);
    const Vec3d H2u = pbc_unwrap_to_ref(pbc, O, H2);

    constexpr double mO = 16.0;
    constexpr double mH = 1.0;
    constexpr double inv = 1.0 / (mO + 2.0 * mH);

    Vec3d com;
    com.v[0] = (mO * O.v[0] + mH * H1u.v[0] + mH * H2u.v[0]) * inv;
    com.v[1] = (mO * O.v[1] + mH * H1u.v[1] + mH * H2u.v[1]) * inv;
    com.v[2] = (mO * O.v[2] + mH * H1u.v[2] + mH * H2u.v[2]) * inv;
    return com;
}

inline Vec3d water_com_from_unwrapped_sites(const Vec3d& O, const Vec3d& H1u, const Vec3d& H2u) {

    constexpr double mO = 16.0;
    constexpr double mH = 1.0;
    constexpr double inv = 1.0 / (mO + 2.0 * mH);

    Vec3d com;
    com.v[0] = (mO * O.v[0] + mH * H1u.v[0] + mH * H2u.v[0]) * inv;
    com.v[1] = (mO * O.v[1] + mH * H1u.v[1] + mH * H2u.v[1]) * inv;
    com.v[2] = (mO * O.v[2] + mH * H1u.v[2] + mH * H2u.v[2]) * inv;
    return com;
}

enum class MolType : uint8_t { Water = 0, Cation = 1, Anion = 2, Other = 3, MOLTYPE_N = 4 };

struct MolSpan {
    int first = 0;
    int natoms = 0;
    MolType type = MolType::Other;
};

struct Topology {
    std::vector<MolSpan> mols;
    std::vector<int> mol_ids_by_type[(int)MolType::MOLTYPE_N];

    void build_type_lists() {
        for (auto& v : mol_ids_by_type) v.clear();
        for (int mid = 0; mid < (int)mols.size(); ++mid) {
            int t = (int)mols[(size_t)mid].type;
            if (t < 0) t = 0;
            if (t >= (int)MolType::MOLTYPE_N) t = (int)MolType::Other;
            mol_ids_by_type[t].push_back(mid);
        }
    }
};

struct AtomSoA {
    std::vector<double> x, y, z;
    std::vector<double> q;
    std::vector<double> m;

    size_t size() const { return x.size(); }

    Vec3d pos(int ai) const { return Vec3d{{x[(size_t)ai], y[(size_t)ai], z[(size_t)ai]}}; }
};

struct Frame {
    Pbc3D pbc;
    int64_t step = 0;
    double time_ps = 0.0;
    AtomSoA atoms;

    struct MolGrid {
        int nx = 1, ny = 1, nz = 1;
        double Lx = 1, Ly = 1, Lz = 1;
        double cx = 1, cy = 1, cz = 1;
        std::vector<std::vector<int>> cells;

        void init(const Pbc3D& pbc, double cell_size) {
            Lx = pbc.L[0];
            Ly = pbc.L[1];
            Lz = pbc.L[2];
            nx = std::max(1, (int)std::floor(Lx / cell_size));
            ny = std::max(1, (int)std::floor(Ly / cell_size));
            nz = std::max(1, (int)std::floor(Lz / cell_size));
            cx = Lx / nx;
            cy = Ly / ny;
            cz = Lz / nz;
            cells.assign((size_t)(nx * ny * nz), {});
        }

        int ncells() const { return nx * ny * nz; }

        void clear() {
            for (auto& v : cells) v.clear();
        }

        static int clampi(int v, int lo, int hi) {
            if (v < lo) return lo;
            if (v > hi) return hi;
            return v;
        }
        static int wrapi(int i, int n) {
            i %= n;
            if (i < 0) i += n;
            return i;
        }

        int cid(int ix, int iy, int iz) const {
            ix = wrapi(ix, nx);
            iy = wrapi(iy, ny);
            iz = wrapi(iz, nz);
            return (iz * ny + iy) * nx + ix;
        }

        int cell_id_from_wrapped(const Vec3d& r) const {
            SIMIO_DEBUG_ASSERT(is_wrapped3_box(r, Lx, Ly, Lz));

            int ix = wrapi((int)std::floor(r.v[0] / cx), nx);
            int iy = wrapi((int)std::floor(r.v[1] / cy), ny);
            int iz = wrapi((int)std::floor(r.v[2] / cz), nz);
            return cid(ix, iy, iz);
        }

        template <class F>
        void for_candidates(const Vec3d& r_wrapped, int ncell, F&& f) const {
            SIMIO_DEBUG_ASSERT(is_wrapped3_box(r_wrapped, Lx, Ly, Lz));

            int ix = wrapi((int)std::floor(r_wrapped.v[0] / cx), nx);
            int iy = wrapi((int)std::floor(r_wrapped.v[1] / cy), ny);
            int iz = wrapi((int)std::floor(r_wrapped.v[2] / cz), nz);

            for (int dz = -ncell; dz <= ncell; ++dz)
                for (int dy = -ncell; dy <= ncell; ++dy)
                    for (int dx = -ncell; dx <= ncell; ++dx) {
                        const auto& bucket = cells[(size_t)cid(ix + dx, iy + dy, iz + dz)];
                        for (int mid : bucket) f(mid);
                    }
        }

        template <class F>
        void for_candidates_box(const Vec3d& r_wrapped, int nx_cell, int ny_cell, int nz_cell,
                                F&& f) const {
            SIMIO_DEBUG_ASSERT(is_wrapped3_box(r_wrapped, Lx, Ly, Lz));

            int ix = wrapi((int)std::floor(r_wrapped.v[0] / cx), nx);
            int iy = wrapi((int)std::floor(r_wrapped.v[1] / cy), ny);
            int iz = wrapi((int)std::floor(r_wrapped.v[2] / cz), nz);

            for (int dz = -nz_cell; dz <= nz_cell; ++dz)
                for (int dy = -ny_cell; dy <= ny_cell; ++dy)
                    for (int dx = -nx_cell; dx <= nx_cell; ++dx) {
                        const auto& bucket = cells[(size_t)cid(ix + dx, iy + dy, iz + dz)];
                        for (int mid : bucket) f(mid);
                    }
        }
    } grid;
};

struct MolCache {
    uint32_t flags = 0;

    Vec3d key_wrapped{};
    Vec3d ref_wrapped{};
    Vec3d com_unwrapped{};
    Vec3d dipole{};
    std::array<Vec3d, 3> sites_u{};

    static constexpr uint32_t HAS_KEY = 1u << 0;
    static constexpr uint32_t HAS_REF = 1u << 1;
    static constexpr uint32_t HAS_SITES = 1u << 2;
    static constexpr uint32_t HAS_COM = 1u << 3;
    static constexpr uint32_t HAS_DIP = 1u << 4;
};

struct MolResults {
    int hba = 0, hbd = 0, hbdmp = 0, hbdhmp = 0;
    int ibc = 0, iba = 0;
    int bnc = 0, bna = 0;

    double mu_x = 0.0, mu_y = 0.0, mu_z = 0.0;
};

struct MolTimeState {
    bool has_prev = false;
    Vec3d prev_key_wrapped{};
    Vec3d prev_key_before_update{};
    Vec3d prev_key_after_update{};
    Vec3d key_cont{};
    Vec3d cum_disp{};
    double max_jump = 0.0;
    double last_step_norm = 0.0;
    Vec3d last_raw_dr{};
    Vec3d last_dr{};
    int pbc_correction_axes = 0;
    bool step_anomaly = false;

    int64_t NL = 0, NR = 0, dN = 0;
};

struct MolState {
    MolCache cache;
    MolResults res;
    MolTimeState time;
};

template <class Func>
inline void parallel_strided(int nthreads, int nitems, Func&& fn) {
    if (nthreads <= 1 || nitems <= 0) {
        for (int i = 0; i < nitems; ++i) fn(i, 0);
        return;
    }
    std::vector<std::thread> th;
    th.reserve((size_t)nthreads);
    for (int tid = 0; tid < nthreads; ++tid) {
        th.emplace_back([&, tid]() {
            for (int i = tid; i < nitems; i += nthreads) fn(i, tid);
        });
    }
    for (auto& t : th) t.join();
}

template <class Func>
inline void parallel_for_molids_strided(int nthreads, const std::vector<int>& mol_ids, Func&& fn) {
    if (nthreads <= 1 || mol_ids.empty()) {
        for (int i = 0; i < (int)mol_ids.size(); ++i) fn(mol_ids[(size_t)i], 0);
        return;
    }
    std::vector<std::thread> th;
    th.reserve((size_t)nthreads);
    for (int tid = 0; tid < nthreads; ++tid) {
        th.emplace_back([&, tid]() {
            for (int i = tid; i < (int)mol_ids.size(); i += nthreads) fn(mol_ids[(size_t)i], tid);
        });
    }
    for (auto& t : th) t.join();
}

enum class Scope { Frame, MolType, Mol, Atom };
enum class Stage { BeginFrame, IntrinsicAndGrid, Neighbor, Time, Reduce, EndFrame };

struct Property {
    virtual ~Property() = default;
    virtual const char* name() const = 0;
    virtual Scope scope() const = 0;
    virtual Stage stage() const = 0;

    virtual void run_frame(const Topology&, Frame&, std::vector<MolState>&) {}
    virtual void run_moltype(const Topology&, Frame&, MolType, std::vector<MolState>&) {}
    virtual void run_mol(const Topology&, Frame&, int mol_id, std::vector<MolState>&, int tid) {}
    virtual void run_atom(const Topology&, Frame&, int atom_id, std::vector<MolState>&, int tid) {}
};

inline void compute_key_and_intrinsics_for_mol(const Topology& topo, Frame& fr, int mol_id,
                                               MolState& st) {
    const MolSpan& m = topo.mols[(size_t)mol_id];
    MolCache& c = st.cache;

    if (m.type == MolType::Water) {
        assert(m.natoms >= 3);
        const int O = m.first + 0;
        const int H1 = m.first + 1;
        const int H2 = m.first + 2;

        Vec3d rO = fr.atoms.pos(O);
        Vec3d rH1 = fr.atoms.pos(H1);
        Vec3d rH2 = fr.atoms.pos(H2);

        fr.pbc.wrap_pos3(rO);
        c.ref_wrapped = rO;
        c.flags |= MolCache::HAS_REF;
        SIMIO_DEBUG_ASSERT(is_wrapped3_box(c.ref_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2]));

        const Vec3d H1u = pbc_unwrap_to_ref(fr.pbc, rO, rH1);
        const Vec3d H2u = pbc_unwrap_to_ref(fr.pbc, rO, rH2);
        c.sites_u[0] = rO;
        c.sites_u[1] = H1u;
        c.sites_u[2] = H2u;
        c.flags |= MolCache::HAS_SITES;

        Vec3d com_u = water_com_from_unwrapped_sites(rO, H1u, H2u);
        c.com_unwrapped = com_u;
        c.flags |= MolCache::HAS_COM;

        Vec3d key = com_u;
        fr.pbc.wrap_pos3(key);
        c.key_wrapped = key;
        c.flags |= MolCache::HAS_KEY;
        SIMIO_DEBUG_ASSERT(is_wrapped3_box(c.key_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2]));

        if (!fr.atoms.q.empty()) {
            auto q = [&](int ai) -> double { return fr.atoms.q[(size_t)ai]; };
            Vec3d mu{};
            mu = mu + (rO - com_u) * q(O);
            mu = mu + (H1u - com_u) * q(H1);
            mu = mu + (H2u - com_u) * q(H2);
            c.dipole = mu;
            c.flags |= MolCache::HAS_DIP;
            st.res.mu_x = mu.v[0];
            st.res.mu_y = mu.v[1];
            st.res.mu_z = mu.v[2];
        }
    } else if (m.type == MolType::Cation || m.type == MolType::Anion) {
        assert(m.natoms >= 1);
        Vec3d r = fr.atoms.pos(m.first);
        fr.pbc.wrap_pos3(r);
        c.key_wrapped = r;
        c.flags |= MolCache::HAS_KEY;
        c.ref_wrapped = r;
        c.flags |= MolCache::HAS_REF;
        SIMIO_DEBUG_ASSERT(is_wrapped3_box(c.key_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2]));
        SIMIO_DEBUG_ASSERT(is_wrapped3_box(c.ref_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2]));
    } else {
        if (m.natoms >= 1) {
            Vec3d r = fr.atoms.pos(m.first);
            fr.pbc.wrap_pos3(r);
            c.key_wrapped = r;
            c.flags |= MolCache::HAS_KEY;
            SIMIO_DEBUG_ASSERT(
                is_wrapped3_box(c.key_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2]));
        }
    }
}

struct BuildIntrinsicsAndGrid : Property {
    int nthreads = 1;
    double grid_cell_size_nm = 0.5;

    explicit BuildIntrinsicsAndGrid(int nthreads_, double cell_size_nm)
        : nthreads(std::max(1, nthreads_)), grid_cell_size_nm(cell_size_nm) {}

    const char* name() const override { return "BuildIntrinsicsAndGrid"; }
    Scope scope() const override { return Scope::Frame; }
    Stage stage() const override { return Stage::IntrinsicAndGrid; }

    void run_frame(const Topology& topo, Frame& fr, std::vector<MolState>& ms) override {
        if (fr.grid.cells.empty()) fr.grid.init(fr.pbc, grid_cell_size_nm);

        const int ncells = fr.grid.ncells();
        const int nmol = (int)topo.mols.size();

        using Buckets = std::vector<std::vector<int>>;
        std::vector<Buckets> tls((size_t)nthreads);
        for (int tid = 0; tid < nthreads; ++tid) tls[(size_t)tid].assign((size_t)ncells, {});

        parallel_strided(nthreads, nmol, [&](int mol_id, int tid) {
            MolState& st = ms[(size_t)mol_id];
            compute_key_and_intrinsics_for_mol(topo, fr, mol_id, st);

            if (!(st.cache.flags & MolCache::HAS_KEY)) return;

            Vec3d key = st.cache.key_wrapped;
            fr.pbc.wrap_pos3(key);
            const int cid = fr.grid.cell_id_from_wrapped(key);
            tls[(size_t)tid][(size_t)cid].push_back(mol_id);
        });

        fr.grid.clear();
        for (int cid = 0; cid < ncells; ++cid) {
            size_t total = 0;
            for (int tid = 0; tid < nthreads; ++tid) total += tls[(size_t)tid][(size_t)cid].size();
            fr.grid.cells[(size_t)cid].reserve(total);
            for (int tid = 0; tid < nthreads; ++tid) {
                auto& src = tls[(size_t)tid][(size_t)cid];
                auto& dst = fr.grid.cells[(size_t)cid];
                dst.insert(dst.end(), src.begin(), src.end());
            }
        }
    }
};

struct BeginFrameClear : Property {
    const char* name() const override { return "BeginFrameClear"; }
    Scope scope() const override { return Scope::Mol; }
    Stage stage() const override { return Stage::BeginFrame; }

    void run_mol(const Topology&, Frame&, int mol_id, std::vector<MolState>& ms, int) override {
        ms[(size_t)mol_id].cache.flags = 0;
        ms[(size_t)mol_id].res = MolResults{};
    }
};

struct TimeDisplacement : Property {
    const char* name() const override { return "TimeDisplacement"; }
    Scope scope() const override { return Scope::Mol; }
    Stage stage() const override { return Stage::Time; }

    void run_mol(const Topology&, Frame& fr, int mol_id, std::vector<MolState>& ms, int) override {
        MolState& st = ms[(size_t)mol_id];

        if (!(st.cache.flags & MolCache::HAS_KEY)) return;
        Vec3d key = st.cache.key_wrapped;
        fr.pbc.wrap_pos3(key);

        if (!st.time.has_prev) {
            st.time.has_prev = true;
            st.time.prev_key_wrapped = key;
            st.time.prev_key_before_update = key;
            st.time.prev_key_after_update = key;
            st.time.key_cont = key;
            st.time.cum_disp = Vec3d{};
            st.time.last_step_norm = 0.0;
            st.time.last_raw_dr = Vec3d{};
            st.time.last_dr = Vec3d{};
            st.time.pbc_correction_axes = 0;
            st.time.step_anomaly = false;
            return;
        }

        Vec3d dr_raw = key - st.time.prev_key_wrapped;
        Vec3d dr = fr.pbc.min_image(dr_raw);
        st.time.prev_key_before_update = st.time.prev_key_wrapped;
        st.time.prev_key_wrapped = key;
        st.time.prev_key_after_update = st.time.prev_key_wrapped;

        st.time.key_cont = st.time.key_cont + dr;
        st.time.cum_disp = st.time.cum_disp + dr;
        st.time.last_raw_dr = dr_raw;
        st.time.last_dr = dr;

        int correction_axes = 0;
        bool anomaly = false;
        for (int a = 0; a < 3; ++a) {
            if (std::abs(dr_raw.v[a] - dr.v[a]) > 1e-10) ++correction_axes;
            if (fr.pbc.L[a] > 0.0 && std::abs(dr.v[a]) > 0.5 * fr.pbc.L[a] + 1e-10) anomaly = true;
        }

        const double step = norm3(dr);
        st.time.last_step_norm = step;
        st.time.pbc_correction_axes = correction_axes;
        st.time.step_anomaly = anomaly;
        if (step > st.time.max_jump) st.time.max_jump = step;
    }
};

struct Pipeline {
    int nthreads = 1;
    std::vector<std::unique_ptr<Property>> props;

    explicit Pipeline(int nthreads_) : nthreads(std::max(1, nthreads_)) {}

    void add(std::unique_ptr<Property> p) { props.emplace_back(std::move(p)); }

    void run_frame_stage(Stage stage, const Topology& topo, Frame& fr, std::vector<MolState>& ms) {
        for (auto& p : props) {
            if (p->stage() == stage && p->scope() == Scope::Frame) {
                p->run_frame(topo, fr, ms);
            }
        }
    }

    void run_mol_stage(Stage stage, MolType type, const Topology& topo, Frame& fr,
                       std::vector<MolState>& ms) {
        const auto& mol_ids = topo.mol_ids_by_type[(int)type];
        parallel_for_molids_strided(nthreads, mol_ids, [&](int mol_id, int tid) {
            for (auto& p : props) {
                if (p->stage() == stage && p->scope() == Scope::Mol) {
                    p->run_mol(topo, fr, mol_id, ms, tid);
                }
            }
        });
    }

    void run_mol_stage_all(Stage stage, const Topology& topo, Frame& fr, std::vector<MolState>& ms) {
        run_mol_stage(stage, MolType::Water, topo, fr, ms);
        run_mol_stage(stage, MolType::Cation, topo, fr, ms);
        run_mol_stage(stage, MolType::Anion, topo, fr, ms);
        run_mol_stage(stage, MolType::Other, topo, fr, ms);
    }

    void run_atom_stage(Stage stage, const Topology& topo, Frame& fr, std::vector<MolState>& ms) {
        const int natoms = (int)fr.atoms.size();
        parallel_strided(nthreads, natoms, [&](int atom_id, int tid) {
            for (auto& p : props) {
                if (p->stage() == stage && p->scope() == Scope::Atom) {
                    p->run_atom(topo, fr, atom_id, ms, tid);
                }
            }
        });
    }

    void process_frame(const Topology& topo, Frame& fr, std::vector<MolState>& ms) {
        assert(ms.size() == topo.mols.size());

        run_mol_stage_all(Stage::BeginFrame, topo, fr, ms);
        run_atom_stage(Stage::BeginFrame, topo, fr, ms);
        run_frame_stage(Stage::IntrinsicAndGrid, topo, fr, ms);
        run_mol_stage_all(Stage::Neighbor, topo, fr, ms);
        run_atom_stage(Stage::Neighbor, topo, fr, ms);
        run_mol_stage_all(Stage::Time, topo, fr, ms);
        run_atom_stage(Stage::Time, topo, fr, ms);
        run_frame_stage(Stage::Reduce, topo, fr, ms);
        run_atom_stage(Stage::Reduce, topo, fr, ms);
        run_frame_stage(Stage::EndFrame, topo, fr, ms);
        run_atom_stage(Stage::EndFrame, topo, fr, ms);
    }
};

inline Pipeline make_default_pipeline(int nthreads, double grid_cell_size_nm = 0.5) {
    Pipeline pipe(nthreads);

    pipe.add(std::make_unique<BeginFrameClear>());
    pipe.add(std::make_unique<BuildIntrinsicsAndGrid>(nthreads, grid_cell_size_nm));
    pipe.add(std::make_unique<TimeDisplacement>());

    return pipe;
}

struct WaterHBondsStub : Property {
    double r_grid = 0.6;

    const char* name() const override { return "WaterHBondsStub"; }
    Scope scope() const override { return Scope::Mol; }
    Stage stage() const override { return Stage::Neighbor; }

    void run_mol(const Topology& topo, Frame& fr, int mol_id, std::vector<MolState>& ms,
                 int) override {
        const MolSpan& m = topo.mols[(size_t)mol_id];
        if (m.type != MolType::Water) return;

        MolState& st = ms[(size_t)mol_id];
        if (!(st.cache.flags & MolCache::HAS_KEY)) return;

        Vec3d key = st.cache.key_wrapped;
        fr.pbc.wrap_pos3(key);

        int nx_cell = std::max(1, (int)std::ceil(r_grid / fr.grid.cx));
        int ny_cell = std::max(1, (int)std::ceil(r_grid / fr.grid.cy));
        int nz_cell = std::max(1, (int)std::ceil(r_grid / fr.grid.cz));

        fr.grid.for_candidates_box(key, nx_cell, ny_cell, nz_cell, [&](int cand_id) {
            if (cand_id == mol_id) return;
        });
    }
};

}  // namespace simio
