#include "nacl_cluster_x.hpp"

#include "simio/analysis/intrinsics/x_grid_cache.hpp"
#include "simio/runtime/cache.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace simio::analysis {
namespace {

int find_root(std::vector<int>& parent, int i) {
    int r = i;
    while (parent[static_cast<size_t>(r)] != r) r = parent[static_cast<size_t>(r)];
    while (parent[static_cast<size_t>(i)] != i) {
        const int p = parent[static_cast<size_t>(i)];
        parent[static_cast<size_t>(i)] = r;
        i = p;
    }
    return r;
}

void unite_roots(std::vector<int>& parent, std::vector<int>& rank, int a, int b) {
    int ra = find_root(parent, a);
    int rb = find_root(parent, b);
    if (ra == rb) return;
    if (rank[static_cast<size_t>(ra)] < rank[static_cast<size_t>(rb)]) std::swap(ra, rb);
    parent[static_cast<size_t>(rb)] = ra;
    if (rank[static_cast<size_t>(ra)] == rank[static_cast<size_t>(rb)]) {
        rank[static_cast<size_t>(ra)] += 1;
    }
}

}  // namespace

NaClClusterXAnalyzer::NaClClusterXAnalyzer(const NaClClusterXConfig& cfg) : cfg_(cfg) {
    if (cfg_.nx <= 0) throw std::runtime_error("NaClClusterXAnalyzer: nx must be > 0");
    if (cfg_.r_nacl <= 0.0) throw std::runtime_error("NaClClusterXAnalyzer: r_nacl must be > 0");

    cluster_count_.init(cfg_.nx);
    ion_count_in_clusters_.init(cfg_.nx);
    na_count_in_clusters_.init(cfg_.nx);
    cl_count_in_clusters_.init(cfg_.nx);
    cluster_size_.init(cfg_.nx);
}

NaClClusterXAnalyzer::NaClClusterXAnalyzer(const NaClClusterXConfig& cfg,
                                           simio::runtime::CacheStore& cache)
    : NaClClusterXAnalyzer(cfg) {
    cache_ = &cache;
    const auto grid = simio::analysis::intrinsics::get_or_build_x_grid(*cache_, 0.0, 1.0, cfg_.nx);
    dx_ = grid.dx;
    x_centers_rel_ = grid.centers_rel;
    has_cached_rel_grid_ = true;
}

void NaClClusterXAnalyzer::process_frame(const Topology& topo,
                                         const Frame& fr,
                                         const std::vector<MolState>& ms) {
    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Lz <= 0.0) {
        throw std::runtime_error("NaClClusterXAnalyzer: invalid box lengths");
    }

    if (!has_ref_box_) {
        has_ref_box_ = true;
        Lx_ref_ = Lx;
        if (has_cached_rel_grid_) {
            dx_ *= Lx_ref_;
            for (double& c : x_centers_rel_) c *= Lx_ref_;
        } else {
            dx_ = Lx_ref_ / static_cast<double>(cfg_.nx);
            x_centers_rel_.resize(static_cast<size_t>(cfg_.nx));
            for (int i = 0; i < cfg_.nx; ++i) {
                x_centers_rel_[static_cast<size_t>(i)] = (static_cast<double>(i) + 0.5) * dx_;
            }
        }
    }
    if (dx_ <= 0.0) throw std::runtime_error("NaClClusterXAnalyzer: invalid dx");

    const double zminw = fr.pbc.wrap_pos(2, cfg_.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg_.zmax);
    const double r2_nacl = cfg_.r_nacl * cfg_.r_nacl;

    std::vector<IonNode> ions;
    ions.reserve(topo.mols.size());
    std::vector<int> na_nodes;
    std::vector<int> cl_nodes;

    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        Species species{};
        if (topo.mols[mid].type == MolType::Cation) {
            species = Species::Na;
        } else if (topo.mols[mid].type == MolType::Anion) {
            species = Species::Cl;
        } else {
            continue;
        }

        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) continue;

        Vec3d key = c.key_wrapped;
        fr.pbc.wrap_pos3(key);
        if (!in_range_pbc(key.v[2], zminw, zmaxw)) continue;

        const int idx = static_cast<int>(ions.size());
        ions.push_back(IonNode{static_cast<int>(mid), species, key});
        if (species == Species::Na) {
            na_nodes.push_back(idx);
        } else {
            cl_nodes.push_back(idx);
        }
    }

    std::vector<int> parent(ions.size());
    std::vector<int> rank(ions.size(), 0);
    std::iota(parent.begin(), parent.end(), 0);
    std::vector<uint8_t> linked(ions.size(), 0);

    for (int ni : na_nodes) {
        const Vec3d& rna = ions[static_cast<size_t>(ni)].pos;
        for (int ci : cl_nodes) {
            const Vec3d dr = fr.pbc.min_image(ions[static_cast<size_t>(ci)].pos - rna);
            if (dot3(dr, dr) > r2_nacl) continue;
            unite_roots(parent, rank, ni, ci);
            linked[static_cast<size_t>(ni)] = 1;
            linked[static_cast<size_t>(ci)] = 1;
        }
    }

    std::vector<double> frame_cluster_count(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_ion_count(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_na_count(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_cl_count(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_size_sum(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<int64_t> frame_size_n(static_cast<size_t>(cfg_.nx), 0);

    std::vector<std::vector<int>> components(ions.size());
    for (int i = 0; i < static_cast<int>(ions.size()); ++i) {
        if (linked[static_cast<size_t>(i)] == 0) continue;
        components[static_cast<size_t>(find_root(parent, i))].push_back(i);
    }

    for (const std::vector<int>& comp : components) {
        if (comp.size() < 2) continue;

        int n_na = 0;
        int n_cl = 0;
        Vec3d ref = ions[static_cast<size_t>(comp.front())].pos;
        Vec3d sum = ref;
        for (size_t k = 0; k < comp.size(); ++k) {
            const IonNode& ion = ions[static_cast<size_t>(comp[k])];
            if (ion.species == Species::Na) {
                ++n_na;
            } else {
                ++n_cl;
            }
            if (k == 0) continue;
            sum = sum + pbc_unwrap_to_ref(fr.pbc, ref, ion.pos);
        }
        if (n_na <= 0 || n_cl <= 0) continue;

        Vec3d center = sum * (1.0 / static_cast<double>(comp.size()));
        fr.pbc.wrap_pos3(center);
        const int ix = static_cast<int>(std::floor(center.v[0] / dx_));
        if (ix < 0 || ix >= cfg_.nx) continue;

        const size_t b = static_cast<size_t>(ix);
        frame_cluster_count[b] += 1.0;
        frame_ion_count[b] += static_cast<double>(comp.size());
        frame_na_count[b] += static_cast<double>(n_na);
        frame_cl_count[b] += static_cast<double>(n_cl);
        frame_size_sum[b] += static_cast<double>(comp.size());
        frame_size_n[b] += 1;
    }

    for (int i = 0; i < cfg_.nx; ++i) {
        const size_t b = static_cast<size_t>(i);
        cluster_count_.add(i, frame_cluster_count[b]);
        ion_count_in_clusters_.add(i, frame_ion_count[b]);
        na_count_in_clusters_.add(i, frame_na_count[b]);
        cl_count_in_clusters_.add(i, frame_cl_count[b]);
        if (frame_size_n[b] > 0) {
            cluster_size_.add(i, frame_size_sum[b] / static_cast<double>(frame_size_n[b]));
        }
    }

    nframes_ += 1;
}

void NaClClusterXAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("NaClClusterXAnalyzer: no frames processed");
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("NaClClusterXAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "x_center_nm,nacl_cluster_count_mean,nacl_cluster_count_sem,"
           "ions_in_nacl_clusters_mean,ions_in_nacl_clusters_sem,"
           "na_in_nacl_clusters_mean,na_in_nacl_clusters_sem,"
           "cl_in_nacl_clusters_mean,cl_in_nacl_clusters_sem,"
           "nacl_cluster_size_mean,nacl_cluster_size_sem\n";

    for (int i = 0; i < cfg_.nx; ++i) {
        const double x_center = (static_cast<double>(i) + 0.5) * dx_;
        ofs << x_center << "," << cluster_count_.mean(i, nframes_) << ","
            << cluster_count_.sem(i, nframes_) << ","
            << ion_count_in_clusters_.mean(i, nframes_) << ","
            << ion_count_in_clusters_.sem(i, nframes_) << ","
            << na_count_in_clusters_.mean(i, nframes_) << ","
            << na_count_in_clusters_.sem(i, nframes_) << ","
            << cl_count_in_clusters_.mean(i, nframes_) << ","
            << cl_count_in_clusters_.sem(i, nframes_) << ","
            << cluster_size_.mean(i) << "," << cluster_size_.sem(i) << "\n";
    }
}

}  // namespace simio::analysis
