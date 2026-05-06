#include "nacl_association_x.hpp"

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

int xbin_from_x(double xw, double dx, int nx) {
    const int ix = static_cast<int>(std::floor(xw / dx));
    if (ix < 0 || ix >= nx) return -1;
    return ix;
}

void add_point(std::vector<double>& bins, double xw, double dx, int nx, double value) {
    const int ix = xbin_from_x(xw, dx, nx);
    if (ix >= 0) bins[static_cast<size_t>(ix)] += value;
}

void add_segment_fractional(std::vector<double>& bins,
                            const Pbc3D& pbc,
                            double xa,
                            double xb_wrapped,
                            double dx,
                            int nx,
                            double value) {
    const double Lx = pbc.L[0];
    const double xb = xa + pbc.min_image(Vec3d{{xb_wrapped - xa, 0.0, 0.0}}).v[0];
    double lo = std::min(xa, xb);
    double hi = std::max(xa, xb);
    const double len = hi - lo;
    if (len <= 1e-12) {
        add_point(bins, pbc.wrap_pos(0, xa), dx, nx, value);
        return;
    }

    int ib0 = static_cast<int>(std::floor(lo / dx));
    int ib1 = static_cast<int>(std::floor((hi - 1e-12) / dx));
    for (int ib = ib0; ib <= ib1; ++ib) {
        const double blo = static_cast<double>(ib) * dx;
        const double bhi = blo + dx;
        const double overlap = std::max(0.0, std::min(hi, bhi) - std::max(lo, blo));
        if (overlap <= 0.0) continue;
        const double mid = 0.5 * (std::max(lo, blo) + std::min(hi, bhi));
        add_point(bins, pbc.wrap_pos(0, mid), dx, nx, value * overlap / len);
    }
}

void set_segment_max(std::vector<double>& bins,
                     const Pbc3D& pbc,
                     double xa,
                     double xb_wrapped,
                     double dx,
                     int nx,
                     double value) {
    const double xb = xa + pbc.min_image(Vec3d{{xb_wrapped - xa, 0.0, 0.0}}).v[0];
    const double lo = std::min(xa, xb);
    const double hi = std::max(xa, xb);
    if (hi - lo <= 1e-12) {
        const int ix = xbin_from_x(pbc.wrap_pos(0, xa), dx, nx);
        if (ix >= 0) bins[static_cast<size_t>(ix)] = std::max(bins[static_cast<size_t>(ix)], value);
        return;
    }

    const int ib0 = static_cast<int>(std::floor(lo / dx));
    const int ib1 = static_cast<int>(std::floor((hi - 1e-12) / dx));
    for (int ib = ib0; ib <= ib1; ++ib) {
        const double blo = static_cast<double>(ib) * dx;
        const double bhi = blo + dx;
        const double overlap = std::max(0.0, std::min(hi, bhi) - std::max(lo, blo));
        if (overlap <= 0.0) continue;
        const double mid = 0.5 * (std::max(lo, blo) + std::min(hi, bhi));
        const int ix = xbin_from_x(pbc.wrap_pos(0, mid), dx, nx);
        if (ix >= 0) bins[static_cast<size_t>(ix)] = std::max(bins[static_cast<size_t>(ix)], value);
    }
}

void set_cluster_span_max(std::vector<double>& bins,
                          const Pbc3D& pbc,
                          const std::vector<Vec3d>& points,
                          double dx,
                          int nx,
                          double value) {
    if (points.empty()) return;
    if (points.size() == 1) {
        const int ix = xbin_from_x(points.front().v[0], dx, nx);
        if (ix >= 0) bins[static_cast<size_t>(ix)] = std::max(bins[static_cast<size_t>(ix)], value);
        return;
    }

    const double x0 = points.front().v[0];
    double lo = x0;
    double hi = x0;
    for (const Vec3d& p : points) {
        const double xu = x0 + pbc.min_image(p - points.front()).v[0];
        lo = std::min(lo, xu);
        hi = std::max(hi, xu);
    }
    set_segment_max(bins, pbc, lo, pbc.wrap_pos(0, hi), dx, nx, value);
}

void add_three_point_span_fractional(std::vector<double>& bins,
                                     const Pbc3D& pbc,
                                     const Vec3d& na,
                                     const Vec3d& ow,
                                     const Vec3d& cl,
                                     double dx,
                                     int nx,
                                     double value) {
    const double x0 = na.v[0];
    const double xw = x0 + pbc.min_image(ow - na).v[0];
    const double xc = x0 + pbc.min_image(cl - na).v[0];
    const double lo = std::min({x0, xw, xc});
    const double hi = std::max({x0, xw, xc});
    add_segment_fractional(bins, pbc, lo, pbc.wrap_pos(0, hi), dx, nx, value);
}

bool in_z_slab(const Pbc3D& pbc, const Vec3d& r, double zminw, double zmaxw) {
    return in_range_pbc(pbc.wrap_pos(2, r.v[2]), zminw, zmaxw);
}

}  // namespace

NaClAssociationXAnalyzer::NaClAssociationXAnalyzer(const NaClAssociationXConfig& cfg) : cfg_(cfg) {
    if (cfg_.nx <= 0) throw std::runtime_error("NaClAssociationXAnalyzer: nx must be > 0");
    if (cfg_.r_cip <= 0.0 || cfg_.r_ssip <= cfg_.r_cip || cfg_.r_naow <= 0.0 || cfg_.r_clow <= 0.0) {
        throw std::runtime_error("NaClAssociationXAnalyzer: invalid cutoffs");
    }
    cip_.init(cfg_.nx);
    ssip_.init(cfg_.nx);
    assoc_.init(cfg_.nx);
    bridge_water_.init(cfg_.nx);
    bridged_pair_.init(cfg_.nx);
    f_cip_.init(cfg_.nx);
    f_ssip_.init(cfg_.nx);
    f_bridge_.init(cfg_.nx);
    cn_shared_.init(cfg_.nx);
    largest_cluster_size_.init(cfg_.nx);
    cn_naow_.init(cfg_.nx);
    cn_clow_.init(cfg_.nx);
}

NaClAssociationXAnalyzer::NaClAssociationXAnalyzer(const NaClAssociationXConfig& cfg,
                                                   simio::runtime::CacheStore& cache)
    : NaClAssociationXAnalyzer(cfg) {
    cache_ = &cache;
    const auto grid = simio::analysis::intrinsics::get_or_build_x_grid(*cache_, 0.0, 1.0, cfg_.nx);
    dx_ = grid.dx;
    x_centers_rel_ = grid.centers_rel;
    has_cached_rel_grid_ = true;
}

void NaClAssociationXAnalyzer::process_frame(const Topology& topo,
                                             const Frame& fr,
                                             const std::vector<MolState>& ms) {
    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Lz <= 0.0) throw std::runtime_error("NaClAssociationXAnalyzer: invalid box lengths");

    if (!has_ref_box_) {
        has_ref_box_ = true;
        if (has_cached_rel_grid_) {
            dx_ *= Lx;
            for (double& c : x_centers_rel_) c *= Lx;
        } else {
            dx_ = Lx / static_cast<double>(cfg_.nx);
        }
    }
    if (dx_ <= 0.0) throw std::runtime_error("NaClAssociationXAnalyzer: invalid dx");

    const double zminw = fr.pbc.wrap_pos(2, cfg_.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg_.zmax);
    const double r2_cip = cfg_.r_cip * cfg_.r_cip;
    const double r2_ssip = cfg_.r_ssip * cfg_.r_ssip;
    const double r2_naow = cfg_.r_naow * cfg_.r_naow;
    const double r2_clow = cfg_.r_clow * cfg_.r_clow;
    const double rmax_hyd = std::max(cfg_.r_naow, cfg_.r_clow);

    const int nx_ssip = std::max(1, static_cast<int>(std::ceil((cfg_.r_ssip + 0.05) / fr.grid.cx)));
    const int ny_ssip = std::max(1, static_cast<int>(std::ceil((cfg_.r_ssip + 0.05) / fr.grid.cy)));
    const int nz_ssip = std::max(1, static_cast<int>(std::ceil((cfg_.r_ssip + 0.05) / fr.grid.cz)));
    const int nx_hyd = std::max(1, static_cast<int>(std::ceil((rmax_hyd + 0.05) / fr.grid.cx)));
    const int ny_hyd = std::max(1, static_cast<int>(std::ceil((rmax_hyd + 0.05) / fr.grid.cy)));
    const int nz_hyd = std::max(1, static_cast<int>(std::ceil((rmax_hyd + 0.05) / fr.grid.cz)));

    std::vector<Vec3d> pos(topo.mols.size(), Vec3d{});
    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        const MolCache& c = ms[mid].cache;
        if (c.flags & MolCache::HAS_REF) {
            pos[mid] = c.ref_wrapped;
        } else if (c.flags & MolCache::HAS_KEY) {
            pos[mid] = c.key_wrapped;
        }
        fr.pbc.wrap_pos3(pos[mid]);
    }

    std::vector<double> frame_cip(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_ssip(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_assoc(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_bridge_water(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_bridged_pair(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_bridged_pair_for_fraction(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_shared_waters_for_ssip(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> cn_na_sum(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> cn_cl_sum(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<int64_t> cn_na_n(static_cast<size_t>(cfg_.nx), 0);
    std::vector<int64_t> cn_cl_n(static_cast<size_t>(cfg_.nx), 0);

    std::vector<int> ion_nodes;
    std::vector<int> parent;
    std::vector<int> rank;
    std::vector<int> node_for_mol(topo.mols.size(), -1);

    auto ensure_node = [&](int mid) {
        int& node = node_for_mol[static_cast<size_t>(mid)];
        if (node >= 0) return node;
        node = static_cast<int>(ion_nodes.size());
        ion_nodes.push_back(mid);
        parent.push_back(node);
        rank.push_back(0);
        return node;
    };

    for (int na_mid : topo.mol_ids_by_type[static_cast<int>(MolType::Cation)]) {
        const Vec3d na = pos[static_cast<size_t>(na_mid)];
        if (!in_z_slab(fr.pbc, na, zminw, zmaxw)) continue;

        int cn_na = 0;
        fr.grid.for_candidates_box(na, nx_hyd, ny_hyd, nz_hyd, [&](int cand_id) {
            if (topo.mols[static_cast<size_t>(cand_id)].type != MolType::Water) return;
            const Vec3d ow = pos[static_cast<size_t>(cand_id)];
            const Vec3d dr = fr.pbc.min_image(ow - na);
            if (dot3(dr, dr) <= r2_naow) ++cn_na;
        });
        const int na_bin = xbin_from_x(na.v[0], dx_, cfg_.nx);
        if (na_bin >= 0) {
            cn_na_sum[static_cast<size_t>(na_bin)] += static_cast<double>(cn_na);
            cn_na_n[static_cast<size_t>(na_bin)] += 1;
        }

        fr.grid.for_candidates_box(na, nx_ssip, ny_ssip, nz_ssip, [&](int cl_mid) {
            if (topo.mols[static_cast<size_t>(cl_mid)].type != MolType::Anion) return;
            const Vec3d cl = pos[static_cast<size_t>(cl_mid)];

            const Vec3d dr_nacl = fr.pbc.min_image(cl - na);
            const double d2_nacl = dot3(dr_nacl, dr_nacl);
            if (d2_nacl > r2_ssip) return;

            ensure_node(na_mid);
            ensure_node(cl_mid);
            unite_roots(parent, rank, node_for_mol[static_cast<size_t>(na_mid)],
                        node_for_mol[static_cast<size_t>(cl_mid)]);

            if (d2_nacl <= r2_cip) {
                add_segment_fractional(frame_cip, fr.pbc, na.v[0], cl.v[0], dx_, cfg_.nx, 1.0);
            } else {
                add_segment_fractional(frame_ssip, fr.pbc, na.v[0], cl.v[0], dx_, cfg_.nx, 1.0);

                int shared_count = 0;
                Vec3d bridge_ow{};
                fr.grid.for_candidates_box(na, nx_hyd, ny_hyd, nz_hyd, [&](int water_mid) {
                    if (topo.mols[static_cast<size_t>(water_mid)].type != MolType::Water) return;
                    const Vec3d ow = pos[static_cast<size_t>(water_mid)];
                    const Vec3d d_naow = fr.pbc.min_image(ow - na);
                    if (dot3(d_naow, d_naow) > r2_naow) return;
                    const Vec3d d_clow = fr.pbc.min_image(ow - cl);
                    if (dot3(d_clow, d_clow) > r2_clow) return;
                    if (shared_count == 0) bridge_ow = ow;
                    ++shared_count;
                });
                if (shared_count > 0) {
                    add_three_point_span_fractional(frame_bridged_pair, fr.pbc, na, bridge_ow, cl,
                                                    dx_, cfg_.nx, 1.0);
                    add_segment_fractional(frame_bridged_pair_for_fraction, fr.pbc, na.v[0], cl.v[0],
                                           dx_, cfg_.nx, 1.0);
                }
                add_segment_fractional(frame_shared_waters_for_ssip, fr.pbc, na.v[0], cl.v[0],
                                       dx_, cfg_.nx, static_cast<double>(shared_count));
            }
            add_segment_fractional(frame_assoc, fr.pbc, na.v[0], cl.v[0], dx_, cfg_.nx, 1.0);
        });
    }

    for (int cl_mid : topo.mol_ids_by_type[static_cast<int>(MolType::Anion)]) {
        const Vec3d cl = pos[static_cast<size_t>(cl_mid)];
        if (!in_z_slab(fr.pbc, cl, zminw, zmaxw)) continue;
        int cn_cl = 0;
        fr.grid.for_candidates_box(cl, nx_hyd, ny_hyd, nz_hyd, [&](int cand_id) {
            if (topo.mols[static_cast<size_t>(cand_id)].type != MolType::Water) return;
            const Vec3d ow = pos[static_cast<size_t>(cand_id)];
            const Vec3d dr = fr.pbc.min_image(ow - cl);
            if (dot3(dr, dr) <= r2_clow) ++cn_cl;
        });
        const int cl_bin = xbin_from_x(cl.v[0], dx_, cfg_.nx);
        if (cl_bin >= 0) {
            cn_cl_sum[static_cast<size_t>(cl_bin)] += static_cast<double>(cn_cl);
            cn_cl_n[static_cast<size_t>(cl_bin)] += 1;
        }
    }

    for (int water_mid : topo.mol_ids_by_type[static_cast<int>(MolType::Water)]) {
        const Vec3d ow = pos[static_cast<size_t>(water_mid)];
        if (!in_z_slab(fr.pbc, ow, zminw, zmaxw)) continue;
        bool near_na = false;
        bool near_cl = false;
        fr.grid.for_candidates_box(ow, nx_hyd, ny_hyd, nz_hyd, [&](int cand_id) {
            const MolType t = topo.mols[static_cast<size_t>(cand_id)].type;
            if (t != MolType::Cation && t != MolType::Anion) return;
            const Vec3d ion = pos[static_cast<size_t>(cand_id)];
            const Vec3d dr = fr.pbc.min_image(ion - ow);
            const double d2 = dot3(dr, dr);
            if (t == MolType::Cation && d2 <= r2_naow) near_na = true;
            if (t == MolType::Anion && d2 <= r2_clow) near_cl = true;
        });
        if (near_na && near_cl) add_point(frame_bridge_water, ow.v[0], dx_, cfg_.nx, 1.0);
    }

    std::vector<std::vector<int>> components(parent.size());
    for (int i = 0; i < static_cast<int>(parent.size()); ++i) {
        const int r = find_root(parent, i);
        components[static_cast<size_t>(r)].push_back(i);
    }
    std::vector<double> frame_largest_cluster_size(static_cast<size_t>(cfg_.nx), 0.0);
    for (const std::vector<int>& comp : components) {
        if (comp.size() < 2) continue;
        std::vector<Vec3d> points;
        points.reserve(comp.size());
        for (int node : comp) points.push_back(pos[static_cast<size_t>(ion_nodes[static_cast<size_t>(node)])]);
        set_cluster_span_max(frame_largest_cluster_size, fr.pbc, points, dx_, cfg_.nx,
                             static_cast<double>(comp.size()));
    }

    for (int i = 0; i < cfg_.nx; ++i) {
        const size_t b = static_cast<size_t>(i);
        cip_.add(i, frame_cip[b]);
        ssip_.add(i, frame_ssip[b]);
        assoc_.add(i, frame_assoc[b]);
        bridge_water_.add(i, frame_bridge_water[b]);
        bridged_pair_.add(i, frame_bridged_pair[b]);
        f_cip_.add(i, frame_assoc[b] > 0.0 ? frame_cip[b] / frame_assoc[b] : 0.0);
        f_ssip_.add(i, frame_assoc[b] > 0.0 ? frame_ssip[b] / frame_assoc[b] : 0.0);
        f_bridge_.add(i, frame_ssip[b] > 0.0 ? frame_bridged_pair_for_fraction[b] / frame_ssip[b] : 0.0);
        cn_shared_.add(i, frame_ssip[b] > 0.0 ? frame_shared_waters_for_ssip[b] / frame_ssip[b] : 0.0);
        largest_cluster_size_.add(i, frame_largest_cluster_size[b]);
        if (cn_na_n[b] > 0) cn_naow_.add(i, cn_na_sum[b] / static_cast<double>(cn_na_n[b]));
        if (cn_cl_n[b] > 0) cn_clow_.add(i, cn_cl_sum[b] / static_cast<double>(cn_cl_n[b]));
    }

    nframes_ += 1;
}

void NaClAssociationXAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("NaClAssociationXAnalyzer: no frames processed");
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("NaClAssociationXAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "x_center_nm,N_CIP_mean,N_CIP_sem,N_SSIP_mean,N_SSIP_sem,N_ASSOC_mean,N_ASSOC_sem,"
           "CN_NaOW_mean,CN_NaOW_sem,CN_ClOW_mean,CN_ClOW_sem,"
           "N_bridge_water_mean,N_bridge_water_sem,N_bridged_pair_mean,N_bridged_pair_sem,"
           "f_CIP_mean,f_CIP_sem,f_SSIP_mean,f_SSIP_sem,f_bridge_mean,f_bridge_sem,"
           "CN_shared_mean,CN_shared_sem,"
           "largest_cluster_size_mean,largest_cluster_size_sem\n";
    for (int i = 0; i < cfg_.nx; ++i) {
        const double x_center = (static_cast<double>(i) + 0.5) * dx_;
        ofs << x_center << "," << cip_.mean(i, nframes_) << "," << cip_.sem(i, nframes_) << ","
            << ssip_.mean(i, nframes_) << "," << ssip_.sem(i, nframes_) << ","
            << assoc_.mean(i, nframes_) << "," << assoc_.sem(i, nframes_) << ","
            << cn_naow_.mean(i) << "," << cn_naow_.sem(i) << ","
            << cn_clow_.mean(i) << "," << cn_clow_.sem(i) << ","
            << bridge_water_.mean(i, nframes_) << "," << bridge_water_.sem(i, nframes_) << ","
            << bridged_pair_.mean(i, nframes_) << "," << bridged_pair_.sem(i, nframes_) << ","
            << f_cip_.mean(i, nframes_) << "," << f_cip_.sem(i, nframes_) << ","
            << f_ssip_.mean(i, nframes_) << "," << f_ssip_.sem(i, nframes_) << ","
            << f_bridge_.mean(i, nframes_) << "," << f_bridge_.sem(i, nframes_) << ","
            << cn_shared_.mean(i, nframes_) << "," << cn_shared_.sem(i, nframes_) << ","
            << largest_cluster_size_.mean(i, nframes_) << "," << largest_cluster_size_.sem(i, nframes_)
            << "\n";
    }
}

}  // namespace simio::analysis
