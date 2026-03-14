#include "water_atom_density_x.hpp"

#include "simio/analysis/intrinsics/x_grid_cache.hpp"
#include "simio/runtime/cache.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace simio::analysis {

WaterAtomDensityXAnalyzer::WaterAtomDensityXAnalyzer(const WaterAtomDensityXConfig& cfg)
    : cfg_(cfg) {
    if (cfg_.nx <= 0) throw std::runtime_error("WaterAtomDensityXAnalyzer: nx must be > 0");
    for (int s = 0; s < SiteN; ++s) {
        rho_[static_cast<size_t>(s)].init(cfg_.nx);
        cnt_[static_cast<size_t>(s)].init(cfg_.nx);
    }
}

WaterAtomDensityXAnalyzer::WaterAtomDensityXAnalyzer(const WaterAtomDensityXConfig& cfg,
                                                     simio::runtime::CacheStore& cache)
    : WaterAtomDensityXAnalyzer(cfg) {
    cache_ = &cache;
    const auto grid = simio::analysis::intrinsics::get_or_build_x_grid(*cache_, 0.0, 1.0, cfg_.nx);
    dx_ = grid.dx;
    x_centers_rel_ = grid.centers_rel;
    has_cached_rel_grid_ = true;
}

void WaterAtomDensityXAnalyzer::process_frame(const Topology& topo,
                                              const Frame& fr,
                                              const std::vector<MolState>& ms) {
    const double Lx = fr.pbc.L[0];
    const double Ly = fr.pbc.L[1];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Ly <= 0.0 || Lz <= 0.0) {
        throw std::runtime_error("WaterAtomDensityXAnalyzer: invalid box lengths");
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
    const double dx = dx_;

    const double zminw = fr.pbc.wrap_pos(2, cfg_.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg_.zmax);
    const double zlen = interval_length_pbc(zminw, zmaxw, Lz);
    if (dx <= 0.0 || zlen <= 0.0) {
        throw std::runtime_error("WaterAtomDensityXAnalyzer: invalid bin/slice dimensions");
    }

    std::array<std::vector<int64_t>, SiteN> frame_counts{};
    for (int s = 0; s < SiteN; ++s) {
        frame_counts[static_cast<size_t>(s)].assign(static_cast<size_t>(cfg_.nx), 0);
    }

    auto count_site = [&](Site site, Vec3d pos_wrapped) {
        fr.pbc.wrap_pos3(pos_wrapped);
        if (!in_range_pbc(pos_wrapped.v[2], zminw, zmaxw)) return;
        const int ix = static_cast<int>(std::floor(pos_wrapped.v[0] / dx));
        if (ix < 0 || ix >= cfg_.nx) return;
        frame_counts[static_cast<size_t>(site)][static_cast<size_t>(ix)] += 1;
    };

    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        const MolSpan& m = topo.mols[mid];
        if (m.type != MolType::Water || m.natoms < 3) continue;

        const MolCache& c = ms[mid].cache;

        Vec3d rO{};
        Vec3d h1w{};
        Vec3d h2w{};
        if ((c.flags & MolCache::HAS_REF) && (c.flags & MolCache::HAS_SITES)) {
            rO = c.ref_wrapped;
            h1w = c.sites_u[1];
            h2w = c.sites_u[2];
        } else {
            rO = fr.atoms.pos(m.first + 0);
            const Vec3d rH1 = fr.atoms.pos(m.first + 1);
            const Vec3d rH2 = fr.atoms.pos(m.first + 2);
            fr.pbc.wrap_pos3(rO);
            h1w = pbc_unwrap_to_ref(fr.pbc, rO, rH1);
            h2w = pbc_unwrap_to_ref(fr.pbc, rO, rH2);
        }

        count_site(Site::Ow, rO);
        count_site(Site::Hw, h1w);
        count_site(Site::Hw, h2w);
    }

    const double bin_volume = dx * Ly * zlen;
    for (int i = 0; i < cfg_.nx; ++i) {
        for (int s = 0; s < SiteN; ++s) {
            const int64_t c = frame_counts[static_cast<size_t>(s)][static_cast<size_t>(i)];
            const double rho = static_cast<double>(c) / bin_volume;
            rho_[static_cast<size_t>(s)].add(i, rho);
            cnt_[static_cast<size_t>(s)].add(i, static_cast<double>(c));
        }
    }

    nframes_ += 1;
}

void WaterAtomDensityXAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("WaterAtomDensityXAnalyzer: no frames processed");
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("WaterAtomDensityXAnalyzer: failed to open output CSV: " + path);
    }

    ofs << std::setprecision(12);
    ofs << "x_center_nm,"
           "rho_ow_mean,rho_ow_sem,rho_hw_mean,rho_hw_sem,"
           "count_ow_mean,count_ow_sem,count_hw_mean,count_hw_sem\n";

    for (int i = 0; i < cfg_.nx; ++i) {
        const double x_center = (static_cast<double>(i) + 0.5) * dx_;
        ofs << x_center << ","
            << rho_[Site::Ow].mean(i, nframes_) << "," << rho_[Site::Ow].sem(i, nframes_) << ","
            << rho_[Site::Hw].mean(i, nframes_) << "," << rho_[Site::Hw].sem(i, nframes_) << ","
            << cnt_[Site::Ow].mean(i, nframes_) << "," << cnt_[Site::Ow].sem(i, nframes_) << ","
            << cnt_[Site::Hw].mean(i, nframes_) << "," << cnt_[Site::Hw].sem(i, nframes_) << "\n";
    }
}

}  // namespace simio::analysis
