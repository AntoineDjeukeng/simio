#include "density_x.hpp"

#include "simio/analysis/intrinsics/x_grid_cache.hpp"
#include "simio/runtime/cache.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace simio::analysis {

DensityXAnalyzer::DensityXAnalyzer(const DensityXConfig& cfg) : cfg_(cfg) {
    if (cfg_.nx <= 0) throw std::runtime_error("DensityXAnalyzer: nx must be > 0");
    for (int s = 0; s < 3; ++s) {
        rho_[static_cast<size_t>(s)].init(cfg_.nx);
        cnt_[static_cast<size_t>(s)].init(cfg_.nx);
    }
}

DensityXAnalyzer::DensityXAnalyzer(const DensityXConfig& cfg, simio::runtime::CacheStore& cache)
    : DensityXAnalyzer(cfg) {
    cache_ = &cache;
    // Build/get a relative x-grid once per run (unit interval [0,1]).
    const auto grid = simio::analysis::intrinsics::get_or_build_x_grid(*cache_, 0.0, 1.0, cfg_.nx);
    dx_ = grid.dx;
    x_centers_rel_ = grid.centers_rel;
    has_cached_rel_grid_ = true;
}

void DensityXAnalyzer::process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms) {
    const double Lx = fr.pbc.L[0];
    const double Ly = fr.pbc.L[1];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Ly <= 0.0 || Lz <= 0.0) {
        throw std::runtime_error("DensityXAnalyzer: invalid box lengths");
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
        throw std::runtime_error("DensityXAnalyzer: invalid bin/slice dimensions");
    }

    std::array<std::vector<int64_t>, 3> frame_counts{};
    for (int s = 0; s < 3; ++s) frame_counts[static_cast<size_t>(s)].assign(static_cast<size_t>(cfg_.nx), 0);

    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        const int sid = species_index_from_type(topo.mols[mid].type);
        if (sid < 0 || sid >= 3) continue;

        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) continue;

        Vec3d key = c.key_wrapped;
        fr.pbc.wrap_pos3(key);
        const double xw = key.v[0];
        const double zw = key.v[2];
        if (!in_range_pbc(zw, zminw, zmaxw)) continue;

        int ix = static_cast<int>(std::floor(xw / dx));
        if (ix < 0 || ix >= cfg_.nx) continue;
        frame_counts[static_cast<size_t>(sid)][static_cast<size_t>(ix)] += 1;
    }

    const double bin_volume = dx * Ly * zlen;
    for (int i = 0; i < cfg_.nx; ++i) {
        for (int s = 0; s < 3; ++s) {
            const int64_t c = frame_counts[static_cast<size_t>(s)][static_cast<size_t>(i)];
            const double rho = static_cast<double>(c) / bin_volume;
            rho_[static_cast<size_t>(s)].add(i, rho);
            cnt_[static_cast<size_t>(s)].add(i, static_cast<double>(c));
        }
    }

    nframes_ += 1;
}

void DensityXAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("DensityXAnalyzer: no frames processed");
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("DensityXAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "x_center_nm,"
           "rho_water_mean,rho_water_sem,rho_na_mean,rho_na_sem,rho_cl_mean,rho_cl_sem,"
           "count_water_mean,count_water_sem,count_na_mean,count_na_sem,count_cl_mean,count_cl_sem\n";

    for (int i = 0; i < cfg_.nx; ++i) {
        const double x_center = (static_cast<double>(i) + 0.5) * dx_;
        ofs << x_center << ","
            << rho_[species_index(Species::Water)].mean(i, nframes_) << ","
            << rho_[species_index(Species::Water)].sem(i, nframes_) << ","
            << rho_[species_index(Species::Na)].mean(i, nframes_) << ","
            << rho_[species_index(Species::Na)].sem(i, nframes_) << ","
            << rho_[species_index(Species::Cl)].mean(i, nframes_) << ","
            << rho_[species_index(Species::Cl)].sem(i, nframes_) << ","
            << cnt_[species_index(Species::Water)].mean(i, nframes_) << ","
            << cnt_[species_index(Species::Water)].sem(i, nframes_) << ","
            << cnt_[species_index(Species::Na)].mean(i, nframes_) << ","
            << cnt_[species_index(Species::Na)].sem(i, nframes_) << ","
            << cnt_[species_index(Species::Cl)].mean(i, nframes_) << ","
            << cnt_[species_index(Species::Cl)].sem(i, nframes_) << "\n";
    }
}

}  // namespace simio::analysis
