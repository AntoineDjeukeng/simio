#include "density_z_in_x_channel.hpp"

#include "simio/analysis/intrinsics/channel_roi.hpp"
#include "simio/analysis/intrinsics/context.hpp"
#include "simio/analysis/intrinsics/z_grid_cache.hpp"
#include "simio/runtime/cache.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace simio::analysis {

DensityZInXChannelAnalyzer::DensityZInXChannelAnalyzer(const DensityZInXChannelConfig& cfg)
    : cfg_(cfg) {
    if (cfg_.nz <= 0) throw std::runtime_error("DensityZInXChannelAnalyzer: nz must be > 0");
    for (int s = 0; s < 3; ++s) {
        rho_[static_cast<size_t>(s)].init(cfg_.nz);
        cnt_[static_cast<size_t>(s)].init(cfg_.nz);
    }
    rho_q_.init(cfg_.nz);
}

DensityZInXChannelAnalyzer::DensityZInXChannelAnalyzer(const DensityZInXChannelConfig& cfg,
                                                       simio::runtime::CacheStore& cache)
    : DensityZInXChannelAnalyzer(cfg) {
    cache_ = &cache;
    // Build/get a relative z-grid once per run (unit interval [0,1]).
    const auto grid = simio::analysis::intrinsics::get_or_build_z_grid(*cache_, 0.0, 1.0, cfg_.nz);
    dz_ = grid.dz;
    has_cached_rel_grid_ = true;
}

void DensityZInXChannelAnalyzer::process_frame(const Topology& topo,
                                               const Frame& fr,
                                               const std::vector<MolState>& ms) {
    const double Lx = fr.pbc.L[0];
    const double Ly = fr.pbc.L[1];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Ly <= 0.0 || Lz <= 0.0) {
        throw std::runtime_error("DensityZInXChannelAnalyzer: invalid box lengths");
    }

    if (!has_ref_box_) {
        has_ref_box_ = true;
        Lz_ref_ = Lz;
        dz_ = has_cached_rel_grid_ ? (dz_ * Lz_ref_) : (Lz_ref_ / static_cast<double>(cfg_.nz));
    }
    if (!has_roi_) {
        if (!cache_) {
            throw std::runtime_error("DensityZInXChannelAnalyzer: cache is not wired");
        }
        simio::analysis::intrinsics::IntrinsicContext ictx{*cache_};
        roi_ = simio::analysis::intrinsics::get_channel_roi_x(ictx, cfg_.xmin, cfg_.xmax, Lx);
        xlen_ = roi_.xlen;
        has_roi_ = true;
    }

    const double dz = dz_;
    const double xlen = xlen_;
    if (dz <= 0.0 || xlen <= 0.0) {
        throw std::runtime_error("DensityZInXChannelAnalyzer: invalid bin/channel dimensions");
    }

    std::array<std::vector<int64_t>, 3> frame_counts{};
    for (int s = 0; s < 3; ++s) {
        frame_counts[static_cast<size_t>(s)].assign(static_cast<size_t>(cfg_.nz), 0);
    }

    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        const int sid = species_index_from_type(topo.mols[mid].type);
        if (sid < 0 || sid >= 3) continue;

        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) continue;

        Vec3d key = c.key_wrapped;
        fr.pbc.wrap_pos3(key);

        const double xw = key.v[0];
        const double zw = key.v[2];
        if (!roi_.contains_x(xw)) continue;

        int iz = static_cast<int>(std::floor(zw / dz));
        if (iz < 0) continue;
        if (iz >= cfg_.nz) iz = cfg_.nz - 1;

        frame_counts[static_cast<size_t>(sid)][static_cast<size_t>(iz)] += 1;
    }

    const double bin_volume = dz * Ly * xlen;
    for (int i = 0; i < cfg_.nz; ++i) {
        double rho_na = 0.0;
        double rho_cl = 0.0;
        for (int s = 0; s < 3; ++s) {
            const int64_t c = frame_counts[static_cast<size_t>(s)][static_cast<size_t>(i)];
            const double rho = static_cast<double>(c) / bin_volume;
            rho_[static_cast<size_t>(s)].add(i, rho);
            cnt_[static_cast<size_t>(s)].add(i, static_cast<double>(c));
            if (s == species_index(Species::Na)) rho_na = rho;
            if (s == species_index(Species::Cl)) rho_cl = rho;
        }
        // Charge density in elementary-charge units: e/nm^3.
        rho_q_.add(i, rho_na - rho_cl);
    }

    nframes_ += 1;
}

void DensityZInXChannelAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("DensityZInXChannelAnalyzer: no frames processed");
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("DensityZInXChannelAnalyzer: failed to open output CSV: " +
                                 path);
    }

    ofs << std::setprecision(12);
    ofs << "z_center_nm,"
           "rho_water_mean,rho_water_sem,rho_na_mean,rho_na_sem,rho_cl_mean,rho_cl_sem,"
           "rho_q_e_mean,rho_q_e_sem,"
           "count_water_mean,count_water_sem,count_na_mean,count_na_sem,count_cl_mean,count_cl_sem\n";

    for (int i = 0; i < cfg_.nz; ++i) {
        const double z_center = (static_cast<double>(i) + 0.5) * dz_;
        ofs << z_center << ","
            << rho_[species_index(Species::Water)].mean(i, nframes_) << ","
            << rho_[species_index(Species::Water)].sem(i, nframes_) << ","
            << rho_[species_index(Species::Na)].mean(i, nframes_) << ","
            << rho_[species_index(Species::Na)].sem(i, nframes_) << ","
            << rho_[species_index(Species::Cl)].mean(i, nframes_) << ","
            << rho_[species_index(Species::Cl)].sem(i, nframes_) << ","
            << rho_q_.mean(i, nframes_) << "," << rho_q_.sem(i, nframes_) << ","
            << cnt_[species_index(Species::Water)].mean(i, nframes_) << ","
            << cnt_[species_index(Species::Water)].sem(i, nframes_) << ","
            << cnt_[species_index(Species::Na)].mean(i, nframes_) << ","
            << cnt_[species_index(Species::Na)].sem(i, nframes_) << ","
            << cnt_[species_index(Species::Cl)].mean(i, nframes_) << ","
            << cnt_[species_index(Species::Cl)].sem(i, nframes_) << "\n";
    }
}

}  // namespace simio::analysis
