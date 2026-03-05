#include "dipole_z_in_x_channel.hpp"

#include "simio/analysis/intrinsics/channel_roi.hpp"
#include "simio/analysis/intrinsics/context.hpp"
#include "simio/analysis/intrinsics/in_channel_mask.hpp"
#include "simio/analysis/intrinsics/z_grid_cache.hpp"
#include "simio/runtime/cache.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace simio::analysis {

DipoleZInXChannelAnalyzer::DipoleZInXChannelAnalyzer(const DipoleZInXChannelConfig& cfg)
    : cfg_(cfg) {
    if (cfg_.nz <= 0) throw std::runtime_error("DipoleZInXChannelAnalyzer: nz must be > 0");
    mux_.init(cfg_.nz);
    mux_fold_.init(cfg_.nz);
    muz_.init(cfg_.nz);
    cnt_.init(cfg_.nz);
}

DipoleZInXChannelAnalyzer::DipoleZInXChannelAnalyzer(const DipoleZInXChannelConfig& cfg,
                                                     simio::runtime::CacheStore& cache)
    : DipoleZInXChannelAnalyzer(cfg) {
    cache_ = &cache;
    // Build/get a relative z-grid once per run (unit interval [0,1]).
    const auto grid = simio::analysis::intrinsics::get_or_build_z_grid(*cache_, 0.0, 1.0, cfg_.nz);
    dz_ = grid.dz;
    has_cached_rel_grid_ = true;
}

void DipoleZInXChannelAnalyzer::process_frame(const Topology& topo,
                                              const Frame& fr,
                                              const std::vector<MolState>& ms) {
    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Lz <= 0.0) {
        throw std::runtime_error("DipoleZInXChannelAnalyzer: invalid box lengths");
    }

    if (!has_ref_box_) {
        has_ref_box_ = true;
        Lz_ref_ = Lz;
        dz_ = has_cached_rel_grid_ ? (dz_ * Lz_ref_) : (Lz_ref_ / static_cast<double>(cfg_.nz));
    }
    if (!has_roi_) {
        if (!cache_) {
            throw std::runtime_error("DipoleZInXChannelAnalyzer: cache is not wired");
        }
        simio::analysis::intrinsics::IntrinsicContext ictx{*cache_};
        roi_ = simio::analysis::intrinsics::get_channel_roi_x(ictx, cfg_.xmin, cfg_.xmax, Lx);
        xlen_ = roi_.xlen;
        has_roi_ = true;
    }

    const double dz = dz_;
    const double xlen = xlen_;
    if (dz <= 0.0 || xlen <= 0.0) {
        throw std::runtime_error("DipoleZInXChannelAnalyzer: invalid bin/channel dimensions");
    }

    std::vector<double> frame_sum_mux(static_cast<size_t>(cfg_.nz), 0.0);
    std::vector<double> frame_sum_mux_fold(static_cast<size_t>(cfg_.nz), 0.0);
    std::vector<double> frame_sum_muz(static_cast<size_t>(cfg_.nz), 0.0);
    std::vector<int64_t> frame_count(static_cast<size_t>(cfg_.nz), 0);

    const size_t nmol = topo.mols.size();
    if (xw_tmp_.size() != nmol) xw_tmp_.assign(nmol, 0.0);
    for (size_t mid = 0; mid < nmol; ++mid) {
        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) {
            xw_tmp_[mid] = 0.0;
            continue;
        }
        Vec3d key = c.key_wrapped;
        fr.pbc.wrap_pos3(key);
        xw_tmp_[mid] = key.v[0];
    }

    const std::int64_t frame_id = frame_counter_++;
    simio::analysis::intrinsics::IntrinsicContext ictx{*cache_};
    const auto mask = simio::analysis::intrinsics::get_in_channel_mask_x(
        ictx,
        frame_id,
        xw_tmp_.data(),
        xw_tmp_.size(),
        cfg_.xmin,
        cfg_.xmax,
        Lx);
    if (mask.in_channel.size() != nmol) {
        throw std::runtime_error("DipoleZInXChannelAnalyzer: in_channel_mask size mismatch");
    }

    constexpr double kNormEps = 1e-14;
    for (size_t mid = 0; mid < nmol; ++mid) {
        const MolSpan& m = topo.mols[mid];
        if (m.type != MolType::Water || m.natoms < 3) continue;

        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) continue;

        Vec3d key = c.key_wrapped;
        fr.pbc.wrap_pos3(key);

        const double xw = key.v[0];
        const double zw = key.v[2];
        if (!std::isfinite(xw) || !std::isfinite(zw)) continue;
        if (mask.in_channel[mid] == 0u) continue;
        const double xu = roi_.map_x_to_channel(xw);

        int iz = static_cast<int>(std::floor(zw / dz));
        if (iz < 0) continue;
        if (iz >= cfg_.nz) iz = cfg_.nz - 1;

        Vec3d rO{};
        Vec3d h1u{};
        Vec3d h2u{};
        if ((c.flags & MolCache::HAS_REF) && (c.flags & MolCache::HAS_SITES)) {
            rO = c.ref_wrapped;
            h1u = c.sites_u[1];
            h2u = c.sites_u[2];
        } else {
            const Vec3d rH1 = fr.atoms.pos(m.first + 1);
            const Vec3d rH2 = fr.atoms.pos(m.first + 2);
            rO = fr.atoms.pos(m.first + 0);
            fr.pbc.wrap_pos3(rO);
            h1u = pbc_unwrap_to_ref(fr.pbc, rO, rH1);
            h2u = pbc_unwrap_to_ref(fr.pbc, rO, rH2);
        }

        Vec3d mu = (h1u - rO) + (h2u - rO);
        const double mnorm = norm3(mu);
        if (mnorm <= kNormEps) continue;
        mu = mu / mnorm;
        const double mux_raw = mu.v[0];
        const bool in_left_half = (xu < 0.5 * xlen);
        const double mux_fold = in_left_half ? -mux_raw : mux_raw;

        frame_sum_mux[static_cast<size_t>(iz)] += mux_raw;
        frame_sum_mux_fold[static_cast<size_t>(iz)] += mux_fold;
        frame_sum_muz[static_cast<size_t>(iz)] += mu.v[2];
        frame_count[static_cast<size_t>(iz)] += 1;
    }

    for (int i = 0; i < cfg_.nz; ++i) {
        const int64_t c = frame_count[static_cast<size_t>(i)];
        const double mean_mux = (c > 0) ? (frame_sum_mux[static_cast<size_t>(i)] / static_cast<double>(c))
                                        : 0.0;
        const double mean_mux_fold =
            (c > 0) ? (frame_sum_mux_fold[static_cast<size_t>(i)] / static_cast<double>(c)) : 0.0;
        const double mean_muz = (c > 0) ? (frame_sum_muz[static_cast<size_t>(i)] / static_cast<double>(c))
                                        : 0.0;

        mux_.add(i, mean_mux);
        mux_fold_.add(i, mean_mux_fold);
        muz_.add(i, mean_muz);
        cnt_.add(i, static_cast<double>(c));
    }

    nframes_ += 1;
}

void DipoleZInXChannelAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("DipoleZInXChannelAnalyzer: no frames processed");
    std::ofstream ofs(path);
    if (!ofs) {
        throw std::runtime_error("DipoleZInXChannelAnalyzer: failed to open output CSV: " + path);
    }

    ofs << std::setprecision(12);
    ofs << "z_center_nm,mux_mean,mux_sem,mux_fold_mean,mux_fold_sem,muz_mean,muz_sem,count_mean,count_sem\n";

    for (int i = 0; i < cfg_.nz; ++i) {
        const double z_center = (static_cast<double>(i) + 0.5) * dz_;
        ofs << z_center << "," << mux_.mean(i, nframes_) << "," << mux_.sem(i, nframes_) << ","
            << mux_fold_.mean(i, nframes_) << "," << mux_fold_.sem(i, nframes_) << ","
            << muz_.mean(i, nframes_) << "," << muz_.sem(i, nframes_) << ","
            << cnt_.mean(i, nframes_) << "," << cnt_.sem(i, nframes_) << "\n";
    }
}

}  // namespace simio::analysis
