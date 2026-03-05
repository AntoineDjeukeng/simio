#include "dipole_x.hpp"

#include "simio/analysis/intrinsics/x_grid_cache.hpp"
#include "simio/runtime/cache.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace simio::analysis {

DipoleXAnalyzer::DipoleXAnalyzer(const DipoleXConfig& cfg) : cfg_(cfg) {
    if (cfg_.nx <= 0) throw std::runtime_error("DipoleXAnalyzer: nx must be > 0");
    x_centers_.assign(static_cast<size_t>(cfg_.nx), 0.0);
    muz_raw_.init(cfg_.nx);
    muz_fold_.init(cfg_.nx);
    mux_.init(cfg_.nx);
    cnt_.init(cfg_.nx);
}

DipoleXAnalyzer::DipoleXAnalyzer(const DipoleXConfig& cfg, simio::runtime::CacheStore& cache)
    : DipoleXAnalyzer(cfg) {
    cache_ = &cache;
    // Build/get a relative x-grid once per run (unit interval [0,1]).
    const auto grid = simio::analysis::intrinsics::get_or_build_x_grid(*cache_, 0.0, 1.0, cfg_.nx);
    x_centers_rel_ = grid.centers_rel;
    dx_ = grid.dx;
    has_cached_rel_grid_ = true;
}

void DipoleXAnalyzer::process_frame(const Topology& topo, const Frame& fr,
                                    const std::vector<MolState>& ms) {
    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Lz <= 0.0) {
        throw std::runtime_error("DipoleXAnalyzer: invalid box lengths");
    }

    const bool whole_x = (cfg_.xmax <= cfg_.xmin);
    const double xminw = whole_x ? 0.0 : fr.pbc.wrap_pos(0, cfg_.xmin);
    const double xmaxw = whole_x ? 0.0 : fr.pbc.wrap_pos(0, cfg_.xmax);
    const double zminw = fr.pbc.wrap_pos(2, cfg_.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg_.zmax);

    const double xlen = whole_x ? Lx : interval_length_pbc(xminw, xmaxw, Lx);
    const double zlen = interval_length_pbc(zminw, zmaxw, Lz);
    if (xlen <= 0.0 || zlen <= 0.0) {
        throw std::runtime_error("DipoleXAnalyzer: invalid ROI length");
    }

    if (!has_x_centers_) {
        if (has_cached_rel_grid_) {
            dx_ *= xlen;
        } else {
            dx_ = xlen / static_cast<double>(cfg_.nx);
        }
        if (dx_ <= 0.0) throw std::runtime_error("DipoleXAnalyzer: invalid x bin width");

        has_x_centers_ = true;
        const double x0 = whole_x ? 0.0 : xminw;
        for (int i = 0; i < cfg_.nx; ++i) {
            const double xc_rel = has_cached_rel_grid_ ? (x_centers_rel_[static_cast<size_t>(i)] * xlen)
                                                       : ((static_cast<double>(i) + 0.5) * dx_);
            x_centers_[static_cast<size_t>(i)] = fr.pbc.wrap_pos(0, x0 + xc_rel);
        }
    }
    const double dx = dx_;

    std::vector<double> frame_sum_muz_raw(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_sum_muz_fold(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<double> frame_sum_mux(static_cast<size_t>(cfg_.nx), 0.0);
    std::vector<int64_t> frame_count(static_cast<size_t>(cfg_.nx), 0);

    constexpr double kNormEps = 1e-14;

    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        const MolSpan& m = topo.mols[mid];
        if (m.type != MolType::Water || m.natoms < 3) continue;

        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) continue;

        Vec3d key = c.key_wrapped;
        fr.pbc.wrap_pos3(key);

        IntervalMap xm{};
        if (whole_x) {
            xm.inside = true;
            xm.u = key.v[0];
            xm.length = Lx;
        } else {
            xm = map_on_pbc_interval(key.v[0], xminw, xmaxw, Lx);
        }
        const IntervalMap zm = map_on_pbc_interval(key.v[2], zminw, zmaxw, Lz);
        if (!(xm.inside && zm.inside)) continue;

        int ix = static_cast<int>(std::floor(xm.u / dx));
        if (ix < 0 || ix >= cfg_.nx) continue;

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

        const double muz_raw = mu.v[2];
        const double mux = mu.v[0];
        const bool lower_half = (zm.u < 0.5 * zlen);
        const double muz_fold = lower_half ? muz_raw : -muz_raw;

        frame_sum_muz_raw[static_cast<size_t>(ix)] += muz_raw;
        frame_sum_muz_fold[static_cast<size_t>(ix)] += muz_fold;
        frame_sum_mux[static_cast<size_t>(ix)] += mux;
        frame_count[static_cast<size_t>(ix)] += 1;
    }

    for (int i = 0; i < cfg_.nx; ++i) {
        const int64_t c = frame_count[static_cast<size_t>(i)];
        const double mean_muz_raw = (c > 0) ? (frame_sum_muz_raw[static_cast<size_t>(i)] / static_cast<double>(c))
                                            : 0.0;
        const double mean_muz_fold =
            (c > 0) ? (frame_sum_muz_fold[static_cast<size_t>(i)] / static_cast<double>(c)) : 0.0;
        const double mean_mux = (c > 0) ? (frame_sum_mux[static_cast<size_t>(i)] / static_cast<double>(c)) : 0.0;

        muz_raw_.add(i, mean_muz_raw);
        muz_fold_.add(i, mean_muz_fold);
        mux_.add(i, mean_mux);
        cnt_.add(i, static_cast<double>(c));
    }

    nframes_ += 1;
}

void DipoleXAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("DipoleXAnalyzer: no frames processed");
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("DipoleXAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "x_center_nm,muz_mean,muz_sem,muz_fold_mean,muz_fold_sem,mux_mean,mux_sem,count_mean,count_sem\n";

    for (size_t i = 0; i < x_centers_.size(); ++i) {
        const int ix = static_cast<int>(i);
        ofs << x_centers_[i] << "," << muz_raw_.mean(ix, nframes_) << "," << muz_raw_.sem(ix, nframes_)
            << "," << muz_fold_.mean(ix, nframes_) << "," << muz_fold_.sem(ix, nframes_) << ","
            << mux_.mean(ix, nframes_) << "," << mux_.sem(ix, nframes_) << "," << cnt_.mean(ix, nframes_)
            << "," << cnt_.sem(ix, nframes_) << "\n";
    }
}

}  // namespace simio::analysis
