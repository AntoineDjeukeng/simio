#include "gating_flux.hpp"

#include "simio/analysis/intrinsics/context.hpp"
#include "simio/analysis/intrinsics/in_channel_mask.hpp"
#include "simio/runtime/cache.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace simio::analysis {

namespace {

constexpr double kCrossEps = 1e-12;

enum class CrossingDir : uint8_t { None = 0, LeftToRight = 1, RightToLeft = 2 };

struct CrossingResult {
    CrossingDir dir = CrossingDir::None;
    double t = 0.0;
};

inline double nearest_plane_image(double x_prev, double x_plane, double Lx) {
    if (Lx <= 0.0) return x_plane;
    const double k = std::round((x_prev - x_plane) / Lx);
    return x_plane + k * Lx;
}

inline CrossingResult detect_crossing(double x_prev, double dx_mi, double x_plane, double Lx) {
    const double x_plane_near = nearest_plane_image(x_prev, x_plane, Lx);
    const double d0 = x_prev - x_plane_near;
    const double d1 = (x_prev + dx_mi) - x_plane_near;

    if (d0 < 0.0 && d1 >= 0.0) {
        const double denom = d1 - d0;
        if (std::abs(denom) <= kCrossEps) return {};
        double t = -d0 / denom;
        if (t < -kCrossEps || t > 1.0 + kCrossEps) return {};
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        return CrossingResult{CrossingDir::LeftToRight, t};
    }

    if (d0 >= 0.0 && d1 < 0.0) {
        const double denom = d1 - d0;
        if (std::abs(denom) <= kCrossEps) return {};
        double t = -d0 / denom;
        if (t < -kCrossEps || t > 1.0 + kCrossEps) return {};
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        return CrossingResult{CrossingDir::RightToLeft, t};
    }

    return {};
}

}  // namespace

GatingSelection parse_gating_selection(const std::string& value) {
    std::string v = value;
    std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (v == "0" || v == "all") return GatingSelection::All;
    if (v == "1" || v == "now" || v == "inchannelnow") return GatingSelection::InChannelNow;
    if (v == "2" || v == "both" || v == "inchannelboth") return GatingSelection::InChannelBoth;

    throw std::runtime_error("Invalid gating selection: " + value + " (use all|now|both or 0|1|2)");
}

GatingFluxAnalyzer::GatingFluxAnalyzer(const GatingFluxConfig& cfg) : cfg_(cfg) {}

GatingFluxAnalyzer::GatingFluxAnalyzer(const GatingFluxConfig& cfg, simio::runtime::CacheStore& cache)
    : GatingFluxAnalyzer(cfg) {
    cache_ = &cache;
}

void GatingFluxAnalyzer::process_frame(const Topology& topo, const Frame& fr,
                                       const std::vector<MolState>& ms, int frame_idx) {
    const size_t nmol = topo.mols.size();
    if (prev_key_.size() != nmol) {
        prev_key_.assign(nmol, Vec3d{});
        has_prev_.assign(nmol, 0);
        prev_in_channel_.assign(nmol, 0);
        xw_tmp_.assign(nmol, 0.0);
        frame_counter_ = 0;
        rows_.clear();
        center_cum_ = PlaneTally{};
        seam_cum_ = PlaneTally{};
        nframes_ = 0;
    }

    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Lz <= 0.0) throw std::runtime_error("GatingFluxAnalyzer: invalid box lengths");

    if (!cache_) {
        static thread_local simio::runtime::CacheStore fallback_cache;
        cache_ = &fallback_cache;
    }

    const double xminw = fr.pbc.wrap_pos(0, cfg_.xmin);
    const double xmaxw = fr.pbc.wrap_pos(0, cfg_.xmax);
    const double xlen = interval_length_pbc(xminw, xmaxw, Lx);
    const bool whole_x = (xlen <= 0.0);
    const double zminw = fr.pbc.wrap_pos(2, cfg_.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg_.zmax);
    const double zlen = interval_length_pbc(zminw, zmaxw, Lz);
    if (zlen <= 0.0) throw std::runtime_error("GatingFluxAnalyzer: invalid z interval");

    const double x_center =
        whole_x ? 0.5 * Lx : wrapped_interval_midpoint(fr.pbc, 0, xminw, xmaxw, Lx);

    if (xw_tmp_.size() != nmol) xw_tmp_.assign(nmol, 0.0);
    for (size_t mid = 0; mid < nmol; ++mid) {
        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) {
            xw_tmp_[mid] = 0.0;
            continue;
        }
        Vec3d key_now = c.key_wrapped;
        fr.pbc.wrap_pos3(key_now);
        const double xw = std::fmod(key_now.v[0], Lx);
        xw_tmp_[mid] = (xw < 0.0) ? (xw + Lx) : xw;
    }
    const std::int64_t frame_id = frame_counter_++;
    simio::analysis::intrinsics::IntrinsicContext mask_ctx{*cache_};
    const auto mask = simio::analysis::intrinsics::get_in_channel_mask_x(
        mask_ctx, frame_id, xw_tmp_.data(), xw_tmp_.size(), cfg_.xmin, cfg_.xmax, Lx);
    if (mask.in_channel.size() != nmol) {
        throw std::runtime_error("GatingFluxAnalyzer: in_channel_mask size mismatch");
    }

    auto in_channel_now = [&](const Vec3d& key, size_t mid) -> bool {
        const bool x_ok = whole_x ? true : (mask.in_channel[mid] != 0u);
        const bool z_ok = map_on_pbc_interval(key.v[2], zminw, zmaxw, Lz).inside;
        return x_ok && z_ok;
    };

    auto in_z_slab = [&](double z_cross) -> bool {
        const double zc = fr.pbc.wrap_pos(2, z_cross);
        return in_range_pbc(zc, zminw, zmaxw);
    };

    auto include_molecule = [&](bool in_now, bool in_prev) -> bool {
        if (cfg_.selection == GatingSelection::All) return true;
        if (cfg_.selection == GatingSelection::InChannelBoth) return in_now && in_prev;
        return in_now;
    };

    PlaneTally center_frame{};
    PlaneTally seam_frame{};
    FrameRow row{};

    auto add_crossing = [](PlaneTally& tally, CrossingDir dir, int sid) {
        if (dir == CrossingDir::LeftToRight) {
            tally.n_right += 1;
            if (sid >= 0 && sid < 3) tally.n_right_species[static_cast<size_t>(sid)] += 1;
        } else if (dir == CrossingDir::RightToLeft) {
            tally.n_left += 1;
            if (sid >= 0 && sid < 3) tally.n_left_species[static_cast<size_t>(sid)] += 1;
        }
    };

    for (size_t mid = 0; mid < nmol; ++mid) {
        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) continue;

        Vec3d key_now = c.key_wrapped;
        fr.pbc.wrap_pos3(key_now);

        const bool in_now = in_channel_now(key_now, mid);
        const bool in_prev = (prev_in_channel_[mid] != 0);
        const bool include = include_molecule(in_now, in_prev);

        if (has_prev_[mid] && include) {
            const Vec3d key_prev = prev_key_[mid];
            const Vec3d dr_mi = fr.pbc.min_image(key_now - key_prev);
            const int sid = species_index_from_type(topo.mols[mid].type);
            const CrossingResult cross_center = detect_crossing(key_prev.v[0], dr_mi.v[0], x_center, Lx);
            if (cross_center.dir != CrossingDir::None) {
                const double z_cross = key_prev.v[2] + cross_center.t * dr_mi.v[2];
                if (in_z_slab(z_cross)) {
                    add_crossing(center_frame, cross_center.dir, sid);
                    add_crossing(center_cum_, cross_center.dir, sid);
                }
            }

            const CrossingResult cross_seam = detect_crossing(key_prev.v[0], dr_mi.v[0], 0.0, Lx);
            if (cross_seam.dir != CrossingDir::None) {
                const double z_cross = key_prev.v[2] + cross_seam.t * dr_mi.v[2];
                if (in_z_slab(z_cross)) {
                    add_crossing(seam_frame, cross_seam.dir, sid);
                    add_crossing(seam_cum_, cross_seam.dir, sid);
                }
            }
        }

        prev_key_[mid] = key_now;
        has_prev_[mid] = 1;
        prev_in_channel_[mid] = in_now ? 1 : 0;
    }

    row.frame_idx = frame_idx;
    row.step = fr.step;
    row.time_ps = fr.time_ps;
    row.center_frame = center_frame;
    row.center_cum = center_cum_;
    row.seam_frame = seam_frame;
    row.seam_cum = seam_cum_;
    rows_.push_back(row);

    nframes_ += 1;
}

void GatingFluxAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("GatingFluxAnalyzer: no frames processed");

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("GatingFluxAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "frame_idx,step,time_ps,"
           "center_water_left,center_water_right,center_water_cum_left,center_water_cum_right,"
           "seam_water_left,seam_water_right,seam_water_cum_left,seam_water_cum_right,"
           "center_na_left,center_na_right,center_na_cum_left,center_na_cum_right,"
           "seam_na_left,seam_na_right,seam_na_cum_left,seam_na_cum_right,"
           "center_cl_left,center_cl_right,center_cl_cum_left,center_cl_cum_right,"
           "seam_cl_left,seam_cl_right,seam_cl_cum_left,seam_cl_cum_right\n";

    for (const FrameRow& r : rows_) {
        ofs << r.frame_idx << "," << r.step << "," << r.time_ps << ","
            << r.center_frame.n_left_species[0] << "," << r.center_frame.n_right_species[0] << ","
            << r.center_cum.n_left_species[0] << "," << r.center_cum.n_right_species[0] << ","
            << r.seam_frame.n_left_species[0] << "," << r.seam_frame.n_right_species[0] << ","
            << r.seam_cum.n_left_species[0] << "," << r.seam_cum.n_right_species[0] << ","
            << r.center_frame.n_left_species[1] << "," << r.center_frame.n_right_species[1] << ","
            << r.center_cum.n_left_species[1] << "," << r.center_cum.n_right_species[1] << ","
            << r.seam_frame.n_left_species[1] << "," << r.seam_frame.n_right_species[1] << ","
            << r.seam_cum.n_left_species[1] << "," << r.seam_cum.n_right_species[1] << ","
            << r.center_frame.n_left_species[2] << "," << r.center_frame.n_right_species[2] << ","
            << r.center_cum.n_left_species[2] << "," << r.center_cum.n_right_species[2] << ","
            << r.seam_frame.n_left_species[2] << "," << r.seam_frame.n_right_species[2] << ","
            << r.seam_cum.n_left_species[2] << "," << r.seam_cum.n_right_species[2] << "\n";
    }
}

}  // namespace simio::analysis
