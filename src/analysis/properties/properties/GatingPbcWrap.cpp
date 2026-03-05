#include "simio/properties/GatingPbcWrap.hpp"

#include <cmath>

namespace simio::properties {

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

inline bool in_channel_z_slab(const Pbc3D& pbc, double z_cross, double z_center, double channel_height_z) {
    if (channel_height_z <= 0.0) return true;
    const double z_rel = pbc.wrap_delta(2, z_cross - z_center);
    return std::abs(z_rel) <= 0.5 * channel_height_z;
}

}  // namespace

GatingPbcWrapProperty::GatingPbcWrapProperty(const GatingPbcWrapConfig& cfg) : cfg_(cfg) {}

const char* GatingPbcWrapProperty::name() const { return "GatingPbcWrapProperty"; }

uint64_t GatingPbcWrapProperty::requires() const {
    return layered::field_mask(layered::Field::Layer0Bins) |
           layered::field_mask(layered::Field::KeyTimeWrapped) |
           layered::field_mask(layered::Field::InChannel);
}

void GatingPbcWrapProperty::compute_frame(const layered::FrameContext& ctx) {
    const size_t nmol = ctx.topo.mols.size();
    if (prev_key_time_wrapped_.size() != nmol) {
        prev_key_time_wrapped_.assign(nmol, Vec3d{});
        has_prev_.assign(nmol, 0);
        prev_in_channel_.assign(nmol, 0);
        frames_.clear();
        cumulative_ = GatingFrameCount{};
    }

    GatingFrameCount out{};

    for (size_t mid = 0; mid < nmol; ++mid) {
        const bool in_now = (ctx.layer1.is_in_channel[mid] != 0);
        const bool in_prev = (prev_in_channel_[mid] != 0);
        const bool include = include_molecule(cfg_.selection, in_now, in_prev);
        const Vec3d key_now = ctx.layer1.key_time_wrapped[mid];

        if (has_prev_[mid] && include) {
            const Vec3d key_prev = prev_key_time_wrapped_[mid];
            const Vec3d dr_raw = key_now - key_prev;
            const Vec3d dr_mi = ctx.frame.pbc.min_image(dr_raw);

            const CrossingResult crossing = detect_crossing(key_prev.v[0], dr_mi.v[0], 0.0, ctx.frame.pbc.L[0]);
            if (crossing.dir != CrossingDir::None) {
                const double z_cross = key_prev.v[2] + crossing.t * dr_mi.v[2];
                if (in_channel_z_slab(ctx.frame.pbc, z_cross, ctx.layer0.z_center,
                                      ctx.layer0.channel_height_z)) {
                    if (crossing.dir == CrossingDir::LeftToRight) {
                        ++out.n_right;
                    } else if (crossing.dir == CrossingDir::RightToLeft) {
                        ++out.n_left;
                    }
                }
            }
        }

        prev_key_time_wrapped_[mid] = key_now;
        has_prev_[mid] = 1;
        prev_in_channel_[mid] = in_now ? 1 : 0;
    }

    out.dn = out.n_right - out.n_left;
    cumulative_.n_left += out.n_left;
    cumulative_.n_right += out.n_right;
    cumulative_.dn = cumulative_.n_right - cumulative_.n_left;
    frames_.push_back(out);
}

void GatingPbcWrapProperty::finalize() {}

const std::vector<GatingFrameCount>& GatingPbcWrapProperty::frames() const { return frames_; }

const GatingFrameCount& GatingPbcWrapProperty::cumulative() const { return cumulative_; }

}  // namespace simio::properties
