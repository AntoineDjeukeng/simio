#include "simio/properties/DensityXZ.hpp"
#include "simio/analysis/intrinsics/context.hpp"
#include "simio/analysis/intrinsics/z_grid_cache.hpp"
#include "simio/runtime/cache.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace simio::properties {

namespace {

struct AxisIntervalMap {
    bool inside = false;
    double u = 0.0;
    double length = 0.0;
};

inline AxisIntervalMap map_on_pbc_interval(double x_wrapped, double mn, double mx, double L) {
    AxisIntervalMap out{};
    if (L <= 0.0) return out;

    if (mn <= mx) {
        out.length = mx - mn;
        if (out.length <= 0.0) return out;
        if (x_wrapped < mn || x_wrapped >= mx) return out;
        out.inside = true;
        out.u = x_wrapped - mn;
        return out;
    }

    out.length = (mx + L) - mn;
    if (out.length <= 0.0) return out;
    if (x_wrapped >= mn) {
        out.inside = true;
        out.u = x_wrapped - mn;
        return out;
    }
    if (x_wrapped < mx) {
        out.inside = true;
        out.u = (x_wrapped + L) - mn;
        return out;
    }
    return out;
}

inline int clamp_bin(double u, double d, int n) {
    int i = (int)std::floor(u / d);
    if (i < 0) i = 0;
    if (i >= n) i = n - 1;
    return i;
}

inline int type_bucket(MolType t) {
    if (t == MolType::Water) return 0;
    if (t == MolType::Cation) return 1;
    if (t == MolType::Anion) return 2;
    return 3;
}

}  // namespace

DensityXZProperty::DensityXZProperty(const DensityXZConfig& cfg) : cfg_(cfg) {}

DensityXZProperty::DensityXZProperty(const DensityXZConfig& cfg, simio::runtime::CacheStore& cache)
    : DensityXZProperty(cfg) {
    cache_ = &cache;
}

const char* DensityXZProperty::name() const { return "DensityXZProperty"; }

uint64_t DensityXZProperty::requires() const {
    return layered::field_mask(layered::Field::KeyWrapped) |
           layered::field_mask(layered::Field::InChannel);
}

void DensityXZProperty::compute_frame(const layered::FrameContext& ctx) {
    const size_t nmol = ctx.topo.mols.size();
    if (prev_in_channel_.size() != nmol) {
        prev_in_channel_.assign(nmol, 0);
    }

    DensityXZFrame out;
    out.nx = std::max(1, cfg_.nx);
    out.nz = std::max(1, cfg_.nz);

    const double Lx = ctx.frame.pbc.L[0];
    const double Ly = ctx.frame.pbc.L[1];
    const double Lz = ctx.frame.pbc.L[2];
    const double x_min_w = ctx.frame.pbc.wrap_pos(0, cfg_.roi.x_min);
    const double x_max_w = ctx.frame.pbc.wrap_pos(0, cfg_.roi.x_max);
    const double z_min_w = ctx.frame.pbc.wrap_pos(2, cfg_.roi.z_min);
    const double z_max_w = ctx.frame.pbc.wrap_pos(2, cfg_.roi.z_max);

    const AxisIntervalMap z_int = map_on_pbc_interval(z_min_w, z_min_w, z_max_w, Lz);

    double x_length = 0.0;
    if (cfg_.roi.x_mode == DensityXMode::FullBox) {
        x_length = Lx;
    } else {
        const AxisIntervalMap x_len = map_on_pbc_interval(x_min_w, x_min_w, x_max_w, Lx);
        x_length = x_len.length;
    }

    out.dx = (out.nx > 0 && x_length > 0.0) ? (x_length / (double)out.nx) : 0.0;
    if (out.nz > 0 && z_int.length > 0.0) {
        // z-grid is relative to the z-window start; zz_int.u uses the same convention.
        if (!cache_) {
            throw std::runtime_error("DensityXZProperty: cache is not wired");
        }
        simio::analysis::intrinsics::IntrinsicContext ictx{*cache_};
        const auto z_grid = simio::analysis::intrinsics::get_z_grid(ictx, 0.0, z_int.length, out.nz);
        out.dz = z_grid.dz;
    } else {
        out.dz = 0.0;
    }
    out.bin_volume = out.dx * Ly * out.dz;

    out.counts.assign((size_t)(out.nx * out.nz), 0);
    out.rho.assign((size_t)(out.nx * out.nz), 0.0);

    if (out.dx <= 0.0 || out.dz <= 0.0 || Ly <= 0.0) {
        frames_.push_back(std::move(out));
        return;
    }

    for (size_t mid = 0; mid < nmol; ++mid) {
        const bool in_now = (ctx.layer1.is_in_channel[mid] != 0);
        const bool in_prev = (prev_in_channel_[mid] != 0);
        const bool include = include_molecule(cfg_.selection, in_now, in_prev);

        if (include) {
            const Vec3d& r = ctx.layer1.key_wrapped[mid];
            const double xw = ctx.frame.pbc.wrap_pos(0, r.v[0]);
            const double zw = ctx.frame.pbc.wrap_pos(2, r.v[2]);

            AxisIntervalMap x_int{};
            if (cfg_.roi.x_mode == DensityXMode::FullBox) {
                x_int.inside = true;
                x_int.length = Lx;
                x_int.u = xw;
            } else {
                x_int = map_on_pbc_interval(xw, x_min_w, x_max_w, Lx);
            }
            const AxisIntervalMap zz_int = map_on_pbc_interval(zw, z_min_w, z_max_w, Lz);

            const bool in_roi = x_int.inside && zz_int.inside;
            if (in_roi) {
                ++out.selected_total;
                switch (type_bucket(ctx.topo.mols[mid].type)) {
                    case 0:
                        ++out.selected_water;
                        break;
                    case 1:
                        ++out.selected_cation;
                        break;
                    case 2:
                        ++out.selected_anion;
                        break;
                    default:
                        ++out.selected_other;
                        break;
                }

                const int ix = clamp_bin(x_int.u, out.dx, out.nx);
                const int iz = clamp_bin(zz_int.u, out.dz, out.nz);
                const int bin = iz * out.nx + ix;
                if (bin >= 0 && bin < out.nx * out.nz) {
                    ++out.counts[(size_t)bin];
                    ++out.binned_total;
                } else {
                    ++out.bin_oob_count;
                }
            }
        }

        prev_in_channel_[mid] = in_now ? 1 : 0;
    }

    if (cfg_.normalize_number_density && out.bin_volume > 0.0) {
        for (size_t i = 0; i < out.counts.size(); ++i) {
            out.rho[i] = (double)out.counts[i] / out.bin_volume;
        }
    }

    frames_.push_back(std::move(out));
}

void DensityXZProperty::finalize() {}

const std::vector<DensityXZFrame>& DensityXZProperty::frames() const { return frames_; }

}  // namespace simio::properties
