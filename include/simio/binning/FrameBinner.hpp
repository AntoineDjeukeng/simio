#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "simio/simio.hpp"

namespace simio::binning {

struct FrameBinnerConfig {
    int nx = 1;
    int nz = 1;
    double x_center_wrapped = std::numeric_limits<double>::quiet_NaN();
    double z_center_wrapped = std::numeric_limits<double>::quiet_NaN();
};

struct FrameBinnerMeta {
    int nx = 1;
    int nz = 1;
    double Lx = 1.0;
    double Lz = 1.0;
    double dx = 1.0;
    double dz = 1.0;
    double x_center = 0.0;
    double z_center = 0.0;
};

class FrameBinner {
  public:
    explicit FrameBinner(const FrameBinnerConfig& cfg = {}) : cfg_(cfg) {}

    void set_config(const FrameBinnerConfig& cfg) { cfg_ = cfg; }
    const FrameBinnerConfig& config() const { return cfg_; }

    FrameBinnerMeta make_meta(const Pbc3D& pbc) const {
        FrameBinnerMeta out;
        out.nx = std::max(1, cfg_.nx);
        out.nz = std::max(1, cfg_.nz);
        out.Lx = pbc.L[0];
        out.Lz = pbc.L[2];
        out.dx = out.Lx / out.nx;
        out.dz = out.Lz / out.nz;
        out.x_center = resolve_center(pbc, 0, cfg_.x_center_wrapped);
        out.z_center = resolve_center(pbc, 2, cfg_.z_center_wrapped);
        return out;
    }

    void map_xz(const Pbc3D& pbc, const FrameBinnerMeta& meta, const Vec3d& r_wrapped, int& ix, int& iz,
                int& bin_id) const {
        const double x_shift = pbc.wrap_pos(0, r_wrapped.v[0] - meta.x_center + 0.5 * meta.Lx);
        const double z_shift = pbc.wrap_pos(2, r_wrapped.v[2] - meta.z_center + 0.5 * meta.Lz);

        ix = (int)std::floor(x_shift / meta.dx);
        iz = (int)std::floor(z_shift / meta.dz);
        if (ix < 0) ix = 0;
        if (iz < 0) iz = 0;
        if (ix >= meta.nx) ix = meta.nx - 1;
        if (iz >= meta.nz) iz = meta.nz - 1;

        bin_id = iz * meta.nx + ix;
    }

  private:
    static double resolve_center(const Pbc3D& pbc, int axis, double c) {
        if (!std::isfinite(c)) return 0.5 * pbc.L[axis];
        return pbc.wrap_pos(axis, c);
    }

    FrameBinnerConfig cfg_{};
};

}  // namespace simio::binning

