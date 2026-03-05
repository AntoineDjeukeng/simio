#pragma once

#include <cmath>
#include <cstdint>

#include "simio/simio.hpp"

namespace simio::regions {

struct ChannelRegionConfig {
    double channel_length_x = 0.0;
    double channel_height_z = 0.0;
    double gate_left_rel_z = 0.0;
    double gate_right_rel_z = 0.0;
};

class ChannelRegion {
  public:
    explicit ChannelRegion(const ChannelRegionConfig& cfg = {}) : cfg_(cfg) {}

    void set_config(const ChannelRegionConfig& cfg) { cfg_ = cfg; }
    const ChannelRegionConfig& config() const { return cfg_; }

    double x_rel(const Pbc3D& pbc, double x_wrapped, double x_center_wrapped) const {
        return pbc.wrap_delta(0, x_wrapped - x_center_wrapped);
    }

    double z_rel(const Pbc3D& pbc, double z_wrapped, double z_center_wrapped) const {
        return pbc.wrap_delta(2, z_wrapped - z_center_wrapped);
    }

    bool in_channel(double x_rel_value, double z_rel_value) const {
        const bool has_x = (cfg_.channel_length_x > 0.0);
        const bool has_z = (cfg_.channel_height_z > 0.0);
        if (!has_x && !has_z) return false;
        if (has_x && std::fabs(x_rel_value) > 0.5 * cfg_.channel_length_x) return false;
        if (has_z && std::fabs(z_rel_value) > 0.5 * cfg_.channel_height_z) return false;
        return true;
    }

    int8_t gate_side(double z_rel_value) const {
        if (cfg_.gate_right_rel_z <= cfg_.gate_left_rel_z) return 0;
        if (z_rel_value < cfg_.gate_left_rel_z) return -1;
        if (z_rel_value > cfg_.gate_right_rel_z) return 1;
        return 0;
    }

  private:
    ChannelRegionConfig cfg_{};
};

}  // namespace simio::regions

