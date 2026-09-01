#pragma once

#include "reactor_setup.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace simio::middle_reservoir {

class MiddleReservoirMonitor {
  public:
    explicit MiddleReservoirMonitor(const ReactorSetup& setup);

    void process_frame(const simio::Frame& frame, int frame_index);
    void write_csv(const std::string& path) const;

  private:
    struct GateTally {
        std::array<int64_t, 3> entered{0, 0, 0};
        std::array<int64_t, 3> exited{0, 0, 0};
    };

    struct FrameRow {
        int frame_index = 0;
        int64_t step = 0;
        double time_ps = 0.0;
        // Regions: left reservoir, middle reservoir, right reservoir, channels.
        std::array<std::array<int64_t, 3>, 4> region_count{};
        GateTally left;
        GateTally right;
        std::array<int64_t, 3> left_cumulative_net{0, 0, 0};
        std::array<int64_t, 3> right_cumulative_net{0, 0, 0};
    };

    simio::Vec3d molecule_key(const simio::Frame& frame, const simio::MolSpan& molecule) const;

    const ReactorSetup& setup_;
    std::vector<simio::Vec3d> previous_key_;
    std::vector<uint8_t> has_previous_;
    std::array<int64_t, 3> left_cumulative_net_{0, 0, 0};
    std::array<int64_t, 3> right_cumulative_net_{0, 0, 0};
    std::vector<FrameRow> rows_;
};

void write_setup_json(const ReactorSetup& setup,
                      const std::string& gro_path,
                      const std::string& report_path,
                      const std::string& output_path);

}  // namespace simio::middle_reservoir
