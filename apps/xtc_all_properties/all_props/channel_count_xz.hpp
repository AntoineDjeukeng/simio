#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "common.hpp"
#include "simio/analysis/intrinsics/channel_roi.hpp"

namespace simio::runtime {
class CacheStore;
}

namespace simio::analysis {

struct ChannelCountXZConfig {
    double xmin = 7.11;
    double xmax = 12.89;
    double zmin = 0.901;
    double zmax = 1.801;
};

// ChannelCountXZAnalyzer writes per-frame molecule counts inside the ROI:
// x in [xmin, xmax) and z in [zmin, zmax), with PBC-safe interval checks.
class ChannelCountXZAnalyzer {
  public:
    explicit ChannelCountXZAnalyzer(const ChannelCountXZConfig& cfg = {});
    ChannelCountXZAnalyzer(const ChannelCountXZConfig& cfg, simio::runtime::CacheStore& cache);

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms,
                       int frame_idx);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }

  private:
    struct FrameRow {
        int frame_idx = 0;
        int64_t step = 0;
        double time_ps = 0.0;
        std::array<int64_t, 3> count{0, 0, 0};  // water, na, cl
    };

    ChannelCountXZConfig cfg_{};
    simio::runtime::CacheStore* cache_ = nullptr;
    simio::analysis::intrinsics::ChannelRoiX roi_{};
    bool has_roi_ = false;
    double xlen_ = 0.0;
    std::vector<double> xw_tmp_{};
    std::int64_t frame_counter_ = 0;
    int nframes_ = 0;
    std::vector<FrameRow> rows_{};
};

}  // namespace simio::analysis
