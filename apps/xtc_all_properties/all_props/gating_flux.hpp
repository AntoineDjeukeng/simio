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

enum class GatingSelection : uint8_t {
    All = 0,
    InChannelNow = 1,
    InChannelBoth = 2,
};

GatingSelection parse_gating_selection(const std::string& v);

struct GatingFluxConfig {
    double xmin = 7.11;
    double xmax = 12.89;
    double zmin = 0.901;
    double zmax = 1.801;
    GatingSelection selection = GatingSelection::All;
};

// GatingFluxAnalyzer detects crossings through center and seam planes
// using minimum-image segment crossings with z-slab filtering.
class GatingFluxAnalyzer {
  public:
    explicit GatingFluxAnalyzer(const GatingFluxConfig& cfg = {});
    GatingFluxAnalyzer(const GatingFluxConfig& cfg, simio::runtime::CacheStore& cache);

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms,
                       int frame_idx);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }

  private:
    struct PlaneTally {
        int64_t n_left = 0;
        int64_t n_right = 0;
        std::array<int64_t, 3> n_left_species{0, 0, 0};
        std::array<int64_t, 3> n_right_species{0, 0, 0};
    };

    struct FrameRow {
        int frame_idx = 0;
        int64_t step = 0;
        double time_ps = 0.0;
        PlaneTally center_frame{};
        PlaneTally center_cum{};
        PlaneTally seam_frame{};
        PlaneTally seam_cum{};
    };

    GatingFluxConfig cfg_{};
    simio::runtime::CacheStore* cache_ = nullptr;
    simio::analysis::intrinsics::ChannelRoiX roi_{};
    bool has_roi_ = false;
    double xlen_ = 0.0;
    std::vector<double> xw_tmp_{};
    std::int64_t frame_counter_ = 0;
    int nframes_ = 0;

    std::vector<Vec3d> prev_key_;
    std::vector<uint8_t> has_prev_;
    std::vector<uint8_t> prev_in_channel_;

    PlaneTally center_cum_{};
    PlaneTally seam_cum_{};
    std::vector<FrameRow> rows_;
};

}  // namespace simio::analysis
