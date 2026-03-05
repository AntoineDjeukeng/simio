#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "common.hpp"
#include "simio/analysis/intrinsics/channel_roi.hpp"

namespace simio::runtime {
class CacheStore;
}

namespace simio::analysis {

struct DipoleZInXChannelConfig {
    double xmin = 7.11;
    double xmax = 12.89;
    int nz = 100;
};

// DipoleZInXChannelAnalyzer computes mean water orientation components
// (mu_x, mu_z) versus z in an x-channel slab [xmin, xmax) (with PBC on x).
class DipoleZInXChannelAnalyzer {
  public:
    explicit DipoleZInXChannelAnalyzer(const DipoleZInXChannelConfig& cfg = {});
    DipoleZInXChannelAnalyzer(const DipoleZInXChannelConfig& cfg, simio::runtime::CacheStore& cache);

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }

  private:
    DipoleZInXChannelConfig cfg_{};
    simio::runtime::CacheStore* cache_ = nullptr;
    simio::analysis::intrinsics::ChannelRoiX roi_{};
    bool has_roi_ = false;
    double xlen_ = 0.0;
    double dz_ = 0.0;
    bool has_cached_rel_grid_ = false;
    int nframes_ = 0;
    double Lz_ref_ = 0.0;
    bool has_ref_box_ = false;
    std::vector<double> xw_tmp_{};
    std::int64_t frame_counter_ = 0;

    RunningStatsAll mux_{};
    RunningStatsAll mux_fold_{};
    RunningStatsAll muz_{};
    RunningStatsAll cnt_{};
};

}  // namespace simio::analysis
