#pragma once

#include <string>

#include "common.hpp"

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

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }

  private:
    DipoleZInXChannelConfig cfg_{};
    int nframes_ = 0;
    double Lz_ref_ = 0.0;
    bool has_ref_box_ = false;

    RunningStatsAll mux_{};
    RunningStatsAll mux_fold_{};
    RunningStatsAll muz_{};
    RunningStatsAll cnt_{};
};

}  // namespace simio::analysis
