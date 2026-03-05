#pragma once

#include <array>
#include <string>

#include "common.hpp"

namespace simio::analysis {

struct DensityZInXChannelConfig {
    double xmin = 7.11;
    double xmax = 12.89;
    int nz = 100;
};

// DensityZInXChannelAnalyzer computes species-resolved number density rho(z)
// in an x-channel slab [xmin, xmax) (with PBC wrapping on x).
class DensityZInXChannelAnalyzer {
  public:
    explicit DensityZInXChannelAnalyzer(const DensityZInXChannelConfig& cfg = {});

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }

  private:
    DensityZInXChannelConfig cfg_{};
    int nframes_ = 0;
    double Lz_ref_ = 0.0;
    bool has_ref_box_ = false;

    std::array<RunningStatsAll, 3> rho_{};
    RunningStatsAll rho_q_{};
    std::array<RunningStatsAll, 3> cnt_{};
};

}  // namespace simio::analysis
