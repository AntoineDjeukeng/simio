#pragma once

#include <array>
#include <string>
#include <vector>

#include "common.hpp"

namespace simio::runtime {
class CacheStore;
}

namespace simio::analysis {

struct DensityXConfig {
    double zmin = 0.901;
    double zmax = 1.801;
    int nx = 100;
};

// DensityXAnalyzer computes species-resolved number density rho(x) in a z slab.
// It assumes per-molecule key positions are already prepared in MolState cache.
class DensityXAnalyzer {
  public:
    explicit DensityXAnalyzer(const DensityXConfig& cfg = {});
    DensityXAnalyzer(const DensityXConfig& cfg, simio::runtime::CacheStore& cache);

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }

  private:
    DensityXConfig cfg_{};
    simio::runtime::CacheStore* cache_ = nullptr;
    double dx_ = 0.0;
    std::vector<double> x_centers_rel_{};
    bool has_cached_rel_grid_ = false;
    int nframes_ = 0;
    double Lx_ref_ = 0.0;
    bool has_ref_box_ = false;

    std::array<RunningStatsAll, 3> rho_{};
    std::array<RunningStatsAll, 3> cnt_{};
};

}  // namespace simio::analysis
