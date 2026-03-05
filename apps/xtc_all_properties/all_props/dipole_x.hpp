#pragma once

#include <string>
#include <vector>

#include "common.hpp"

namespace simio::runtime {
class CacheStore;
}

namespace simio::analysis {

struct DipoleXConfig {
    double xmin = 7.11;
    double xmax = 12.89;
    double zmin = 0.901;
    double zmax = 1.801;
    int nx = 100;
};

// DipoleXAnalyzer computes mean water orientation components versus x.
// It reuses cached per-water geometry from MolState when available.
class DipoleXAnalyzer {
  public:
    explicit DipoleXAnalyzer(const DipoleXConfig& cfg = {});
    DipoleXAnalyzer(const DipoleXConfig& cfg, simio::runtime::CacheStore& cache);

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }

  private:
    DipoleXConfig cfg_{};
    simio::runtime::CacheStore* cache_ = nullptr;
    int nframes_ = 0;
    double dx_ = 0.0;
    bool has_x_centers_ = false;
    bool has_cached_rel_grid_ = false;
    std::vector<double> x_centers_rel_;
    std::vector<double> x_centers_;

    RunningStatsAll muz_raw_{};
    RunningStatsAll muz_fold_{};
    RunningStatsAll mux_{};
    RunningStatsAll cnt_{};
};

}  // namespace simio::analysis
