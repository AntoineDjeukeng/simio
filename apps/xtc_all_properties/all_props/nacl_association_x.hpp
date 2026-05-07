#pragma once

#include <array>
#include <string>
#include <vector>

#include "common.hpp"

namespace simio::runtime {
class CacheStore;
}

namespace simio::analysis {

struct NaClAssociationXConfig {
    double zmin = 0.901;
    double zmax = 1.801;
    int nx = 100;
    double r_cip = 0.35;
    double r_ssip = 0.55;
    double r_naow = 0.35;
    double r_clow = 0.38;
};

class NaClAssociationXAnalyzer {
  public:
    explicit NaClAssociationXAnalyzer(const NaClAssociationXConfig& cfg = {});
    NaClAssociationXAnalyzer(const NaClAssociationXConfig& cfg, simio::runtime::CacheStore& cache);

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }

  private:
    NaClAssociationXConfig cfg_{};
    simio::runtime::CacheStore* cache_ = nullptr;
    double dx_ = 0.0;
    bool has_cached_rel_grid_ = false;
    bool has_ref_box_ = false;
    int nframes_ = 0;
    std::vector<double> x_centers_rel_;

    RunningStatsAll cip_{};
    RunningStatsAll ssip_{};
    RunningStatsAll assoc_{};
    RunningStatsAll bridge_water_{};
    RunningStatsAll bridged_pair_{};
    RunningStatsAll f_cip_{};
    RunningStatsAll f_ssip_{};
    RunningStatsAll f_bridge_{};
    RunningStatsAll cn_shared_{};
    RunningStatsAll largest_ssip_component_size_{};
    RunningStatsAll largest_cip_component_size_{};
    RunningStatsNonEmpty ssip_mean_degree_{};
    RunningStatsNonEmpty cip_mean_degree_{};
    RunningStatsNonEmpty cn_naow_{};
    RunningStatsNonEmpty cn_clow_{};
};

}  // namespace simio::analysis
