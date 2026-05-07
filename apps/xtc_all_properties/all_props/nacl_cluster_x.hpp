#pragma once

#include <array>
#include <string>
#include <vector>

#include "common.hpp"

namespace simio::runtime {
class CacheStore;
}

namespace simio::analysis {

struct NaClClusterXConfig {
    double zmin = 0.901;
    double zmax = 1.801;
    int nx = 100;
    double r_nacl = 0.59;
};

class NaClClusterXAnalyzer {
  public:
    explicit NaClClusterXAnalyzer(const NaClClusterXConfig& cfg = {});
    NaClClusterXAnalyzer(const NaClClusterXConfig& cfg, simio::runtime::CacheStore& cache);

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }
    int nx() const { return cfg_.nx; }

  private:
    struct IonNode {
        int mol_id = -1;
        Species species = Species::Na;
        Vec3d pos{};
    };

    NaClClusterXConfig cfg_{};
    simio::runtime::CacheStore* cache_ = nullptr;
    double dx_ = 0.0;
    std::vector<double> x_centers_rel_{};
    bool has_cached_rel_grid_ = false;
    int nframes_ = 0;
    double Lx_ref_ = 0.0;
    bool has_ref_box_ = false;

    RunningStatsAll cluster_count_{};
    RunningStatsAll ion_count_in_clusters_{};
    RunningStatsAll na_count_in_clusters_{};
    RunningStatsAll cl_count_in_clusters_{};
    RunningStatsNonEmpty cluster_size_{};
};

}  // namespace simio::analysis
