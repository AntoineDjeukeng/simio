#pragma once

#include <array>
#include <string>
#include <vector>

#include "common.hpp"

namespace simio::runtime {
class CacheStore;
}

namespace simio::analysis {

struct CoordXConfig {
    double xmin = 7.11;
    double xmax = 12.89;
    double zmin = 0.901;
    double zmax = 1.801;
    int nx = 100;
    double r_cw = 0.35;  // Na-O and Na<-W cutoff
    double r_aw = 0.38;  // Cl-O and Cl<-W cutoff
    double r_oo = 0.35;  // O-O cutoff
};

// CoordXAnalyzer computes coordination and HB observables vs x.
// It consumes cached molecule geometry from MolState and only performs per-frame reductions.
class CoordXAnalyzer {
  public:
    explicit CoordXAnalyzer(const CoordXConfig& cfg = {});
    CoordXAnalyzer(const CoordXConfig& cfg, simio::runtime::CacheStore& cache);

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }
    int nx() const { return cfg_.nx; }
    double x_center(int i) const { return x_centers_.at(static_cast<size_t>(i)); }
    double ibc_mean(int i) const { return stats_[IBC].mean(i); }
    double ibc_sem(int i) const { return stats_[IBC].sem(i); }
    double hbww_don_mean(int i) const { return stats_[HBWW_DON].mean(i); }
    double hbww_don_sem(int i) const { return stats_[HBWW_DON].sem(i); }
    double hbww_acc_mean(int i) const { return stats_[HBWW_ACC].mean(i); }
    double hbww_acc_sem(int i) const { return stats_[HBWW_ACC].sem(i); }
    double hbwcl_don_mean(int i) const { return stats_[HBWCL_DON].mean(i); }
    double hbwcl_don_sem(int i) const { return stats_[HBWCL_DON].sem(i); }

  private:
    enum Metric : int {
        IBC = 0,
        IBA = 1,
        BNC = 2,
        BNA = 3,
        BWW = 4,
        HBWW_DON = 5,
        HBWW_ACC = 6,
        HBWW_TOT = 7,
        HBWCL_DON = 8,
        MetricN = 9,
    };

    CoordXConfig cfg_{};
    simio::runtime::CacheStore* cache_ = nullptr;
    int nframes_ = 0;
    double dx_ = 0.0;
    bool has_x_centers_ = false;
    bool has_cached_rel_grid_ = false;
    std::vector<double> x_centers_rel_;
    std::vector<double> x_centers_;
    std::array<RunningStatsNonEmpty, MetricN> stats_{};
};

}  // namespace simio::analysis
