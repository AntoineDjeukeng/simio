#pragma once

#include <array>
#include <string>
#include <vector>

#include "common.hpp"

namespace simio::runtime {
class CacheStore;
}

namespace simio::analysis {

struct WaterAtomDensityXConfig {
    double zmin = 0.901;
    double zmax = 1.801;
    int nx = 100;
};

// WaterAtomDensityXAnalyzer computes atom-resolved water density rho(x) in a z slab.
// O_w uses the oxygen position; H_w counts both hydrogen atoms explicitly.
class WaterAtomDensityXAnalyzer {
  public:
    explicit WaterAtomDensityXAnalyzer(const WaterAtomDensityXConfig& cfg = {});
    WaterAtomDensityXAnalyzer(const WaterAtomDensityXConfig& cfg, simio::runtime::CacheStore& cache);

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms);
    void write_csv(const std::string& path) const;

    int nframes() const { return nframes_; }
    int nx() const { return cfg_.nx; }
    double x_center(int i) const { return (static_cast<double>(i) + 0.5) * dx_; }
    double rho_ow_mean(int i) const { return rho_[Ow].mean(i, nframes_); }
    double rho_ow_sem(int i) const { return rho_[Ow].sem(i, nframes_); }
    double rho_hw_mean(int i) const { return rho_[Hw].mean(i, nframes_); }
    double rho_hw_sem(int i) const { return rho_[Hw].sem(i, nframes_); }

  private:
    enum Site : int {
        Ow = 0,
        Hw = 1,
        SiteN = 2,
    };

    WaterAtomDensityXConfig cfg_{};
    simio::runtime::CacheStore* cache_ = nullptr;
    double dx_ = 0.0;
    std::vector<double> x_centers_rel_{};
    bool has_cached_rel_grid_ = false;
    int nframes_ = 0;
    double Lx_ref_ = 0.0;
    bool has_ref_box_ = false;

    std::array<RunningStatsAll, SiteN> rho_{};
    std::array<RunningStatsAll, SiteN> cnt_{};
};

}  // namespace simio::analysis
