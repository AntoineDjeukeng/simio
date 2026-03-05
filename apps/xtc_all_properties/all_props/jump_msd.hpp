#pragma once

#include <array>
#include <cstdint>
#include <deque>
#include <string>
#include <vector>

#include "common.hpp"

namespace simio::analysis {

struct JumpMsdConfig {
    double zmin = 0.901;
    double zmax = 1.801;
    double x_channel_min = 7.11;
    double x_channel_max = 12.89;
    double bound_layer_nm = 0.30;
    int keep_frames = 50;  // Maximum lag in frames (rows 1..keep_frames).
};

// JumpMsdAnalyzer computes lag-resolved sliding-window statistics from continuous
// trajectories (MolState::time.key_cont), with full x/y and z-slab restriction.
// Max-jump is reported for x only using min-image displacement between two frames.
class JumpMsdAnalyzer {
  public:
    explicit JumpMsdAnalyzer(const JumpMsdConfig& cfg = {});

    void process_frame(const Topology& topo, const Frame& fr, const std::vector<MolState>& ms,
                       int frame_idx);
    void write_csv(const std::string& path) const;
    void write_channel_msd_x_csv(const std::string& path) const;
    void write_channel_msd_z_csv(const std::string& path) const;
    void write_state_channel_csv(const std::string& path) const;
    void write_drift_channel_csv(const std::string& path) const;
    void write_vacf_y_csv(const std::string& path,
                          std::array<double, 3>* plateau_last10_nm2_per_ps = nullptr) const;
    void write_vacf_x_channel_csv(
        const std::string& path,
        std::array<double, 3>* plateau_last10_nm2_per_ps = nullptr) const;
    void write_vacf_z_channel_raw_csv(
        const std::string& path,
        std::array<double, 3>* plateau_last10_nm2_per_ps = nullptr) const;

    int nframes() const { return nframes_; }

  private:
    struct Snapshot {
        int frame_idx = -1;
        int64_t step = 0;
        double time_ps = 0.0;
        std::vector<Vec3d> key_cont;
        std::vector<double> x_raw;
        std::vector<double> z_raw;
        std::vector<double> vx_raw;
        std::vector<double> vy;
        std::vector<double> vz_raw;
        std::vector<uint8_t> valid;
        std::vector<uint8_t> vel_valid;
        std::vector<uint8_t> in_slab;
        std::vector<uint8_t> in_channel;
        std::vector<uint8_t> in_bound;
        std::vector<uint8_t> in_core;
    };

    struct StateFrameRow {
        int frame_idx = -1;
        int64_t step = 0;
        double time_ps = 0.0;
        std::array<int32_t, 3> n_bound{0, 0, 0};
        std::array<int32_t, 3> n_core{0, 0, 0};
    };

    struct DriftAccum {
        double sum = 0.0;
        double sumsq = 0.0;
        int64_t n = 0;
        void add(double v) {
            sum += v;
            sumsq += v * v;
            n += 1;
        }
    };

    JumpMsdConfig cfg_{};
    int nframes_ = 0;
    size_t nmol_ = 0;
    std::vector<int> species_idx_;
    std::deque<Snapshot> history_;

    std::vector<int64_t> n_samples_by_lag_;
    RunningStatsNonEmpty lag_time_ps_;

    RunningStatsNonEmpty msd_x_all_;
    RunningStatsNonEmpty msd_y_all_;
    RunningStatsNonEmpty msd_z_all_;
    RunningStatsNonEmpty msd_z_raw_all_;
    RunningStatsNonEmpty msd_r2_all_;

    std::array<RunningStatsNonEmpty, 3> msd_x_sp_{};
    std::array<RunningStatsNonEmpty, 3> msd_y_sp_{};
    std::array<RunningStatsNonEmpty, 3> msd_z_sp_{};
    std::array<RunningStatsNonEmpty, 3> msd_z_raw_sp_{};
    std::array<RunningStatsNonEmpty, 3> msd_r2_sp_{};

    RunningStatsNonEmpty max_jump_x_all_;
    std::array<RunningStatsNonEmpty, 3> max_jump_x_sp_{};

    std::vector<int32_t> runlen_in_channel_;
    std::array<std::vector<double>, 3> msd_x_channel_strict_sum_dx2_{};
    std::array<std::vector<int64_t>, 3> msd_x_channel_strict_n_{};
    std::array<std::vector<double>, 3> msd_x_channel_loose_sum_dx2_{};
    std::array<std::vector<int64_t>, 3> msd_x_channel_loose_n_{};
    std::array<std::vector<double>, 3> msd_z_channel_unwrapped_sum_dz2_{};
    std::array<std::vector<double>, 3> msd_z_channel_raw_sum_dz2_{};
    std::array<std::vector<int64_t>, 3> msd_z_channel_n_{};
    std::array<double, 3> max_abs_dz_lag1_sp_{0.0, 0.0, 0.0};
    std::vector<double> msd_x_channel_strict_lag_ps_sum_{};
    std::vector<int64_t> msd_x_channel_strict_lag_ps_n_{};
    std::vector<double> msd_x_channel_loose_lag_ps_sum_{};
    std::vector<int64_t> msd_x_channel_loose_lag_ps_n_{};

    std::vector<StateFrameRow> state_rows_{};
    std::array<DriftAccum, 3> drift_vx_core_{};
    std::array<DriftAccum, 3> drift_vx_bound_{};
    double drift_dt_ps_sum_ = 0.0;
    int64_t drift_dt_ps_n_ = 0;

    std::array<std::vector<double>, 3> vacf_y_sum_{};
    std::array<std::vector<int64_t>, 3> vacf_y_n_{};
    std::vector<double> vacf_y_lag_ps_sum_{};
    std::vector<int64_t> vacf_y_lag_ps_n_{};

    std::array<std::vector<double>, 3> vacf_x_channel_sum_{};
    std::array<std::vector<int64_t>, 3> vacf_x_channel_n_{};
    std::vector<double> vacf_x_channel_lag_ps_sum_{};
    std::vector<int64_t> vacf_x_channel_lag_ps_n_{};

    std::array<std::vector<double>, 3> vacf_z_channel_raw_sum_{};
    std::array<std::vector<int64_t>, 3> vacf_z_channel_raw_n_{};
    std::vector<double> vacf_z_channel_raw_lag_ps_sum_{};
    std::vector<int64_t> vacf_z_channel_raw_lag_ps_n_{};
};

}  // namespace simio::analysis
