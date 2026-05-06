#pragma once

#include <string>

namespace simio::analysis {

struct ReplicaFeatureExportConfig {
    double carbon_left_x_nm = 0.0;
    double carbon_right_x_nm = 0.0;
    double mouth_width_nm = 0.62;
    double slope_tail_ns = 45.0;
    std::string metadata_path;
};

void write_replica_features_csv(const std::string& out_csv,
                                const ReplicaFeatureExportConfig& cfg,
                                const std::string& density_x_csv,
                                const std::string& gating_flux_csv,
                                const std::string& nacl_cluster_x_csv,
                                const std::string& nacl_association_x_csv);

}  // namespace simio::analysis
