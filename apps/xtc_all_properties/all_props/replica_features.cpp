#include "replica_features.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace simio::analysis {
namespace {

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
};

struct SpeciesSpec {
    const char* label;
    const char* internal;
};

struct ClusterSpec {
    const char* label;
    const char* column;
};

const SpeciesSpec kSpecies[] = {
    {"SOL", "water"},
    {"NA", "na"},
    {"CL", "cl"},
};

const char* kRegions[] = {
    "center",
    "res_left",
    "res_right",
    "mouth_left",
    "mouth_right",
};

const ClusterSpec kClusterColumns[] = {
    {"cluster_count", "nacl_cluster_count_mean"},
    {"ions_in_clusters", "ions_in_nacl_clusters_mean"},
    {"na_in_clusters", "na_in_nacl_clusters_mean"},
    {"cl_in_clusters", "cl_in_nacl_clusters_mean"},
    {"cluster_size", "nacl_cluster_size_mean"},
};

const ClusterSpec kAssociationColumns[] = {
    {"CIP", "N_CIP_mean"},
    {"SSIP", "N_SSIP_mean"},
    {"ASSOC", "N_ASSOC_mean"},
    {"CN_NaOW", "CN_NaOW_mean"},
    {"CN_ClOW", "CN_ClOW_mean"},
    {"bridge_water", "N_bridge_water_mean"},
    {"bridged_pair", "N_bridged_pair_mean"},
    {"f_CIP", "f_CIP_mean"},
    {"f_SSIP", "f_SSIP_mean"},
    {"f_bridge", "f_bridge_mean"},
    {"CN_shared", "CN_shared_mean"},
    {"largest_ssip_component_size", "largest_ssip_component_size_mean"},
    {"largest_cip_component_size", "largest_cip_component_size_mean"},
    {"ssip_mean_degree", "ssip_mean_degree_mean"},
    {"cip_mean_degree", "cip_mean_degree_mean"},
};

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cell;
    std::istringstream iss(line);
    while (std::getline(iss, cell, ',')) out.push_back(cell);
    if (!line.empty() && line.back() == ',') out.push_back("");
    return out;
}

CsvTable read_csv(const std::string& path, bool required) {
    std::ifstream ifs(path);
    if (!ifs) {
        if (required) throw std::runtime_error("ReplicaFeatureExport: failed to open CSV: " + path);
        return {};
    }

    CsvTable t;
    std::string line;
    if (!std::getline(ifs, line)) return t;
    t.header = split_csv_line(line);
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        t.rows.push_back(split_csv_line(line));
    }
    return t;
}

int column_index(const CsvTable& t, const std::string& name, bool required) {
    for (size_t i = 0; i < t.header.size(); ++i) {
        if (t.header[i] == name) return static_cast<int>(i);
    }
    if (required) throw std::runtime_error("ReplicaFeatureExport: missing CSV column: " + name);
    return -1;
}

double as_double(const std::vector<std::string>& row, int idx) {
    if (idx < 0 || static_cast<size_t>(idx) >= row.size()) return std::numeric_limits<double>::quiet_NaN();
    try {
        return std::stod(row[static_cast<size_t>(idx)]);
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

std::string region_for_x(double x, double x_left, double x_right, double mouth_width_nm) {
    if (x >= x_left && x < x_right) return "center";
    if (x >= x_left - mouth_width_nm && x < x_left) return "mouth_left";
    if (x >= x_right && x < x_right + mouth_width_nm) return "mouth_right";
    if (x < x_left - mouth_width_nm) return "res_left";
    if (x >= x_right + mouth_width_nm) return "res_right";
    return "";
}

double mean(const std::vector<double>& xs) {
    double sum = 0.0;
    int n = 0;
    for (double x : xs) {
        if (std::isnan(x)) continue;
        sum += x;
        ++n;
    }
    if (n <= 0) return std::numeric_limits<double>::quiet_NaN();
    return sum / static_cast<double>(n);
}

double linear_slope(const std::vector<double>& xs, const std::vector<double>& ys) {
    if (xs.size() != ys.size() || xs.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    double sx = 0.0;
    double sy = 0.0;
    int n = 0;
    for (size_t i = 0; i < xs.size(); ++i) {
        if (std::isnan(xs[i]) || std::isnan(ys[i])) continue;
        sx += xs[i];
        sy += ys[i];
        ++n;
    }
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();

    const double xbar = sx / static_cast<double>(n);
    const double ybar = sy / static_cast<double>(n);
    double num = 0.0;
    double den = 0.0;
    for (size_t i = 0; i < xs.size(); ++i) {
        if (std::isnan(xs[i]) || std::isnan(ys[i])) continue;
        const double dx = xs[i] - xbar;
        num += dx * (ys[i] - ybar);
        den += dx * dx;
    }
    if (den == 0.0) return std::numeric_limits<double>::quiet_NaN();
    return num / den;
}

void add_feature(std::vector<std::pair<std::string, std::string>>& row,
                 const std::string& key,
                 const std::string& value) {
    row.push_back({key, value});
}

void add_feature(std::vector<std::pair<std::string, std::string>>& row,
                 const std::string& key,
                 double value) {
    std::ostringstream oss;
    oss << std::setprecision(12) << value;
    row.push_back({key, oss.str()});
}

std::vector<std::string> path_parts(const std::string& path) {
    std::vector<std::string> parts;
    std::string cur;
    for (char c : path) {
        if (c == '/' || c == '\\') {
            if (!cur.empty()) parts.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) parts.push_back(cur);
    return parts;
}

std::map<std::string, std::string> infer_metadata(const std::string& path) {
    std::map<std::string, std::string> meta;
    const std::regex h_re("^H[_-]?([0-9.]+)$", std::regex::icase);
    const std::regex l_re("^L[_-]?([0-9.]+)$", std::regex::icase);
    const std::regex field_re("^FIELD[_-]?([0-9]+)$", std::regex::icase);
    const std::regex rep_re("^rep[_-]?([0-9]+)$", std::regex::icase);

    for (const std::string& p : path_parts(path)) {
        std::smatch m;
        if (std::regex_match(p, m, h_re)) meta["H"] = m[1].str();
        if (std::regex_match(p, m, l_re)) meta["L"] = m[1].str();
        if (std::regex_match(p, m, field_re)) meta["field"] = m[1].str();
        if (std::regex_match(p, m, rep_re)) meta["replica"] = m[1].str();
        std::string lo = p;
        for (char& c : lo) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (lo == "neutral" || lo == "negative" || lo == "positive") meta["charge"] = lo;
    }
    return meta;
}

void append_density_features(std::vector<std::pair<std::string, std::string>>& out,
                             const std::string& density_csv,
                             double x_left,
                             double x_right,
                             double mouth_width_nm) {
    const CsvTable t = read_csv(density_csv, true);
    const int x_col = column_index(t, "x_center_nm", true);
    std::map<std::string, std::vector<double>> buckets;

    for (const auto& r : t.rows) {
        const std::string region = region_for_x(as_double(r, x_col), x_left, x_right, mouth_width_nm);
        if (region.empty()) continue;
        for (const SpeciesSpec& sp : kSpecies) {
            const std::string col = std::string("rho_") + sp.internal + "_mean";
            const int ci = column_index(t, col, false);
            if (ci < 0) continue;
            buckets[std::string(sp.label) + "_" + region + "_mean"].push_back(as_double(r, ci));
        }
    }

    for (const SpeciesSpec& sp : kSpecies) {
        for (const char* region : kRegions) {
            const std::string key = std::string(sp.label) + "_" + region + "_mean";
            add_feature(out, key, mean(buckets[key]));
        }
    }
}

void append_cluster_features(std::vector<std::pair<std::string, std::string>>& out,
                             const std::string& cluster_csv,
                             double x_left,
                             double x_right,
                             double mouth_width_nm) {
    const CsvTable t = read_csv(cluster_csv, false);
    if (t.header.empty()) return;
    const int x_col = column_index(t, "x_center_nm", true);
    std::map<std::string, std::vector<double>> buckets;

    for (const auto& r : t.rows) {
        const std::string region = region_for_x(as_double(r, x_col), x_left, x_right, mouth_width_nm);
        if (region.empty()) continue;
        for (const ClusterSpec& spec : kClusterColumns) {
            const int ci = column_index(t, spec.column, false);
            if (ci < 0) continue;
            buckets[std::string("NACL_") + spec.label + "_" + region + "_mean"].push_back(as_double(r, ci));
        }
    }

    for (const ClusterSpec& spec : kClusterColumns) {
        for (const char* region : kRegions) {
            const std::string key = std::string("NACL_") + spec.label + "_" + region + "_mean";
            add_feature(out, key, mean(buckets[key]));
        }
    }
}

void append_association_features(std::vector<std::pair<std::string, std::string>>& out,
                                 const std::string& association_csv,
                                 double x_left,
                                 double x_right,
                                 double mouth_width_nm) {
    const CsvTable t = read_csv(association_csv, false);
    if (t.header.empty()) return;
    const int x_col = column_index(t, "x_center_nm", true);
    std::map<std::string, std::vector<double>> buckets;

    for (const auto& r : t.rows) {
        const std::string region = region_for_x(as_double(r, x_col), x_left, x_right, mouth_width_nm);
        if (region.empty()) continue;
        for (const ClusterSpec& spec : kAssociationColumns) {
            const int ci = column_index(t, spec.column, false);
            if (ci < 0) continue;
            buckets[std::string("NACL_assoc_") + spec.label + "_" + region + "_mean"].push_back(as_double(r, ci));
        }
    }

    for (const ClusterSpec& spec : kAssociationColumns) {
        for (const char* region : kRegions) {
            const std::string key = std::string("NACL_assoc_") + spec.label + "_" + region + "_mean";
            add_feature(out, key, mean(buckets[key]));
        }
    }
}

void append_gating_features(std::vector<std::pair<std::string, std::string>>& out,
                            const std::string& gating_csv,
                            double slope_tail_ns) {
    const CsvTable t = read_csv(gating_csv, true);
    const int time_col = column_index(t, "time_ps", true);
    std::vector<double> time_ns;
    time_ns.reserve(t.rows.size());
    double t_end = 0.0;
    for (const auto& r : t.rows) {
        const double t_ns = as_double(r, time_col) / 1000.0;
        time_ns.push_back(t_ns);
        if (!std::isnan(t_ns)) t_end = std::max(t_end, t_ns);
    }
    const double t_start = t_end - slope_tail_ns;

    for (const SpeciesSpec& sp : kSpecies) {
        const std::string left_name = std::string("center_") + sp.internal + "_left";
        const std::string right_name = std::string("center_") + sp.internal + "_right";
        const std::string cum_right_name = std::string("center_") + sp.internal + "_cum_right";
        const int left_col = column_index(t, left_name, true);
        const int right_col = column_index(t, right_name, true);
        const int cum_right_col = column_index(t, cum_right_name, false);

        std::vector<double> left_cum;
        std::vector<double> right_cum;
        std::vector<double> net_cum;
        left_cum.reserve(t.rows.size());
        right_cum.reserve(t.rows.size());
        net_cum.reserve(t.rows.size());

        double left_sum = 0.0;
        double right_sum = 0.0;
        for (const auto& r : t.rows) {
            left_sum += as_double(r, left_col);
            if (cum_right_col >= 0) {
                right_sum = as_double(r, cum_right_col);
            } else {
                right_sum += as_double(r, right_col);
            }
            left_cum.push_back(left_sum);
            right_cum.push_back(right_sum);
            net_cum.push_back(right_sum - left_sum);
        }

        std::vector<double> fit_t;
        std::vector<double> fit_left;
        std::vector<double> fit_right;
        std::vector<double> fit_net;
        for (size_t i = 0; i < time_ns.size(); ++i) {
            if (std::isnan(time_ns[i]) || time_ns[i] < t_start) continue;
            fit_t.push_back(time_ns[i]);
            fit_left.push_back(left_cum[i]);
            fit_right.push_back(right_cum[i]);
            fit_net.push_back(net_cum[i]);
        }

        add_feature(out, std::string("center_") + sp.label + "_cum_left_final", left_cum.empty() ? 0.0 : left_cum.back());
        add_feature(out, std::string("center_") + sp.label + "_cum_right_final", right_cum.empty() ? 0.0 : right_cum.back());
        add_feature(out, std::string(sp.label) + "_net_slope", linear_slope(fit_t, fit_net));
        add_feature(out, std::string(sp.label) + "_left_slope", linear_slope(fit_t, fit_left));
        add_feature(out, std::string(sp.label) + "_right_slope", linear_slope(fit_t, fit_right));
    }
}

}  // namespace

void write_replica_features_csv(const std::string& out_csv,
                                const ReplicaFeatureExportConfig& cfg,
                                const std::string& density_x_csv,
                                const std::string& gating_flux_csv,
                                const std::string& nacl_cluster_x_csv,
                                const std::string& nacl_association_x_csv) {
    if (cfg.carbon_right_x_nm <= cfg.carbon_left_x_nm) {
        throw std::runtime_error("ReplicaFeatureExport: invalid carbon/channel x bounds");
    }
    if (cfg.mouth_width_nm <= 0.0 || cfg.slope_tail_ns <= 0.0) {
        throw std::runtime_error("ReplicaFeatureExport: mouth_width_nm and slope_tail_ns must be > 0");
    }

    std::vector<std::pair<std::string, std::string>> row;
    const auto meta = infer_metadata(cfg.metadata_path);
    for (const char* key : {"H", "L", "charge", "field", "replica"}) {
        const auto it = meta.find(key);
        add_feature(row, key, it == meta.end() ? "" : it->second);
    }
    add_feature(row, "carbon_left_x_nm", cfg.carbon_left_x_nm);
    add_feature(row, "carbon_right_x_nm", cfg.carbon_right_x_nm);
    add_feature(row, "mouth_width_nm", cfg.mouth_width_nm);
    add_feature(row, "slope_tail_ns", cfg.slope_tail_ns);

    append_density_features(row, density_x_csv, cfg.carbon_left_x_nm, cfg.carbon_right_x_nm, cfg.mouth_width_nm);
    append_gating_features(row, gating_flux_csv, cfg.slope_tail_ns);
    append_cluster_features(row, nacl_cluster_x_csv, cfg.carbon_left_x_nm, cfg.carbon_right_x_nm, cfg.mouth_width_nm);
    append_association_features(row, nacl_association_x_csv, cfg.carbon_left_x_nm, cfg.carbon_right_x_nm,
                                cfg.mouth_width_nm);

    std::ofstream ofs(out_csv);
    if (!ofs) throw std::runtime_error("ReplicaFeatureExport: failed to open output CSV: " + out_csv);
    for (size_t i = 0; i < row.size(); ++i) {
        if (i) ofs << ",";
        ofs << row[i].first;
    }
    ofs << "\n";
    for (size_t i = 0; i < row.size(); ++i) {
        if (i) ofs << ",";
        ofs << row[i].second;
    }
    ofs << "\n";
}

}  // namespace simio::analysis
