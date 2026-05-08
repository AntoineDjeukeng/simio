#include "axial_profile_x_export.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace simio::analysis {

namespace {

double wrap_pos_1d(double x, double L) {
    if (L <= 0.0) return x;
    const double y = std::fmod(x, L);
    return (y < 0.0) ? (y + L) : y;
}

double min_image_1d(double dx, double L) {
    if (L <= 0.0) return dx;
    return dx - L * std::round(dx / L);
}

double interval_midpoint_1d(double xminw, double xmaxw, double L) {
    if (L <= 0.0) return 0.0;
    if (xminw <= xmaxw) return xminw + 0.5 * (xmaxw - xminw);
    return wrap_pos_1d(xminw + 0.5 * ((xmaxw + L) - xminw), L);
}

struct AxialProfileRow {
    double x_center = 0.0;
    double x_relative = 0.0;
    int in_channel = 0;
    double rho_ow_mean = 0.0;
    double rho_ow_sem = 0.0;
    double rho_hw_mean = 0.0;
    double rho_hw_sem = 0.0;
    double rho_na_mean = 0.0;
    double rho_na_sem = 0.0;
    double rho_cl_mean = 0.0;
    double rho_cl_sem = 0.0;
    double mux_mean = 0.0;
    double mux_sem = 0.0;
    double muz_mean = 0.0;
    double muz_sem = 0.0;
    double muz_fold_mean = 0.0;
    double muz_fold_sem = 0.0;
    double ac_mean = 0.0;
    double ac_sem = 0.0;
    double na_bond_mean = 0.0;
    double na_bond_sem = 0.0;
    double dn_mean = 0.0;
    double dn_sem = 0.0;
    double cl_bond_mean = 0.0;
    double cl_bond_sem = 0.0;
};

void require_same_x(double ref, double got, double tol, const char* label, int i) {
    if (std::abs(ref - got) > tol) {
        std::ostringstream oss;
        oss << std::setprecision(17)
            << "axial_profile_x: x grid mismatch for " << label
            << " at bin " << i
            << " ref=" << ref
            << " got=" << got
            << " diff=" << (got - ref)
            << " tol=" << tol;
        throw std::runtime_error(oss.str());
    }
}

}  // namespace

void write_axial_profile_x_csv(const std::string& path,
                               const AxialProfileXExportConfig& cfg,
                               const DensityXAnalyzer& density,
                               const WaterAtomDensityXAnalyzer& water_atom_density,
                               const DipoleXAnalyzer& dipole,
                               const CoordXAnalyzer& coord) {
    const int nx = density.nx();
    if (nx <= 0) throw std::runtime_error("axial_profile_x: invalid nx");
    if (water_atom_density.nx() != nx || dipole.nx() != nx || coord.nx() != nx) {
        throw std::runtime_error("axial_profile_x: analyzer nx mismatch");
    }
    if (cfg.box_length_x_nm <= 0.0) {
        throw std::runtime_error("axial_profile_x: invalid box_length_x_nm");
    }

    const double xminw = wrap_pos_1d(cfg.xmin, cfg.box_length_x_nm);
    const double xmaxw = wrap_pos_1d(cfg.xmax, cfg.box_length_x_nm);
    const double x_center_channel = interval_midpoint_1d(xminw, xmaxw, cfg.box_length_x_nm);
    const double x_edge_left_rel = min_image_1d(xminw - x_center_channel, cfg.box_length_x_nm);
    const double x_edge_right_rel = min_image_1d(xmaxw - x_center_channel, cfg.box_length_x_nm);
    const double xtol = 1e-9 * std::max(1.0, cfg.box_length_x_nm);

    std::vector<AxialProfileRow> rows;
    rows.reserve(static_cast<size_t>(nx));
    for (int i = 0; i < nx; ++i) {
        const double x_center = density.x_center(i);
        require_same_x(x_center, water_atom_density.x_center(i), xtol, "water_atom_density", i);
        require_same_x(x_center, dipole.x_center(i), xtol, "dipole", i);
        require_same_x(x_center, coord.x_center(i), xtol, "coord", i);

        AxialProfileRow row{};
        row.x_center = x_center;
        row.x_relative = min_image_1d(x_center - x_center_channel, cfg.box_length_x_nm);
        row.in_channel = in_range_pbc(x_center, xminw, xmaxw) ? 1 : 0;
        row.rho_ow_mean = water_atom_density.rho_ow_mean(i);
        row.rho_ow_sem = water_atom_density.rho_ow_sem(i);
        row.rho_hw_mean = water_atom_density.rho_hw_mean(i);
        row.rho_hw_sem = water_atom_density.rho_hw_sem(i);
        row.rho_na_mean = density.rho_mean(Species::Na, i);
        row.rho_na_sem = density.rho_sem(Species::Na, i);
        row.rho_cl_mean = density.rho_mean(Species::Cl, i);
        row.rho_cl_sem = density.rho_sem(Species::Cl, i);
        row.mux_mean = dipole.mux_mean(i);
        row.mux_sem = dipole.mux_sem(i);
        row.muz_mean = dipole.muz_mean(i);
        row.muz_sem = dipole.muz_sem(i);
        row.muz_fold_mean = dipole.muz_fold_mean(i);
        row.muz_fold_sem = dipole.muz_fold_sem(i);
        row.ac_mean = coord.hbww_acc_mean(i);
        row.ac_sem = coord.hbww_acc_sem(i);
        row.na_bond_mean = coord.ibc_mean(i);
        row.na_bond_sem = coord.ibc_sem(i);
        row.dn_mean = coord.hbww_don_mean(i);
        row.dn_sem = coord.hbww_don_sem(i);
        row.cl_bond_mean = coord.hbwcl_don_mean(i);
        row.cl_bond_sem = coord.hbwcl_don_sem(i);
        rows.push_back(row);
    }

    std::sort(rows.begin(), rows.end(), [](const AxialProfileRow& a, const AxialProfileRow& b) {
        return a.x_relative < b.x_relative;
    });

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("axial_profile_x: failed to open output CSV: " + path);
    ofs << std::setprecision(12);
    ofs << "x_center_nm,x_relative_nm,in_channel,"
           "channel_center_nm,channel_edge_left_nm,channel_edge_right_nm,"
           "channel_edge_left_relative_nm,channel_edge_right_relative_nm,"
           "rho_ow_mean,rho_ow_sem,rho_hw_mean,rho_hw_sem,rho_na_mean,rho_na_sem,rho_cl_mean,rho_cl_sem,"
           "mux_mean,mux_sem,muz_mean,muz_sem,muz_fold_mean,muz_fold_sem,"
           "ac_mean,ac_sem,na_bond_mean,na_bond_sem,ac_plus_na_mean,"
           "dn_mean,dn_sem,cl_bond_mean,cl_bond_sem,dn_plus_cl_mean\n";

    for (const AxialProfileRow& row : rows) {
        ofs << row.x_center << "," << row.x_relative << "," << row.in_channel << ","
            << x_center_channel << "," << xminw << "," << xmaxw << "," << x_edge_left_rel << ","
            << x_edge_right_rel << ","
            << row.rho_ow_mean << "," << row.rho_ow_sem << "," << row.rho_hw_mean << ","
            << row.rho_hw_sem << "," << row.rho_na_mean << "," << row.rho_na_sem << ","
            << row.rho_cl_mean << "," << row.rho_cl_sem << ","
            << row.mux_mean << "," << row.mux_sem << "," << row.muz_mean << "," << row.muz_sem
            << "," << row.muz_fold_mean << "," << row.muz_fold_sem << ","
            << row.ac_mean << "," << row.ac_sem << "," << row.na_bond_mean << ","
            << row.na_bond_sem << "," << (row.ac_mean + row.na_bond_mean) << ","
            << row.dn_mean << "," << row.dn_sem << "," << row.cl_bond_mean << ","
            << row.cl_bond_sem << "," << (row.dn_mean + row.cl_bond_mean) << "\n";
    }
}

}  // namespace simio::analysis
