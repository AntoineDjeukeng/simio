#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "simio/layered_pipeline.hpp"
#include "simio/properties/GatingCenterPlane.hpp"
#include "simio/properties/GatingPbcWrap.hpp"
#include "simio/simio.hpp"
#include "xtc.h"

namespace {

struct CliConfig {
    const char* xtc_path = nullptr;
    int max_frames = 100;
    int nthreads = 4;
    int nsol = 5896;
    int nna = 110;
    int ncl = 110;
    int nx = 100;
    int nz = 100;
    double x_center_wrapped = std::numeric_limits<double>::quiet_NaN();
    double z_center_wrapped = std::numeric_limits<double>::quiet_NaN();
    double channel_length_x = 0.0;
    double channel_height_z = 0.0;
    simio::properties::MoleculeSelection selection = simio::properties::MoleculeSelection::All;
    std::string out_csv = "gating_flux.csv";
};

int parse_int_arg(const char* s, const char* name) {
    try {
        return std::stoi(std::string(s));
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Invalid integer for ") + name + ": " + s);
    }
}

double parse_double_arg(const char* s, const char* name) {
    try {
        return std::stod(std::string(s));
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Invalid number for ") + name + ": " + s);
    }
}

simio::properties::MoleculeSelection parse_selection_arg(const char* s) {
    std::string v = s;
    std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return (char)std::tolower(c); });

    if (v == "0" || v == "all") return simio::properties::MoleculeSelection::All;
    if (v == "1" || v == "now" || v == "inchannelnow") {
        return simio::properties::MoleculeSelection::InChannelNow;
    }
    if (v == "2" || v == "both" || v == "inchannelboth") {
        return simio::properties::MoleculeSelection::InChannelBoth;
    }
    throw std::runtime_error("Invalid selection (use all|now|both or 0|1|2): " + std::string(s));
}

void print_usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <trajectory.xtc> [max_frames=100] [threads=4] [nsol=5896] [nna=110] [ncl=110] [nx=100]"
                 " [nz=100] [x_center_wrapped=nan] [z_center_wrapped=nan] [channel_length_x=0.0]"
                 " [channel_height_z=0.0] [selection=all] [out_csv=gating_flux.csv]\n";
}

CliConfig parse_cli(int argc, char** argv) {
    if (argc < 2) throw std::runtime_error("Missing trajectory path.");

    CliConfig c;
    c.xtc_path = argv[1];
    if (argc > 2) c.max_frames = parse_int_arg(argv[2], "max_frames");
    if (argc > 3) c.nthreads = parse_int_arg(argv[3], "threads");
    if (argc > 4) c.nsol = parse_int_arg(argv[4], "nsol");
    if (argc > 5) c.nna = parse_int_arg(argv[5], "nna");
    if (argc > 6) c.ncl = parse_int_arg(argv[6], "ncl");
    if (argc > 7) c.nx = parse_int_arg(argv[7], "nx");
    if (argc > 8) c.nz = parse_int_arg(argv[8], "nz");
    if (argc > 9) c.x_center_wrapped = parse_double_arg(argv[9], "x_center_wrapped");
    if (argc > 10) c.z_center_wrapped = parse_double_arg(argv[10], "z_center_wrapped");
    if (argc > 11) c.channel_length_x = parse_double_arg(argv[11], "channel_length_x");
    if (argc > 12) c.channel_height_z = parse_double_arg(argv[12], "channel_height_z");
    if (argc > 13) c.selection = parse_selection_arg(argv[13]);
    if (argc > 14) c.out_csv = argv[14];

    if (c.max_frames <= 0 || c.nthreads <= 0 || c.nsol < 0 || c.nna < 0 || c.ncl < 0 || c.nx <= 0 ||
        c.nz <= 0 || c.channel_length_x < 0.0 || c.channel_height_z < 0.0) {
        throw std::runtime_error("Invalid non-positive numeric argument.");
    }
    return c;
}

simio::Topology build_topology(int nsol, int nna, int ncl) {
    simio::Topology topo;
    int atom_cursor = 0;

    for (int i = 0; i < nsol; ++i) {
        topo.mols.push_back(simio::MolSpan{atom_cursor, 3, simio::MolType::Water});
        atom_cursor += 3;
    }
    for (int i = 0; i < nna; ++i) {
        topo.mols.push_back(simio::MolSpan{atom_cursor, 1, simio::MolType::Cation});
        atom_cursor += 1;
    }
    for (int i = 0; i < ncl; ++i) {
        topo.mols.push_back(simio::MolSpan{atom_cursor, 1, simio::MolType::Anion});
        atom_cursor += 1;
    }

    topo.build_type_lists();
    return topo;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            print_usage(argv[0]);
            return 1;
        }

        const CliConfig cfg = parse_cli(argc, argv);
        const int expected_natoms = cfg.nsol * 3 + cfg.nna + cfg.ncl;
        const int expected_nmols = cfg.nsol + cfg.nna + cfg.ncl;

        simio::Topology topo = build_topology(cfg.nsol, cfg.nna, cfg.ncl);
        if ((int)topo.mols.size() != expected_nmols) {
            throw std::runtime_error("Topology molecule count mismatch.");
        }

        XtcTraj traj{};
        int rc = xtc_open(&traj, cfg.xtc_path, /*cap=*/1);
        if (rc != exdrOK) {
            std::cerr << "Failed to open XTC: rc=" << rc << "\n";
            return 1;
        }
        if (traj.natoms != expected_natoms) {
            std::cerr << "Natoms mismatch: file=" << traj.natoms << " expected=" << expected_natoms << "\n";
            xtc_close(&traj);
            return 1;
        }

        simio::Frame fr;
        fr.atoms.x.resize((size_t)expected_natoms);
        fr.atoms.y.resize((size_t)expected_natoms);
        fr.atoms.z.resize((size_t)expected_natoms);

        simio::layered::Layer0Config l0;
        l0.nx = cfg.nx;
        l0.nz = cfg.nz;
        l0.x_center_wrapped = cfg.x_center_wrapped;
        l0.z_center_wrapped = cfg.z_center_wrapped;
        l0.channel_length_x = cfg.channel_length_x;
        l0.channel_height_z = cfg.channel_height_z;
        simio::layered::Layer1Config l1;
        simio::layered::Pipeline pipe(l0, l1);

        simio::properties::GatingCenterPlaneConfig center_cfg;
        center_cfg.selection = cfg.selection;
        auto center = std::make_unique<simio::properties::GatingCenterPlaneProperty>(center_cfg);
        simio::properties::GatingCenterPlaneProperty* center_ptr = center.get();
        pipe.add_property(std::move(center));

        simio::properties::GatingPbcWrapConfig wrap_cfg;
        wrap_cfg.selection = cfg.selection;
        auto wrap = std::make_unique<simio::properties::GatingPbcWrapProperty>(wrap_cfg);
        simio::properties::GatingPbcWrapProperty* wrap_ptr = wrap.get();
        pipe.add_property(std::move(wrap));

        std::ofstream ofs(cfg.out_csv);
        if (!ofs) {
            xtc_close(&traj);
            throw std::runtime_error("Failed to open output CSV: " + cfg.out_csv);
        }
        ofs << std::setprecision(12);
        ofs << "frame_idx,step,time_ps,center_n_left,center_n_right,center_dn,center_cum_left,center_cum_right,"
               "center_cum_dn,seam_n_left,seam_n_right,seam_dn,seam_cum_left,seam_cum_right,seam_cum_dn\n";

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[gating-config] file=" << cfg.xtc_path << " frames=" << cfg.max_frames
                  << " threads=" << cfg.nthreads << " nx/nz=" << cfg.nx << "/" << cfg.nz
                  << " channel_len_x=" << cfg.channel_length_x << " channel_h_z=" << cfg.channel_height_z
                  << " selection=" << (int)cfg.selection << " out_csv=" << cfg.out_csv << "\n";
        std::cout << "[topology] nmol=" << expected_nmols << " water=" << cfg.nsol << " na=" << cfg.nna
                  << " cl=" << cfg.ncl << " natoms=" << expected_natoms << "\n";

        int64_t center_cum_left = 0;
        int64_t center_cum_right = 0;
        int64_t seam_cum_left = 0;
        int64_t seam_cum_right = 0;

        int frames_done = 0;
        while (frames_done < cfg.max_frames && (rc = xtc_read_next(&traj)) == exdrOK) {
            const XtcFrame* xr = xtc_tail(&traj);
            if (!xr) {
                xtc_close(&traj);
                throw std::runtime_error("Null frame from reader.");
            }

            fr.step = xr->step;
            fr.time_ps = xr->time;
            fr.pbc = simio::Pbc3D(xr->box[0][0], xr->box[1][1], xr->box[2][2]);
            if (fr.pbc.L[0] <= 0.0 || fr.pbc.L[1] <= 0.0 || fr.pbc.L[2] <= 0.0) {
                xtc_close(&traj);
                throw std::runtime_error("Invalid box lengths.");
            }

            for (int ai = 0; ai < expected_natoms; ++ai) {
                fr.atoms.x[(size_t)ai] = xr->x[ai][0];
                fr.atoms.y[(size_t)ai] = xr->x[ai][1];
                fr.atoms.z[(size_t)ai] = xr->x[ai][2];
            }

            pipe.process_frame(topo, fr);

            const auto& c = center_ptr->frames().back();
            const auto& s = wrap_ptr->frames().back();

            center_cum_left += c.n_left;
            center_cum_right += c.n_right;
            seam_cum_left += s.n_left;
            seam_cum_right += s.n_right;

            const int64_t center_cum_dn = center_cum_right - center_cum_left;
            const int64_t seam_cum_dn = seam_cum_right - seam_cum_left;

            ofs << frames_done << "," << fr.step << "," << fr.time_ps << "," << c.n_left << "," << c.n_right
                << "," << c.dn << "," << center_cum_left << "," << center_cum_right << "," << center_cum_dn
                << "," << s.n_left << "," << s.n_right << "," << s.dn << "," << seam_cum_left << ","
                << seam_cum_right << "," << seam_cum_dn << "\n";

            std::cout << "[gating frame " << frames_done << "] step=" << fr.step << " time_ps=" << fr.time_ps
                      << " center(L/R/dN)=(" << c.n_left << "/" << c.n_right << "/" << c.dn << ")"
                      << " seam(L/R/dN)=(" << s.n_left << "/" << s.n_right << "/" << s.dn << ")\n";

            ++frames_done;
        }

        xtc_close(&traj);
        if (rc != exdrOK && rc != exdrENDOFFILE) {
            throw std::runtime_error("Reader error rc=" + std::to_string(rc));
        }
        if (frames_done == 0) {
            throw std::runtime_error("No frames were processed.");
        }

        pipe.finalize();
        const auto& center_total = center_ptr->cumulative();
        const auto& seam_total = wrap_ptr->cumulative();
        std::cout << "[gating cumulative] center(L/R/dN)=(" << center_total.n_left << "/" << center_total.n_right
                  << "/" << center_total.dn << ") seam(L/R/dN)=(" << seam_total.n_left << "/"
                  << seam_total.n_right << "/" << seam_total.dn << ")\n";
        std::cout << "Processed " << frames_done << " frame(s). Wrote " << cfg.out_csv << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

