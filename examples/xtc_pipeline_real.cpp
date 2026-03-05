#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "simio/simio.hpp"
#include "xtc.h"

namespace {

constexpr double kEps = 1e-10;
constexpr std::int64_t kTlsWarnThreshold = 5000000;

struct GridValidation {
    int entries = 0;
    int duplicates = 0;
    int missing = 0;
    int out_of_range = 0;
};

struct GridBucketStats {
    int nonempty = 0;
    int max_bucket = 0;
    double mean_nonempty = 0.0;
};

struct StepStats {
    double min_step_mi = 0.0;
    double max_step_mi = 0.0;
    double mean_step_mi = 0.0;
    int pbc_cross_total = 0;
    int pbc_cross_x = 0;
    int pbc_cross_y = 0;
    int pbc_cross_z = 0;
    double max_corr_component = 0.0;
    int max_corr_mol = -1;
    int max_corr_axis = -1;
    simio::Vec3d max_corr_raw{};
    simio::Vec3d max_corr_mi{};
    simio::Vec3d max_corr_diff{};
    double max_corr_norm = 0.0;
    int max_step_mol = -1;
    simio::Vec3d max_step_vec{};
    int k_hist_x_m1 = 0;
    int k_hist_x_0 = 0;
    int k_hist_x_p1 = 0;
    int k_hist_x_other = 0;
    int k_hist_y_m1 = 0;
    int k_hist_y_0 = 0;
    int k_hist_y_p1 = 0;
    int k_hist_y_other = 0;
    int k_hist_z_m1 = 0;
    int k_hist_z_0 = 0;
    int k_hist_z_p1 = 0;
    int k_hist_z_other = 0;
};

struct WrappedViolations {
    int key = 0;
    int ref = 0;
};

struct InjectConfig {
    int mol_id = -1;
    int axis = 0;
    double frac_L = 0.8;
    int frame_id = 1;
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

GridValidation validate_grid(const simio::Frame::MolGrid& grid, int nmol) {
    GridValidation out;
    std::vector<int> seen((size_t)nmol, 0);

    for (const auto& bucket : grid.cells) {
        out.entries += (int)bucket.size();
        for (int mid : bucket) {
            if (mid < 0 || mid >= nmol) {
                ++out.out_of_range;
                continue;
            }
            ++seen[(size_t)mid];
        }
    }

    for (int c : seen) {
        if (c == 0) ++out.missing;
        if (c > 1) out.duplicates += (c - 1);
    }

    return out;
}

GridBucketStats compute_grid_bucket_stats(const simio::Frame::MolGrid& grid) {
    GridBucketStats out;
    int sum_nonempty = 0;

    for (const auto& bucket : grid.cells) {
        const int sz = (int)bucket.size();
        if (sz > out.max_bucket) out.max_bucket = sz;
        if (sz > 0) {
            ++out.nonempty;
            sum_nonempty += sz;
        }
    }

    if (out.nonempty > 0) {
        out.mean_nonempty = (double)sum_nonempty / (double)out.nonempty;
    }
    return out;
}

StepStats compute_step_stats(const simio::Pbc3D& pbc, const std::vector<simio::MolState>& ms) {
    StepStats out;
    if (ms.empty()) return out;

    out.min_step_mi = std::numeric_limits<double>::max();
    double sum = 0.0;

    for (int mid = 0; mid < (int)ms.size(); ++mid) {
        const auto& st = ms[(size_t)mid];
        const double step = st.time.last_step_norm;
        if (step < out.min_step_mi) out.min_step_mi = step;
        if (step > out.max_step_mi) {
            out.max_step_mi = step;
            out.max_step_mol = mid;
            out.max_step_vec = st.time.last_dr;
        }
        sum += step;

        bool crossed_mol = false;
        for (int a = 0; a < 3; ++a) {
            const double diff = st.time.last_raw_dr.v[a] - st.time.last_dr.v[a];
            const double corr = std::abs(diff);

            long long k_round = 0;
            if (pbc.L[a] > 0.0) {
                k_round = std::llround(diff / pbc.L[a]);
            }
            if (a == 0) {
                if (k_round == -1) ++out.k_hist_x_m1;
                else if (k_round == 0) ++out.k_hist_x_0;
                else if (k_round == 1) ++out.k_hist_x_p1;
                else ++out.k_hist_x_other;
            } else if (a == 1) {
                if (k_round == -1) ++out.k_hist_y_m1;
                else if (k_round == 0) ++out.k_hist_y_0;
                else if (k_round == 1) ++out.k_hist_y_p1;
                else ++out.k_hist_y_other;
            } else {
                if (k_round == -1) ++out.k_hist_z_m1;
                else if (k_round == 0) ++out.k_hist_z_0;
                else if (k_round == 1) ++out.k_hist_z_p1;
                else ++out.k_hist_z_other;
            }

            if (corr > out.max_corr_component) {
                out.max_corr_component = corr;
                out.max_corr_mol = mid;
                out.max_corr_axis = a;
                out.max_corr_raw = st.time.last_raw_dr;
                out.max_corr_mi = st.time.last_dr;
                out.max_corr_diff =
                    simio::Vec3d{{st.time.last_raw_dr.v[0] - st.time.last_dr.v[0],
                                  st.time.last_raw_dr.v[1] - st.time.last_dr.v[1],
                                  st.time.last_raw_dr.v[2] - st.time.last_dr.v[2]}};
                out.max_corr_norm = simio::norm3(out.max_corr_diff);
            }

            const bool corrected_axis = (pbc.L[a] > 0.0) && (corr > 0.5 * pbc.L[a]);
            if (corrected_axis) {
                crossed_mol = true;
                if (a == 0) ++out.pbc_cross_x;
                if (a == 1) ++out.pbc_cross_y;
                if (a == 2) ++out.pbc_cross_z;
            }
        }
        if (crossed_mol) ++out.pbc_cross_total;
    }

    out.mean_step_mi = sum / (double)ms.size();
    return out;
}

WrappedViolations count_wrapped_violations(const simio::Frame& fr, const std::vector<simio::MolState>& ms) {
    WrappedViolations out;

    for (const auto& st : ms) {
        if ((st.cache.flags & simio::MolCache::HAS_KEY) &&
            !simio::is_wrapped3_box(st.cache.key_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2])) {
            ++out.key;
        }
        if ((st.cache.flags & simio::MolCache::HAS_REF) &&
            !simio::is_wrapped3_box(st.cache.ref_wrapped, fr.pbc.L[0], fr.pbc.L[1], fr.pbc.L[2])) {
            ++out.ref;
        }
    }

    return out;
}

void inject_pbc_shift_if_requested(const InjectConfig& inj, const simio::Topology& topo, simio::Frame& fr,
                                   int frame_id) {
    if (inj.mol_id < 0 || frame_id != inj.frame_id) return;

    const simio::MolSpan& m = topo.mols[(size_t)inj.mol_id];
    const double L = fr.pbc.L[inj.axis];
    const double delta = inj.frac_L * L;

    for (int k = 0; k < m.natoms; ++k) {
        const int ai = m.first + k;
        if (inj.axis == 0) {
            fr.atoms.x[(size_t)ai] = fr.pbc.wrap_pos(0, fr.atoms.x[(size_t)ai] + delta);
        } else if (inj.axis == 1) {
            fr.atoms.y[(size_t)ai] = fr.pbc.wrap_pos(1, fr.atoms.y[(size_t)ai] + delta);
        } else {
            fr.atoms.z[(size_t)ai] = fr.pbc.wrap_pos(2, fr.atoms.z[(size_t)ai] + delta);
        }
    }

    std::cout << "[inject] frame=" << frame_id << " mol_id=" << inj.mol_id << " axis=" << inj.axis
              << " frac_L=" << inj.frac_L << " applied_to_atoms=" << m.natoms << "\n";
}

const char* axis_name(int a) {
    if (a == 0) return "x";
    if (a == 1) return "y";
    if (a == 2) return "z";
    return "?";
}

void print_vec3(const simio::Vec3d& v) {
    std::cout << "(" << v.v[0] << "," << v.v[1] << "," << v.v[2] << ")";
}

void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " <trajectory.xtc> [max_frames=10] [threads=4] [grid_cell_nm=0.5] [nsol=5896] [nna=110] [ncl=110]"
        << " [inject_mol_id=-1] [inject_axis=0] [inject_fracL=0.8] [inject_frame=1]\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char* path = argv[1];
    const int max_frames = (argc > 2) ? parse_int_arg(argv[2], "max_frames") : 10;
    const int nthreads = (argc > 3) ? parse_int_arg(argv[3], "threads") : 4;
    const double grid_cell_nm = (argc > 4) ? parse_double_arg(argv[4], "grid_cell_nm") : 0.5;
    const int nsol = (argc > 5) ? parse_int_arg(argv[5], "nsol") : 5896;
    const int nna = (argc > 6) ? parse_int_arg(argv[6], "nna") : 110;
    const int ncl = (argc > 7) ? parse_int_arg(argv[7], "ncl") : 110;

    InjectConfig inj;
    inj.mol_id = (argc > 8) ? parse_int_arg(argv[8], "inject_mol_id") : -1;
    inj.axis = (argc > 9) ? parse_int_arg(argv[9], "inject_axis") : 0;
    inj.frac_L = (argc > 10) ? parse_double_arg(argv[10], "inject_fracL") : 0.8;
    inj.frame_id = (argc > 11) ? parse_int_arg(argv[11], "inject_frame") : 1;

    if (max_frames <= 0 || nthreads <= 0 || grid_cell_nm <= 0.0 || nsol < 0 || nna < 0 || ncl < 0) {
        throw std::runtime_error("All numeric inputs must be positive (or non-negative for counts).");
    }
    if (inj.axis < 0 || inj.axis > 2) {
        throw std::runtime_error("inject_axis must be 0, 1, or 2.");
    }
    if (inj.frac_L <= 0.0) {
        throw std::runtime_error("inject_fracL must be > 0.");
    }
    if (inj.frame_id < 0) {
        throw std::runtime_error("inject_frame must be >= 0.");
    }

    const int expected_natoms = nsol * 3 + nna + ncl;
    const int expected_nmols = nsol + nna + ncl;

    simio::Topology topo = build_topology(nsol, nna, ncl);
    if ((int)topo.mols.size() != expected_nmols) {
        throw std::runtime_error("Topology molecule count mismatch.");
    }
    if (inj.mol_id >= expected_nmols) {
        throw std::runtime_error("inject_mol_id out of range for topology.");
    }

    XtcTraj traj{};
    int rc = xtc_open(&traj, path, /*cap=*/1);
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

    std::vector<simio::MolState> ms((size_t)expected_nmols);
    auto pipe = simio::make_default_pipeline(nthreads, grid_cell_nm);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[topology] nmol=" << expected_nmols << " water=" << nsol << " na=" << nna
              << " cl=" << ncl << " expected_natoms=" << expected_natoms << " file_natoms=" << traj.natoms
              << "\n";
    std::cout << "[keys] grid_key=COM time_key=COM ref_key=O\n";
    std::cout << "[config] file=" << path << " max_frames=" << max_frames << " threads=" << nthreads
              << " grid_cell_nm=" << grid_cell_nm;
    if (inj.mol_id >= 0) {
        std::cout << " inject=(mol=" << inj.mol_id << ",axis=" << inj.axis << ",fracL=" << inj.frac_L
                  << ",frame=" << inj.frame_id << ")";
    }
    std::cout << "\n";

    int frames_done = 0;
    bool non_orth_warned = false;
    bool printed_grid_once = false;
    bool has_prev_time = false;
    double prev_time_ps = 0.0;
    double dt_sum_ps = 0.0;
    int dt_samples = 0;

    while (frames_done < max_frames && (rc = xtc_read_next(&traj)) == exdrOK) {
        const XtcFrame* xr = xtc_tail(&traj);
        if (!xr) {
            std::cerr << "Null frame returned by reader\n";
            xtc_close(&traj);
            return 1;
        }

        fr.step = xr->step;
        fr.time_ps = xr->time;

        if (!non_orth_warned) {
            const bool offdiag = (std::abs(xr->box[0][1]) > 1e-8f || std::abs(xr->box[0][2]) > 1e-8f ||
                                  std::abs(xr->box[1][0]) > 1e-8f || std::abs(xr->box[1][2]) > 1e-8f ||
                                  std::abs(xr->box[2][0]) > 1e-8f || std::abs(xr->box[2][1]) > 1e-8f);
            if (offdiag) {
                std::cout << "[warn] non-orthogonal box detected; using diagonal lengths only for orthorhombic PBC\n";
            }
            non_orth_warned = true;
        }

        fr.pbc = simio::Pbc3D(xr->box[0][0], xr->box[1][1], xr->box[2][2]);
        if (fr.pbc.L[0] <= 0.0 || fr.pbc.L[1] <= 0.0 || fr.pbc.L[2] <= 0.0) {
            std::cerr << "Invalid box lengths in frame " << frames_done << "\n";
            xtc_close(&traj);
            return 1;
        }

        fr.grid.init(fr.pbc, grid_cell_nm);
        if (!printed_grid_once) {
            const std::int64_t tls_vectors = (std::int64_t)fr.grid.ncells() * (std::int64_t)nthreads;
            std::cout << "[grid] L=(" << fr.pbc.L[0] << "," << fr.pbc.L[1] << "," << fr.pbc.L[2]
                      << ") cell=" << grid_cell_nm << " n=(" << fr.grid.nx << "," << fr.grid.ny << ","
                      << fr.grid.nz << ") ncells=" << fr.grid.ncells() << " tls_vectors=" << tls_vectors
                      << "\n";
            if (tls_vectors > kTlsWarnThreshold) {
                std::cout << "[warn] large TLS footprint nthreads*ncells=" << tls_vectors << "\n";
            }
            printed_grid_once = true;
        }

        for (int ai = 0; ai < expected_natoms; ++ai) {
            fr.atoms.x[(size_t)ai] = xr->x[ai][0];
            fr.atoms.y[(size_t)ai] = xr->x[ai][1];
            fr.atoms.z[(size_t)ai] = xr->x[ai][2];
        }

        inject_pbc_shift_if_requested(inj, topo, fr, frames_done);

        pipe.process_frame(topo, fr, ms);

        int key_count = 0;
        for (const auto& st : ms) {
            if (st.cache.flags & simio::MolCache::HAS_KEY) ++key_count;
        }

        const WrappedViolations wv = count_wrapped_violations(fr, ms);
        const GridValidation gv = validate_grid(fr.grid, expected_nmols);
        const GridBucketStats gb = compute_grid_bucket_stats(fr.grid);
        const StepStats ss = compute_step_stats(fr.pbc, ms);

        double last_dt_ps = 0.0;
        if (has_prev_time) {
            last_dt_ps = fr.time_ps - prev_time_ps;
            dt_sum_ps += last_dt_ps;
            ++dt_samples;
        }
        prev_time_ps = fr.time_ps;
        has_prev_time = true;

        const bool ok = (key_count == expected_nmols) && (wv.key == 0) && (wv.ref == 0) &&
                        (gv.entries == expected_nmols) && (gv.duplicates == 0) && (gv.missing == 0) &&
                        (gv.out_of_range == 0);

        std::cout << "[frame " << frames_done << "] step=" << fr.step << " time_ps=" << fr.time_ps << "\n";
        std::cout << "  keys=" << key_count << "/" << expected_nmols << " key_range_viol=" << wv.key
                  << " ref_range_viol=" << wv.ref << "\n";
        std::cout << "  grid_entries=" << gv.entries << " dup=" << gv.duplicates << " miss=" << gv.missing
                  << " out_of_range=" << gv.out_of_range << " status=" << (ok ? "OK" : "FAIL") << "\n";
        std::cout << "  pbc_cross: x=" << ss.pbc_cross_x << " y=" << ss.pbc_cross_y
                  << " z=" << ss.pbc_cross_z << " total=" << ss.pbc_cross_total << "\n";
        std::cout << "  max_corr: mol=" << ss.max_corr_mol << " axis=" << axis_name(ss.max_corr_axis)
                  << " corr=" << ss.max_corr_component << " dr_raw=";
        print_vec3(ss.max_corr_raw);
        std::cout << " dr_mi=";
        print_vec3(ss.max_corr_mi);
        std::cout << " corr_vec=";
        print_vec3(ss.max_corr_diff);
        std::cout << " |corr|=" << ss.max_corr_norm << "\n";
        std::cout << "  max_step_mi: mol=" << ss.max_step_mol << " value=" << ss.max_step_mi << " dr_mi=";
        print_vec3(ss.max_step_vec);
        std::cout << "\n";
        std::cout << "  step_mi(min/max/mean)=" << ss.min_step_mi << "/" << ss.max_step_mi << "/"
                  << ss.mean_step_mi << "\n";
        std::cout << "  k_hist: x(-1/0/+1/other)=" << ss.k_hist_x_m1 << "/" << ss.k_hist_x_0 << "/"
                  << ss.k_hist_x_p1 << "/" << ss.k_hist_x_other
                  << " y(-1/0/+1/other)=" << ss.k_hist_y_m1 << "/" << ss.k_hist_y_0 << "/"
                  << ss.k_hist_y_p1 << "/" << ss.k_hist_y_other
                  << " z(-1/0/+1/other)=" << ss.k_hist_z_m1 << "/" << ss.k_hist_z_0 << "/"
                  << ss.k_hist_z_p1 << "/" << ss.k_hist_z_other << "\n";
        std::cout << "  grid_stats: nonempty=" << gb.nonempty << " max_bucket=" << gb.max_bucket
                  << " mean_nonempty=" << gb.mean_nonempty << "\n";
        if (inj.mol_id >= 0 && inj.mol_id < expected_nmols) {
            const auto& im = ms[(size_t)inj.mol_id];
            std::cout << "  inject_mol " << inj.mol_id << ": dr_raw=";
            print_vec3(im.time.last_raw_dr);
            std::cout << " dr_mi=";
            print_vec3(im.time.last_dr);
            std::cout << "\n";
        }
        if (dt_samples > 0) {
            std::cout << "  time_delta_ps: last=" << last_dt_ps
                      << " mean=" << (dt_sum_ps / (double)dt_samples) << "\n";
        }

        if (!ok) {
            std::cerr << "Validation failed on frame " << frames_done << "\n";
            xtc_close(&traj);
            return 1;
        }

        ++frames_done;
    }

    xtc_close(&traj);

    if (rc != exdrOK && rc != exdrENDOFFILE) {
        std::cerr << "Reader error: rc=" << rc << "\n";
        return 1;
    }

    std::cout << "Processed " << frames_done << " frame(s) successfully.\n";
    return 0;
}
