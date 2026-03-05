#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "simio/simio.hpp"
#include "xtc.h"

namespace {

struct CliConfig {
    const char* xtc_path = nullptr;
    int max_frames = 100;
    int nthreads = 4;
    double grid_cell_nm = 0.5;
    int nsol = 5896;
    int nna = 110;
    int ncl = 110;
    double xmin = 7.11;
    double xmax = 12.89;
    double zmin = 0.901;
    double zmax = 1.801;
    int nx = 100;
    std::string out_csv = "dipole_x_mean.csv";
};

struct RunningStats1D {
    std::vector<double> sum;
    std::vector<double> sumsq;

    void init(int n) {
        sum.assign((size_t)n, 0.0);
        sumsq.assign((size_t)n, 0.0);
    }

    void add(int i, double v) {
        sum[(size_t)i] += v;
        sumsq[(size_t)i] += v * v;
    }
};

struct FrameDipoleProfile {
    std::vector<double> sum_muz_fold;
    std::vector<double> sum_muz_raw;
    std::vector<double> sum_mux;
    std::vector<int64_t> count;
    int64_t selected_water = 0;
    int64_t binned_water = 0;
    int64_t oob = 0;
    int64_t zero_norm = 0;
};

struct IntervalMap {
    bool inside = false;
    double u = 0.0;
    double length = 0.0;
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

void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " <trajectory.xtc> [max_frames=100] [threads=4] [grid_cell_nm=0.5] [nsol=5896] [nna=110] [ncl=110]"
        << " [xmin=7.11] [xmax=12.89] [zmin=0.901] [zmax=1.801] [nx=100] [out_csv=dipole_x_mean.csv]\n"
        << "       Set xmin==xmax (e.g. 0 0) to use whole x [0,Lx).\n";
}

CliConfig parse_cli(int argc, char** argv) {
    if (argc < 2) throw std::runtime_error("Missing trajectory path.");

    CliConfig c;
    c.xtc_path = argv[1];
    if (argc > 2) c.max_frames = parse_int_arg(argv[2], "max_frames");
    if (argc > 3) c.nthreads = parse_int_arg(argv[3], "threads");
    if (argc > 4) c.grid_cell_nm = parse_double_arg(argv[4], "grid_cell_nm");
    if (argc > 5) c.nsol = parse_int_arg(argv[5], "nsol");
    if (argc > 6) c.nna = parse_int_arg(argv[6], "nna");
    if (argc > 7) c.ncl = parse_int_arg(argv[7], "ncl");
    if (argc > 8) c.xmin = parse_double_arg(argv[8], "xmin");
    if (argc > 9) c.xmax = parse_double_arg(argv[9], "xmax");
    if (argc > 10) c.zmin = parse_double_arg(argv[10], "zmin");
    if (argc > 11) c.zmax = parse_double_arg(argv[11], "zmax");
    if (argc > 12) c.nx = parse_int_arg(argv[12], "nx");
    if (argc > 13) c.out_csv = argv[13];

    if (c.max_frames <= 0 || c.nthreads <= 0 || c.grid_cell_nm <= 0.0 || c.nsol < 0 || c.nna < 0 ||
        c.ncl < 0 || c.nx <= 0) {
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

IntervalMap map_on_pbc_interval(double xw, double mnw, double mxw, double L) {
    IntervalMap out{};
    if (L <= 0.0) return out;
    if (mnw <= mxw) {
        out.length = mxw - mnw;
        if (out.length <= 0.0) return out;
        if (xw < mnw || xw >= mxw) return out;
        out.inside = true;
        out.u = xw - mnw;
        return out;
    }
    out.length = (mxw + L) - mnw;
    if (out.length <= 0.0) return out;
    if (xw >= mnw) {
        out.inside = true;
        out.u = xw - mnw;
        return out;
    }
    if (xw < mxw) {
        out.inside = true;
        out.u = (xw + L) - mnw;
        return out;
    }
    return out;
}

double interval_length_pbc(double mnw, double mxw, double L) {
    if (mnw <= mxw) return mxw - mnw;
    return (mxw + L) - mnw;
}

double sem_from(const RunningStats1D& s, int i, int n) {
    if (n <= 1) return 0.0;
    const double sum = s.sum[(size_t)i];
    const double sumsq = s.sumsq[(size_t)i];
    const double mean = sum / (double)n;
    double var = (sumsq - (double)n * mean * mean) / (double)(n - 1);
    if (var < 0.0) var = 0.0;
    return std::sqrt(var / (double)n);
}

FrameDipoleProfile compute_frame_profile(const simio::Topology& topo, const simio::Frame& fr,
                                         const std::vector<simio::MolState>& ms, double xmin, double xmax,
                                         double zmin, double zmax, int nx) {
    FrameDipoleProfile out;
    out.sum_muz_fold.assign((size_t)nx, 0.0);
    out.sum_muz_raw.assign((size_t)nx, 0.0);
    out.sum_mux.assign((size_t)nx, 0.0);
    out.count.assign((size_t)nx, 0);

    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    const bool whole_x = (xmax <= xmin);
    const double xminw = whole_x ? 0.0 : fr.pbc.wrap_pos(0, xmin);
    const double xmaxw = whole_x ? 0.0 : fr.pbc.wrap_pos(0, xmax);
    const double zminw = fr.pbc.wrap_pos(2, zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, zmax);
    const double xlen = whole_x ? Lx : interval_length_pbc(xminw, xmaxw, Lx);
    const double zlen = interval_length_pbc(zminw, zmaxw, Lz);
    const double dx = xlen / (double)nx;
    const double z_center = fr.pbc.wrap_pos(2, zminw + 0.5 * zlen);
    constexpr double kNormEps = 1e-14;

    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        const simio::MolSpan& m = topo.mols[mid];
        if (m.type != simio::MolType::Water || m.natoms < 3) continue;

        const simio::MolState& st = ms[mid];
        if (!(st.cache.flags & simio::MolCache::HAS_KEY)) continue;

        simio::Vec3d key = st.cache.key_wrapped;
        fr.pbc.wrap_pos3(key);
        const double xw = key.v[0];
        const double zw = key.v[2];

        IntervalMap xm{};
        if (whole_x) {
            xm.inside = true;
            xm.u = xw;
            xm.length = Lx;
        } else {
            xm = map_on_pbc_interval(xw, xminw, xmaxw, Lx);
        }
        const IntervalMap zm = map_on_pbc_interval(zw, zminw, zmaxw, Lz);
        if (!(xm.inside && zm.inside)) continue;

        ++out.selected_water;

        int ix = (int)std::floor(xm.u / dx);
        if (ix < 0 || ix >= nx) {
            ++out.oob;
            continue;
        }

        simio::Vec3d rO{};
        simio::Vec3d h1u{};
        simio::Vec3d h2u{};
        if ((st.cache.flags & simio::MolCache::HAS_REF) && (st.cache.flags & simio::MolCache::HAS_SITES)) {
            // Reuse geometry prepared once in BuildIntrinsicsAndGrid.
            rO = st.cache.ref_wrapped;
            h1u = st.cache.sites_u[1];
            h2u = st.cache.sites_u[2];
        } else {
            // Fallback for robustness if cache flags are unavailable.
            const int O = m.first + 0;
            const int H1 = m.first + 1;
            const int H2 = m.first + 2;
            rO = fr.atoms.pos(O);
            const simio::Vec3d rH1 = fr.atoms.pos(H1);
            const simio::Vec3d rH2 = fr.atoms.pos(H2);
            fr.pbc.wrap_pos3(rO);
            h1u = simio::pbc_unwrap_to_ref(fr.pbc, rO, rH1);
            h2u = simio::pbc_unwrap_to_ref(fr.pbc, rO, rH2);
        }
        simio::Vec3d mu = (h1u - rO) + (h2u - rO);

        const double mnorm = simio::norm3(mu);
        if (mnorm <= kNormEps) {
            ++out.zero_norm;
            mu = simio::Vec3d{};
        } else {
            mu = mu / mnorm;
        }

        const double mux = mu.v[0];
        const double muz_raw = mu.v[2];
        const bool lower_half = (zm.u < 0.5 * zlen);
        const double muz_fold = lower_half ? muz_raw : -muz_raw;

        out.sum_mux[(size_t)ix] += mux;
        out.sum_muz_raw[(size_t)ix] += muz_raw;
        out.sum_muz_fold[(size_t)ix] += muz_fold;
        ++out.count[(size_t)ix];
        ++out.binned_water;
    }

    return out;
}

void write_csv(const std::string& path, const std::vector<double>& x_centers, const RunningStats1D& muz_raw,
               const RunningStats1D& muz_fold,
               const RunningStats1D& mux, const RunningStats1D& cnt, int nframes) {
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("Failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "x_center_nm,muz_mean,muz_sem,muz_fold_mean,muz_fold_sem,mux_mean,mux_sem,count_mean,count_sem\n";

    for (size_t i = 0; i < x_centers.size(); ++i) {
        const double m_muz_raw = muz_raw.sum[i] / (double)nframes;
        const double m_muz_fold = muz_fold.sum[i] / (double)nframes;
        const double m_mux = mux.sum[i] / (double)nframes;
        const double m_cnt = cnt.sum[i] / (double)nframes;
        ofs << x_centers[i] << "," << m_muz_raw << "," << sem_from(muz_raw, (int)i, nframes) << ","
            << m_muz_fold << "," << sem_from(muz_fold, (int)i, nframes) << "," << m_mux << ","
            << sem_from(mux, (int)i, nframes) << "," << m_cnt << "," << sem_from(cnt, (int)i, nframes)
            << "\n";
    }
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

        std::vector<simio::MolState> ms((size_t)expected_nmols);
        auto pipe = simio::make_default_pipeline(cfg.nthreads, cfg.grid_cell_nm);

        RunningStats1D st_muz;
        RunningStats1D st_muz_raw;
        RunningStats1D st_mux;
        RunningStats1D st_cnt;
        st_muz_raw.init(cfg.nx);
        st_muz.init(cfg.nx);
        st_mux.init(cfg.nx);
        st_cnt.init(cfg.nx);

        std::vector<double> x_centers;
        x_centers.assign((size_t)cfg.nx, 0.0);
        bool has_ref_x_centers = false;
        bool warned_x_grid_change = false;
        double ref_xlen = 0.0;
        double ref_xminw = 0.0;

        std::cout << std::fixed << std::setprecision(6);
        const bool whole_x = (cfg.xmax <= cfg.xmin);
        std::cout << "[dipole-config] file=" << cfg.xtc_path << " frames=" << cfg.max_frames
                  << " threads=" << cfg.nthreads << " nx=" << cfg.nx << " x=[" << cfg.xmin << "," << cfg.xmax
                  << "] x_mode=" << (whole_x ? "whole_box" : "roi") << " z=[" << cfg.zmin << "," << cfg.zmax
                  << "] out_csv=" << cfg.out_csv << "\n";
        std::cout << "[topology] nmol=" << expected_nmols << " water=" << cfg.nsol << " na=" << cfg.nna
                  << " cl=" << cfg.ncl << " natoms=" << expected_natoms << "\n";

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

            pipe.process_frame(topo, fr, ms);

            const double xminw = whole_x ? 0.0 : fr.pbc.wrap_pos(0, cfg.xmin);
            const double xmaxw = whole_x ? 0.0 : fr.pbc.wrap_pos(0, cfg.xmax);
            const double xlen = whole_x ? fr.pbc.L[0] : interval_length_pbc(xminw, xmaxw, fr.pbc.L[0]);
            if (xlen <= 0.0) {
                xtc_close(&traj);
                throw std::runtime_error("Invalid x interval length.");
            }

            if (!has_ref_x_centers) {
                has_ref_x_centers = true;
                ref_xlen = xlen;
                ref_xminw = xminw;
                const double dx = ref_xlen / (double)cfg.nx;
                for (int i = 0; i < cfg.nx; ++i) {
                    x_centers[(size_t)i] = fr.pbc.wrap_pos(0, ref_xminw + ((double)i + 0.5) * dx);
                }
            } else if (!warned_x_grid_change && std::abs(xlen - ref_xlen) > 1e-6) {
                std::cout << "[warn] x interval length changed from " << ref_xlen << " to " << xlen
                          << "; CSV x centers use first-frame mapping\n";
                warned_x_grid_change = true;
            }

            const FrameDipoleProfile fd =
                compute_frame_profile(topo, fr, ms, cfg.xmin, cfg.xmax, cfg.zmin, cfg.zmax, cfg.nx);
            const bool ok = (fd.selected_water == fd.binned_water) && (fd.oob == 0);

            for (int i = 0; i < cfg.nx; ++i) {
                const int64_t c = fd.count[(size_t)i];
                const double mean_muz_raw = (c > 0) ? (fd.sum_muz_raw[(size_t)i] / (double)c) : 0.0;
                const double mean_muz_fold = (c > 0) ? (fd.sum_muz_fold[(size_t)i] / (double)c) : 0.0;
                const double mean_mux = (c > 0) ? (fd.sum_mux[(size_t)i] / (double)c) : 0.0;

                st_muz_raw.add(i, mean_muz_raw);
                st_muz.add(i, mean_muz_fold);
                st_mux.add(i, mean_mux);
                st_cnt.add(i, (double)c);
            }

            std::cout << "[dipole frame " << frames_done << "] step=" << fr.step << " time_ps=" << fr.time_ps
                      << " selected_water=" << fd.selected_water << " binned_water=" << fd.binned_water
                      << " oob=" << fd.oob << " zero_norm=" << fd.zero_norm << " status="
                      << (ok ? "OK" : "FAIL") << "\n";

            if (!ok) {
                xtc_close(&traj);
                throw std::runtime_error("Dipole invariant failed.");
            }

            ++frames_done;
        }

        xtc_close(&traj);
        if (rc != exdrOK && rc != exdrENDOFFILE) {
            throw std::runtime_error("Reader error rc=" + std::to_string(rc));
        }
        if (frames_done == 0) {
            throw std::runtime_error("No frames were processed.");
        }

        write_csv(cfg.out_csv, x_centers, st_muz_raw, st_muz, st_mux, st_cnt, frames_done);
        std::cout << "Processed " << frames_done << " frame(s). Wrote " << cfg.out_csv << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
