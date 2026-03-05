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

enum Species : int { Water = 0, Na = 1, Cl = 2, SpeciesN = 3 };

struct CliConfig {
    const char* xtc_path = nullptr;
    int max_frames = 100;
    int nthreads = 4;
    double grid_cell_nm = 0.5;
    int nsol = 5896;
    int nna = 110;
    int ncl = 110;
    double zmin = 0.901;
    double zmax = 1.801;
    int nx = 100;
    std::string out_csv = "density_x_mean.csv";
};

struct FrameDensity {
    std::vector<int64_t> counts[(size_t)SpeciesN];
    int64_t selected[(size_t)SpeciesN]{0, 0, 0};
    int64_t binned[(size_t)SpeciesN]{0, 0, 0};
    int64_t oob = 0;
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

struct DensityAccumulator {
    int nx = 0;
    int nframes = 0;
    double Lx_ref = 0.0;
    bool has_ref_box = false;
    bool warned_variable_box = false;

    RunningStats1D rho[(size_t)SpeciesN];
    RunningStats1D cnt[(size_t)SpeciesN];

    void init(int nx_) {
        nx = nx_;
        nframes = 0;
        has_ref_box = false;
        warned_variable_box = false;
        for (int s = 0; s < SpeciesN; ++s) {
            rho[(size_t)s].init(nx);
            cnt[(size_t)s].init(nx);
        }
    }
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
        << " [zmin=0.901] [zmax=1.801] [nx=100] [out_csv=density_x_mean.csv]\n";
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
    if (argc > 8) c.zmin = parse_double_arg(argv[8], "zmin");
    if (argc > 9) c.zmax = parse_double_arg(argv[9], "zmax");
    if (argc > 10) c.nx = parse_int_arg(argv[10], "nx");
    if (argc > 11) c.out_csv = argv[11];

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

inline int species_from_type(simio::MolType t) {
    if (t == simio::MolType::Water) return Species::Water;
    if (t == simio::MolType::Cation) return Species::Na;
    if (t == simio::MolType::Anion) return Species::Cl;
    return -1;
}

inline bool in_range_pbc(double xw, double mnw, double mxw) {
    if (mnw <= mxw) return (xw >= mnw && xw < mxw);
    return (xw >= mnw) || (xw < mxw);
}

inline double interval_length_pbc(double mnw, double mxw, double L) {
    if (mnw <= mxw) return mxw - mnw;
    return (mxw + L) - mnw;
}

FrameDensity compute_density_frame(const simio::Topology& topo, const simio::Frame& fr,
                                   const std::vector<simio::MolState>& ms, double zmin, double zmax, int nx) {
    FrameDensity out;
    for (int s = 0; s < SpeciesN; ++s) out.counts[(size_t)s].assign((size_t)nx, 0);

    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    const double dx = Lx / (double)nx;

    const double zminw = fr.pbc.wrap_pos(2, zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, zmax);

    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        const int sid = species_from_type(topo.mols[mid].type);
        if (sid < 0) continue;

        const simio::MolState& st = ms[mid];
        if (!(st.cache.flags & simio::MolCache::HAS_KEY)) continue;

        simio::Vec3d key = st.cache.key_wrapped;
        fr.pbc.wrap_pos3(key);
        const double xw = key.v[0];
        const double zw = key.v[2];

        if (!in_range_pbc(zw, zminw, zmaxw)) continue;

        ++out.selected[(size_t)sid];

        int ix = (int)std::floor(xw / dx);
        if (ix < 0 || ix >= nx) {
            ++out.oob;
            continue;
        }

        ++out.counts[(size_t)sid][(size_t)ix];
        ++out.binned[(size_t)sid];
    }

    return out;
}

double mean_from(const RunningStats1D& st, int i, int n) {
    return st.sum[(size_t)i] / (double)n;
}

double sem_from(const RunningStats1D& st, int i, int n) {
    if (n <= 1) return 0.0;
    const double sum = st.sum[(size_t)i];
    const double sumsq = st.sumsq[(size_t)i];
    const double mean = sum / (double)n;
    double var = (sumsq - (double)n * mean * mean) / (double)(n - 1);
    if (var < 0.0) var = 0.0;
    return std::sqrt(var / (double)n);
}

void write_density_csv(const std::string& path, const DensityAccumulator& acc) {
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("Failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "x_center_nm,"
           "rho_water_mean,rho_water_sem,rho_na_mean,rho_na_sem,rho_cl_mean,rho_cl_sem,"
           "count_water_mean,count_water_sem,count_na_mean,count_na_sem,count_cl_mean,count_cl_sem\n";

    const double dx_ref = acc.Lx_ref / (double)acc.nx;
    for (int i = 0; i < acc.nx; ++i) {
        const double x_center = ((double)i + 0.5) * dx_ref;
        ofs << x_center << ","
            << mean_from(acc.rho[(size_t)Species::Water], i, acc.nframes) << ","
            << sem_from(acc.rho[(size_t)Species::Water], i, acc.nframes) << ","
            << mean_from(acc.rho[(size_t)Species::Na], i, acc.nframes) << ","
            << sem_from(acc.rho[(size_t)Species::Na], i, acc.nframes) << ","
            << mean_from(acc.rho[(size_t)Species::Cl], i, acc.nframes) << ","
            << sem_from(acc.rho[(size_t)Species::Cl], i, acc.nframes) << ","
            << mean_from(acc.cnt[(size_t)Species::Water], i, acc.nframes) << ","
            << sem_from(acc.cnt[(size_t)Species::Water], i, acc.nframes) << ","
            << mean_from(acc.cnt[(size_t)Species::Na], i, acc.nframes) << ","
            << sem_from(acc.cnt[(size_t)Species::Na], i, acc.nframes) << ","
            << mean_from(acc.cnt[(size_t)Species::Cl], i, acc.nframes) << ","
            << sem_from(acc.cnt[(size_t)Species::Cl], i, acc.nframes) << "\n";
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

        DensityAccumulator acc;
        acc.init(cfg.nx);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[density-config] file=" << cfg.xtc_path << " frames=" << cfg.max_frames
                  << " threads=" << cfg.nthreads << " nx=" << cfg.nx << " zmin=" << cfg.zmin
                  << " zmax=" << cfg.zmax << " out_csv=" << cfg.out_csv << "\n";
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

            const FrameDensity fd = compute_density_frame(topo, fr, ms, cfg.zmin, cfg.zmax, cfg.nx);
            const bool ok = (fd.selected[(size_t)Species::Water] == fd.binned[(size_t)Species::Water]) &&
                            (fd.selected[(size_t)Species::Na] == fd.binned[(size_t)Species::Na]) &&
                            (fd.selected[(size_t)Species::Cl] == fd.binned[(size_t)Species::Cl]) &&
                            (fd.oob == 0);

            const double zminw = fr.pbc.wrap_pos(2, cfg.zmin);
            const double zmaxw = fr.pbc.wrap_pos(2, cfg.zmax);
            const double zlen = interval_length_pbc(zminw, zmaxw, fr.pbc.L[2]);
            const double dx = fr.pbc.L[0] / (double)cfg.nx;
            const double denom = dx * fr.pbc.L[1] * zlen;

            if (denom <= 0.0) {
                xtc_close(&traj);
                throw std::runtime_error("Non-positive normalization denominator.");
            }

            if (!acc.has_ref_box) {
                acc.has_ref_box = true;
                acc.Lx_ref = fr.pbc.L[0];
            } else if (!acc.warned_variable_box) {
                if (std::abs(fr.pbc.L[0] - acc.Lx_ref) > 1e-6) {
                    std::cout << "[warn] Lx changed from " << acc.Lx_ref << " to " << fr.pbc.L[0]
                              << "; x_center in CSV uses first-frame Lx\n";
                    acc.warned_variable_box = true;
                }
            }

            for (int i = 0; i < cfg.nx; ++i) {
                for (int s = 0; s < SpeciesN; ++s) {
                    const double c = (double)fd.counts[(size_t)s][(size_t)i];
                    acc.cnt[(size_t)s].add(i, c);
                    acc.rho[(size_t)s].add(i, c / denom);
                }
            }
            ++acc.nframes;

            std::cout << "[density frame " << frames_done << "] step=" << fr.step << " time_ps=" << fr.time_ps
                      << " selected(water/na/cl)=" << fd.selected[(size_t)Species::Water] << "/"
                      << fd.selected[(size_t)Species::Na] << "/" << fd.selected[(size_t)Species::Cl]
                      << " binned(water/na/cl)=" << fd.binned[(size_t)Species::Water] << "/"
                      << fd.binned[(size_t)Species::Na] << "/" << fd.binned[(size_t)Species::Cl]
                      << " oob=" << fd.oob << " status=" << (ok ? "OK" : "FAIL") << "\n";

            if (!ok) {
                xtc_close(&traj);
                throw std::runtime_error("Density invariant failed.");
            }

            ++frames_done;
        }

        xtc_close(&traj);
        if (rc != exdrOK && rc != exdrENDOFFILE) {
            throw std::runtime_error("Reader error rc=" + std::to_string(rc));
        }
        if (acc.nframes == 0) {
            throw std::runtime_error("No frames were processed.");
        }

        write_density_csv(cfg.out_csv, acc);
        std::cout << "Processed " << acc.nframes << " frame(s). Wrote " << cfg.out_csv << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

