#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "simio/simio.hpp"
#include "xtc.h"

namespace {

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
    MetricN = 9
};

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
    double r_cw = 0.35;  // Na-O and Na<-W coordination cutoff
    double r_aw = 0.38;  // Cl-O and Cl<-W coordination cutoff
    double r_oo = 0.35;  // O-O cutoff for water-water coordination
    std::string out_csv = "coord_x_mean.csv";
};

struct RunningStats1D {
    std::vector<double> sum;
    std::vector<double> sumsq;
    std::vector<int64_t> n_nonempty;

    void init(int n) {
        sum.assign((size_t)n, 0.0);
        sumsq.assign((size_t)n, 0.0);
        n_nonempty.assign((size_t)n, 0);
    }

    void add(int i, double v) {
        sum[(size_t)i] += v;
        sumsq[(size_t)i] += v * v;
        n_nonempty[(size_t)i] += 1;
    }
};

struct IntervalMap {
    bool inside = false;
    double u = 0.0;
    double length = 0.0;
};

struct FrameCoord {
    std::vector<double> sum_count[(size_t)MetricN];
    std::vector<int64_t> n_central[(size_t)MetricN];
    int64_t selected_water = 0;
    int64_t selected_na = 0;
    int64_t selected_cl = 0;
    int64_t binned_water = 0;
    int64_t binned_na = 0;
    int64_t binned_cl = 0;
    int64_t oob = 0;
    int64_t cand_visits_ibc = 0;
    int64_t cand_visits_iba = 0;
    int64_t cand_visits_bnc = 0;
    int64_t cand_visits_bna = 0;
    int64_t cand_visits_bww = 0;
    int max_ibc = 0;
    int max_iba = 0;
    int max_bnc = 0;
    int max_bna = 0;
    int max_bww = 0;
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
        << " [xmin=7.11] [xmax=12.89] [zmin=0.901] [zmax=1.801] [nx=100] [r_cw=0.35] [r_aw=0.38]"
        << " [r_oo=0.35]"
        << " [out_csv=coord_x_mean.csv]\n";
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
    if (argc > 13) c.r_cw = parse_double_arg(argv[13], "r_cw");
    if (argc > 14) c.r_aw = parse_double_arg(argv[14], "r_aw");
    if (argc > 15) c.r_oo = parse_double_arg(argv[15], "r_oo");
    if (argc > 16) c.out_csv = argv[16];

    if (c.max_frames <= 0 || c.nthreads <= 0 || c.grid_cell_nm <= 0.0 || c.nsol < 0 || c.nna < 0 ||
        c.ncl < 0 || c.nx <= 0 || c.r_cw <= 0.0 || c.r_aw <= 0.0 || c.r_oo <= 0.0) {
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

inline int metric_idx(Metric m) { return (int)m; }

inline bool in_type(simio::MolType got, simio::MolType want) { return got == want; }

FrameCoord compute_frame_coord(const CliConfig& cfg, const simio::Topology& topo, const simio::Frame& fr,
                               const std::vector<simio::MolState>& ms, int frame_id) {
    FrameCoord out;
    for (int m = 0; m < MetricN; ++m) {
        out.sum_count[(size_t)m].assign((size_t)cfg.nx, 0.0);
        out.n_central[(size_t)m].assign((size_t)cfg.nx, 0);
    }

    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    const double xminw = fr.pbc.wrap_pos(0, cfg.xmin);
    const double xmaxw = fr.pbc.wrap_pos(0, cfg.xmax);
    const double zminw = fr.pbc.wrap_pos(2, cfg.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg.zmax);

    const double xlen = interval_length_pbc(xminw, xmaxw, Lx);
    const double dx = xlen / (double)cfg.nx;
    if (dx <= 0.0) throw std::runtime_error("Invalid x bin width.");

    const int nx_cw = std::max(1, (int)std::ceil((cfg.r_cw + 0.05) / fr.grid.cx));
    const int ny_cw = std::max(1, (int)std::ceil((cfg.r_cw + 0.05) / fr.grid.cy));
    const int nz_cw = std::max(1, (int)std::ceil((cfg.r_cw + 0.05) / fr.grid.cz));
    const int nx_aw = std::max(1, (int)std::ceil((cfg.r_aw + 0.05) / fr.grid.cx));
    const int ny_aw = std::max(1, (int)std::ceil((cfg.r_aw + 0.05) / fr.grid.cy));
    const int nz_aw = std::max(1, (int)std::ceil((cfg.r_aw + 0.05) / fr.grid.cz));
    const int nx_oo = std::max(1, (int)std::ceil((cfg.r_oo + 0.05) / fr.grid.cx));
    const int ny_oo = std::max(1, (int)std::ceil((cfg.r_oo + 0.05) / fr.grid.cy));
    const int nz_oo = std::max(1, (int)std::ceil((cfg.r_oo + 0.05) / fr.grid.cz));

    const double r2_cw = cfg.r_cw * cfg.r_cw;
    const double r2_aw = cfg.r_aw * cfg.r_aw;
    const double r2_oo = cfg.r_oo * cfg.r_oo;
    constexpr double kCosThetaHb = 0.8660254037844386;  // cos(30 deg)
    constexpr double kRhcl = 0.30;
    const double r2_hcl = kRhcl * kRhcl;
    constexpr double kNormEps = 1e-12;

    std::vector<simio::Vec3d> pos((size_t)topo.mols.size(), simio::Vec3d{});
    std::vector<uint8_t> has_water_sites((size_t)topo.mols.size(), 0);
    std::vector<simio::Vec3d> water_h1u((size_t)topo.mols.size(), simio::Vec3d{});
    std::vector<simio::Vec3d> water_h2u((size_t)topo.mols.size(), simio::Vec3d{});
    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        const simio::MolSpan& m = topo.mols[mid];
        const simio::MolCache& c = ms[mid].cache;
        if (m.natoms < 1) continue;
        if (m.type == simio::MolType::Water && m.natoms >= 3) {
            if (c.flags & simio::MolCache::HAS_REF) {
                // Coordination/HB observables are O-centric for water.
                pos[mid] = c.ref_wrapped;
            } else {
                simio::Vec3d rO = fr.atoms.pos(m.first + 0);
                fr.pbc.wrap_pos3(rO);
                pos[mid] = rO;
            }
            if (c.flags & simio::MolCache::HAS_SITES) {
                has_water_sites[mid] = 1;
                water_h1u[mid] = c.sites_u[1];
                water_h2u[mid] = c.sites_u[2];
            }
            continue;
        }

        if (c.flags & simio::MolCache::HAS_KEY) {
            pos[mid] = c.key_wrapped;
            continue;
        }

        simio::Vec3d r = fr.atoms.pos(m.first);
        fr.pbc.wrap_pos3(r);
        pos[mid] = r;
    }

    int probe_water = -1;
    int probe_na = -1;
    int probe_cl = -1;
    if (!topo.mol_ids_by_type[(int)simio::MolType::Water].empty()) {
        probe_water = topo.mol_ids_by_type[(int)simio::MolType::Water][0];
    }
    if (!topo.mol_ids_by_type[(int)simio::MolType::Cation].empty()) {
        probe_na = topo.mol_ids_by_type[(int)simio::MolType::Cation][0];
    }
    if (!topo.mol_ids_by_type[(int)simio::MolType::Anion].empty()) {
        probe_cl = topo.mol_ids_by_type[(int)simio::MolType::Anion][0];
    }

    auto xbin_from = [&](const simio::Vec3d& r, int& ix) -> bool {
        const IntervalMap xm = map_on_pbc_interval(r.v[0], xminw, xmaxw, Lx);
        const IntervalMap zm = map_on_pbc_interval(r.v[2], zminw, zmaxw, Lz);
        if (!(xm.inside && zm.inside)) return false;
        ix = (int)std::floor(xm.u / dx);
        if (ix < 0 || ix >= cfg.nx) return false;
        return true;
    };

    for (int mid = 0; mid < (int)topo.mols.size(); ++mid) {
        const simio::MolSpan& m = topo.mols[(size_t)mid];
        const simio::Vec3d center = pos[(size_t)mid];
        int ix = -1;
        const bool in_roi = xbin_from(center, ix);
        if (!in_roi) continue;

        if (m.type == simio::MolType::Water) {
            ++out.selected_water;
            if (ix < 0 || ix >= cfg.nx) {
                ++out.oob;
                continue;
            }
            ++out.binned_water;

            simio::Vec3d rH1u{};
            simio::Vec3d rH2u{};
            if (has_water_sites[(size_t)mid]) {
                rH1u = water_h1u[(size_t)mid];
                rH2u = water_h2u[(size_t)mid];
            } else {
                const int H1 = m.first + 1;
                const int H2 = m.first + 2;
                const simio::Vec3d rH1 = fr.atoms.pos(H1);
                const simio::Vec3d rH2 = fr.atoms.pos(H2);
                rH1u = simio::pbc_unwrap_to_ref(fr.pbc, center, rH1);
                rH2u = simio::pbc_unwrap_to_ref(fr.pbc, center, rH2);
            }
            const simio::Vec3d OH1 = rH1u - center;
            const simio::Vec3d OH2 = rH2u - center;
            const double oh1_norm = simio::norm3(OH1);
            const double oh2_norm = simio::norm3(OH2);

            int ibc_count = 0;
            int iba_count = 0;
            int bww_count = 0;
            int hbww_don_count = 0;
            int hbww_acc_count = 0;
            int hbwcl_don_count = 0;
            double dmin_c = std::numeric_limits<double>::infinity();
            double dmin_a = std::numeric_limits<double>::infinity();
            int64_t cand_c = 0;
            int64_t cand_a = 0;
            int64_t cand_w2 = 0;
            bool h1_don_ww = false;
            bool h2_don_ww = false;
            bool h1_don_cl = false;
            bool h2_don_cl = false;

            fr.grid.for_candidates_box(center, nx_cw, ny_cw, nz_cw, [&](int cand_id) {
                if (!in_type(topo.mols[(size_t)cand_id].type, simio::MolType::Cation)) return;
                ++cand_c;
                const simio::Vec3d dr = fr.pbc.min_image(pos[(size_t)cand_id] - center);
                const double d2 = simio::dot3(dr, dr);
                if (d2 < dmin_c) dmin_c = d2;
                if (d2 <= r2_cw) ++ibc_count;
            });

            fr.grid.for_candidates_box(center, nx_aw, ny_aw, nz_aw, [&](int cand_id) {
                if (!in_type(topo.mols[(size_t)cand_id].type, simio::MolType::Anion)) return;
                ++cand_a;
                const simio::Vec3d dr = fr.pbc.min_image(pos[(size_t)cand_id] - center);
                const double d2 = simio::dot3(dr, dr);
                if (d2 < dmin_a) dmin_a = d2;
                if (d2 <= r2_aw) ++iba_count;

                if (d2 > r2_aw) return;
                const double ocl_norm = std::sqrt(d2);
                if (ocl_norm <= kNormEps) return;

                if (!h1_don_cl && oh1_norm > kNormEps) {
                    const double cosang = simio::dot3(OH1, dr) / (oh1_norm * ocl_norm);
                    const simio::Vec3d hcl = fr.pbc.min_image(pos[(size_t)cand_id] - rH1u);
                    if (cosang >= kCosThetaHb && simio::dot3(hcl, hcl) <= r2_hcl) {
                        h1_don_cl = true;
                        ++hbwcl_don_count;
                    }
                }
                if (!h2_don_cl && oh2_norm > kNormEps) {
                    const double cosang = simio::dot3(OH2, dr) / (oh2_norm * ocl_norm);
                    const simio::Vec3d hcl = fr.pbc.min_image(pos[(size_t)cand_id] - rH2u);
                    if (cosang >= kCosThetaHb && simio::dot3(hcl, hcl) <= r2_hcl) {
                        h2_don_cl = true;
                        ++hbwcl_don_count;
                    }
                }
            });

            fr.grid.for_candidates_box(center, nx_oo, ny_oo, nz_oo, [&](int cand_id) {
                if (!in_type(topo.mols[(size_t)cand_id].type, simio::MolType::Water)) return;
                if (cand_id == mid) return;
                ++cand_w2;
                const simio::Vec3d dr = fr.pbc.min_image(pos[(size_t)cand_id] - center);
                const double d2 = simio::dot3(dr, dr);
                if (d2 <= r2_oo) ++bww_count;
                if (d2 > r2_oo) return;

                const double oo_norm = std::sqrt(d2);
                if (oo_norm <= kNormEps) return;

                if (!h1_don_ww && oh1_norm > kNormEps) {
                    const double cosang = simio::dot3(OH1, dr) / (oh1_norm * oo_norm);
                    if (cosang >= kCosThetaHb) {
                        h1_don_ww = true;
                        ++hbww_don_count;
                    }
                }
                if (!h2_don_ww && oh2_norm > kNormEps) {
                    const double cosang = simio::dot3(OH2, dr) / (oh2_norm * oo_norm);
                    if (cosang >= kCosThetaHb) {
                        h2_don_ww = true;
                        ++hbww_don_count;
                    }
                }

                if (hbww_acc_count >= 2) return;
                const simio::MolSpan& wn = topo.mols[(size_t)cand_id];
                if (wn.natoms < 3) return;

                simio::Vec3d nH1u{};
                simio::Vec3d nH2u{};
                if (has_water_sites[(size_t)cand_id]) {
                    nH1u = water_h1u[(size_t)cand_id];
                    nH2u = water_h2u[(size_t)cand_id];
                } else {
                    const simio::Vec3d nH1 = fr.atoms.pos(wn.first + 1);
                    const simio::Vec3d nH2 = fr.atoms.pos(wn.first + 2);
                    nH1u = simio::pbc_unwrap_to_ref(fr.pbc, pos[(size_t)cand_id], nH1);
                    nH2u = simio::pbc_unwrap_to_ref(fr.pbc, pos[(size_t)cand_id], nH2);
                }
                const simio::Vec3d OtoCenter = fr.pbc.min_image(center - pos[(size_t)cand_id]);
                const double oc_norm = simio::norm3(OtoCenter);
                if (oc_norm <= kNormEps) return;

                bool donated_to_center = false;
                const simio::Vec3d nOH1 = nH1u - pos[(size_t)cand_id];
                const simio::Vec3d nOH2 = nH2u - pos[(size_t)cand_id];
                const double noh1_norm = simio::norm3(nOH1);
                const double noh2_norm = simio::norm3(nOH2);
                if (noh1_norm > kNormEps) {
                    const double cosang = simio::dot3(nOH1, OtoCenter) / (noh1_norm * oc_norm);
                    if (cosang >= kCosThetaHb) donated_to_center = true;
                }
                if (!donated_to_center && noh2_norm > kNormEps) {
                    const double cosang = simio::dot3(nOH2, OtoCenter) / (noh2_norm * oc_norm);
                    if (cosang >= kCosThetaHb) donated_to_center = true;
                }
                if (donated_to_center) ++hbww_acc_count;
            });

            const int hbww_tot_count = hbww_don_count + hbww_acc_count;

            out.sum_count[(size_t)metric_idx(IBC)][(size_t)ix] += (double)ibc_count;
            out.sum_count[(size_t)metric_idx(IBA)][(size_t)ix] += (double)iba_count;
            out.sum_count[(size_t)metric_idx(BWW)][(size_t)ix] += (double)bww_count;
            out.sum_count[(size_t)metric_idx(HBWW_DON)][(size_t)ix] += (double)hbww_don_count;
            out.sum_count[(size_t)metric_idx(HBWW_ACC)][(size_t)ix] += (double)hbww_acc_count;
            out.sum_count[(size_t)metric_idx(HBWW_TOT)][(size_t)ix] += (double)hbww_tot_count;
            out.sum_count[(size_t)metric_idx(HBWCL_DON)][(size_t)ix] += (double)hbwcl_don_count;
            out.n_central[(size_t)metric_idx(IBC)][(size_t)ix] += 1;
            out.n_central[(size_t)metric_idx(IBA)][(size_t)ix] += 1;
            out.n_central[(size_t)metric_idx(BWW)][(size_t)ix] += 1;
            out.n_central[(size_t)metric_idx(HBWW_DON)][(size_t)ix] += 1;
            out.n_central[(size_t)metric_idx(HBWW_ACC)][(size_t)ix] += 1;
            out.n_central[(size_t)metric_idx(HBWW_TOT)][(size_t)ix] += 1;
            out.n_central[(size_t)metric_idx(HBWCL_DON)][(size_t)ix] += 1;
            out.cand_visits_ibc += cand_c;
            out.cand_visits_iba += cand_a;
            out.cand_visits_bww += cand_w2;
            if (ibc_count > out.max_ibc) out.max_ibc = ibc_count;
            if (iba_count > out.max_iba) out.max_iba = iba_count;
            if (bww_count > out.max_bww) out.max_bww = bww_count;

            if (frame_id < 3 && mid == probe_water) {
                const double dmin_c_nm = std::isfinite(dmin_c) ? std::sqrt(dmin_c) : -1.0;
                const double dmin_a_nm = std::isfinite(dmin_a) ? std::sqrt(dmin_a) : -1.0;
                std::cout << "[coord dbg][frame " << frame_id << "][water " << mid << "] cand(na/cl/w)="
                          << cand_c << "/" << cand_a << "/" << cand_w2 << " IBC=" << ibc_count
                          << " IBA=" << iba_count << " BWW=" << bww_count << " HBww(d/a/t)="
                          << hbww_don_count << "/" << hbww_acc_count << "/" << hbww_tot_count
                          << " HBwcl_don=" << hbwcl_don_count
                          << " dmin(na/cl)=" << dmin_c_nm << "/" << dmin_a_nm << "\n";
            }
        } else if (m.type == simio::MolType::Cation) {
            ++out.selected_na;
            if (ix < 0 || ix >= cfg.nx) {
                ++out.oob;
                continue;
            }
            ++out.binned_na;

            int bnc_count = 0;
            double dmin_w = std::numeric_limits<double>::infinity();
            int64_t cand_w = 0;
            fr.grid.for_candidates_box(center, nx_cw, ny_cw, nz_cw, [&](int cand_id) {
                if (!in_type(topo.mols[(size_t)cand_id].type, simio::MolType::Water)) return;
                ++cand_w;
                const simio::Vec3d dr = fr.pbc.min_image(pos[(size_t)cand_id] - center);
                const double d2 = simio::dot3(dr, dr);
                if (d2 < dmin_w) dmin_w = d2;
                if (d2 <= r2_cw) ++bnc_count;
            });
            out.sum_count[(size_t)metric_idx(BNC)][(size_t)ix] += (double)bnc_count;
            out.n_central[(size_t)metric_idx(BNC)][(size_t)ix] += 1;
            out.cand_visits_bnc += cand_w;
            if (bnc_count > out.max_bnc) out.max_bnc = bnc_count;

            if (frame_id < 3 && mid == probe_na) {
                const double dmin_nm = std::isfinite(dmin_w) ? std::sqrt(dmin_w) : -1.0;
                std::cout << "[coord dbg][frame " << frame_id << "][na " << mid << "] cand(w)=" << cand_w
                          << " BNC=" << bnc_count << " dmin(w)=" << dmin_nm << "\n";
            }
        } else if (m.type == simio::MolType::Anion) {
            ++out.selected_cl;
            if (ix < 0 || ix >= cfg.nx) {
                ++out.oob;
                continue;
            }
            ++out.binned_cl;

            int bna_count = 0;
            double dmin_w = std::numeric_limits<double>::infinity();
            int64_t cand_w = 0;
            fr.grid.for_candidates_box(center, nx_aw, ny_aw, nz_aw, [&](int cand_id) {
                if (!in_type(topo.mols[(size_t)cand_id].type, simio::MolType::Water)) return;
                ++cand_w;
                const simio::Vec3d dr = fr.pbc.min_image(pos[(size_t)cand_id] - center);
                const double d2 = simio::dot3(dr, dr);
                if (d2 < dmin_w) dmin_w = d2;
                if (d2 <= r2_aw) ++bna_count;
            });
            out.sum_count[(size_t)metric_idx(BNA)][(size_t)ix] += (double)bna_count;
            out.n_central[(size_t)metric_idx(BNA)][(size_t)ix] += 1;
            out.cand_visits_bna += cand_w;
            if (bna_count > out.max_bna) out.max_bna = bna_count;

            if (frame_id < 3 && mid == probe_cl) {
                const double dmin_nm = std::isfinite(dmin_w) ? std::sqrt(dmin_w) : -1.0;
                std::cout << "[coord dbg][frame " << frame_id << "][cl " << mid << "] cand(w)=" << cand_w
                          << " BNA=" << bna_count << " dmin(w)=" << dmin_nm << "\n";
            }
        }
    }

    return out;
}

double mean_from_nonempty(const RunningStats1D& st, int i) {
    const int64_t n = st.n_nonempty[(size_t)i];
    if (n <= 0) return 0.0;
    return st.sum[(size_t)i] / (double)n;
}

double sem_from_nonempty(const RunningStats1D& st, int i) {
    const int64_t n = st.n_nonempty[(size_t)i];
    if (n <= 1) return 0.0;
    const double sum = st.sum[(size_t)i];
    const double sumsq = st.sumsq[(size_t)i];
    const double mean = sum / (double)n;
    double var = (sumsq - (double)n * mean * mean) / (double)(n - 1);
    if (var < 0.0) var = 0.0;
    return std::sqrt(var / (double)n);
}

void write_csv(const std::string& out_csv, const std::vector<double>& x_centers, const RunningStats1D& ibc,
               const RunningStats1D& iba, const RunningStats1D& bnc, const RunningStats1D& bna,
               const RunningStats1D& bww, const RunningStats1D& hbww_don, const RunningStats1D& hbww_acc,
               const RunningStats1D& hbww_tot, const RunningStats1D& hbwcl_don) {
    std::ofstream ofs(out_csv);
    if (!ofs) throw std::runtime_error("Failed to open output CSV: " + out_csv);

    ofs << std::setprecision(12);
    ofs << "x_center_nm,IBC_mean,IBC_sem,IBA_mean,IBA_sem,BNC_mean,BNC_sem,BNA_mean,BNA_sem,BWW_mean,BWW_sem,"
           "HBww_don_mean,HBww_don_sem,HBww_acc_mean,HBww_acc_sem,HBww_tot_mean,HBww_tot_sem,"
           "HBwcl_don_mean,HBwcl_don_sem\n";
    for (size_t i = 0; i < x_centers.size(); ++i) {
        const double ibc_m = mean_from_nonempty(ibc, (int)i);
        const double iba_m = mean_from_nonempty(iba, (int)i);
        const double bnc_m = mean_from_nonempty(bnc, (int)i);
        const double bna_m = mean_from_nonempty(bna, (int)i);
        const double bww_m = mean_from_nonempty(bww, (int)i);
        const double hbww_don_m = mean_from_nonempty(hbww_don, (int)i);
        const double hbww_acc_m = mean_from_nonempty(hbww_acc, (int)i);
        const double hbww_tot_m = mean_from_nonempty(hbww_tot, (int)i);
        const double hbwcl_don_m = mean_from_nonempty(hbwcl_don, (int)i);
        ofs << x_centers[i] << "," << ibc_m << "," << sem_from_nonempty(ibc, (int)i) << "," << iba_m << ","
            << sem_from_nonempty(iba, (int)i) << "," << bnc_m << "," << sem_from_nonempty(bnc, (int)i)
            << "," << bna_m << "," << sem_from_nonempty(bna, (int)i) << "," << bww_m << ","
            << sem_from_nonempty(bww, (int)i) << "," << hbww_don_m << ","
            << sem_from_nonempty(hbww_don, (int)i) << "," << hbww_acc_m << ","
            << sem_from_nonempty(hbww_acc, (int)i) << "," << hbww_tot_m << ","
            << sem_from_nonempty(hbww_tot, (int)i) << "," << hbwcl_don_m << ","
            << sem_from_nonempty(hbwcl_don, (int)i) << "\n";
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

        RunningStats1D st_ibc, st_iba, st_bnc, st_bna, st_bww;
        RunningStats1D st_hbww_don, st_hbww_acc, st_hbww_tot, st_hbwcl_don;
        st_ibc.init(cfg.nx);
        st_iba.init(cfg.nx);
        st_bnc.init(cfg.nx);
        st_bna.init(cfg.nx);
        st_bww.init(cfg.nx);
        st_hbww_don.init(cfg.nx);
        st_hbww_acc.init(cfg.nx);
        st_hbww_tot.init(cfg.nx);
        st_hbwcl_don.init(cfg.nx);

        std::vector<double> x_centers((size_t)cfg.nx, 0.0);
        bool has_ref_xgrid = false;
        bool warned_xgrid_change = false;
        double ref_xlen = 0.0;
        double ref_xminw = 0.0;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[coord-config] file=" << cfg.xtc_path << " frames=" << cfg.max_frames
                  << " threads=" << cfg.nthreads << " nx=" << cfg.nx << " x=[" << cfg.xmin << "," << cfg.xmax
                  << "] z=[" << cfg.zmin << "," << cfg.zmax << "] r_cw=" << cfg.r_cw << " r_aw=" << cfg.r_aw
                  << " r_oo=" << cfg.r_oo << " out_csv=" << cfg.out_csv << "\n";
        std::cout << "[topology] nmol=" << expected_nmols << " water=" << cfg.nsol << " na=" << cfg.nna
                  << " cl=" << cfg.ncl << " natoms=" << expected_natoms << "\n";
        std::cout << "[coord-def] distance atom=water oxygen only\n";

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

            const double xminw = fr.pbc.wrap_pos(0, cfg.xmin);
            const double xmaxw = fr.pbc.wrap_pos(0, cfg.xmax);
            const double xlen = interval_length_pbc(xminw, xmaxw, fr.pbc.L[0]);
            if (xlen <= 0.0) {
                xtc_close(&traj);
                throw std::runtime_error("Invalid x interval length.");
            }

            if (!has_ref_xgrid) {
                has_ref_xgrid = true;
                ref_xlen = xlen;
                ref_xminw = xminw;
                const double dx = ref_xlen / (double)cfg.nx;
                for (int i = 0; i < cfg.nx; ++i) {
                    x_centers[(size_t)i] = fr.pbc.wrap_pos(0, ref_xminw + ((double)i + 0.5) * dx);
                }
            } else if (!warned_xgrid_change && std::abs(xlen - ref_xlen) > 1e-6) {
                std::cout << "[warn] x interval length changed from " << ref_xlen << " to " << xlen
                          << "; CSV x centers use first-frame mapping\n";
                warned_xgrid_change = true;
            }

            const FrameCoord fc = compute_frame_coord(cfg, topo, fr, ms, frames_done);
            const bool ok = (fc.selected_water == fc.binned_water) && (fc.selected_na == fc.binned_na) &&
                            (fc.selected_cl == fc.binned_cl) && (fc.oob == 0);

            for (int i = 0; i < cfg.nx; ++i) {
                auto add_if_nonempty = [&](RunningStats1D& st, Metric m) {
                    const int mi = metric_idx(m);
                    const int64_t n = fc.n_central[(size_t)mi][(size_t)i];
                    if (n <= 0) return;
                    const double frame_mean = fc.sum_count[(size_t)mi][(size_t)i] / (double)n;
                    st.add(i, frame_mean);
                };
                add_if_nonempty(st_ibc, IBC);
                add_if_nonempty(st_iba, IBA);
                add_if_nonempty(st_bnc, BNC);
                add_if_nonempty(st_bna, BNA);
                add_if_nonempty(st_bww, BWW);
                add_if_nonempty(st_hbww_don, HBWW_DON);
                add_if_nonempty(st_hbww_acc, HBWW_ACC);
                add_if_nonempty(st_hbww_tot, HBWW_TOT);
                add_if_nonempty(st_hbwcl_don, HBWCL_DON);
            }

            std::cout << "[coord frame " << frames_done << "] step=" << fr.step << " time_ps=" << fr.time_ps
                      << " selected(w/na/cl)=" << fc.selected_water << "/" << fc.selected_na << "/"
                      << fc.selected_cl << " binned(w/na/cl)=" << fc.binned_water << "/" << fc.binned_na << "/"
                      << fc.binned_cl << " oob=" << fc.oob << " status=" << (ok ? "OK" : "FAIL") << "\n";
            if (frames_done < 3) {
                std::cout << "  candidates(IBC/IBA/BNC/BNA/BWW)=" << fc.cand_visits_ibc << "/"
                          << fc.cand_visits_iba << "/" << fc.cand_visits_bnc << "/" << fc.cand_visits_bna
                          << "/" << fc.cand_visits_bww << "\n";
            }
            constexpr int kSuspiciousNeighborCount = 50;
            if (fc.max_ibc >= kSuspiciousNeighborCount || fc.max_iba >= kSuspiciousNeighborCount ||
                fc.max_bnc >= kSuspiciousNeighborCount || fc.max_bna >= kSuspiciousNeighborCount ||
                fc.max_bww >= kSuspiciousNeighborCount) {
                std::cout << "[warn] frame " << frames_done
                          << " unusually high coordination count observed: max(IBC/IBA/BNC/BNA/BWW)="
                          << fc.max_ibc << "/" << fc.max_iba << "/" << fc.max_bnc << "/" << fc.max_bna << "/"
                          << fc.max_bww << " (heuristic threshold=" << kSuspiciousNeighborCount << ")\n";
            }
            if (!ok) {
                xtc_close(&traj);
                throw std::runtime_error("Coordination invariant failed.");
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

        write_csv(cfg.out_csv, x_centers, st_ibc, st_iba, st_bnc, st_bna, st_bww, st_hbww_don, st_hbww_acc,
                  st_hbww_tot, st_hbwcl_don);
        std::cout << "Processed " << frames_done << " frame(s). Wrote " << cfg.out_csv << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
