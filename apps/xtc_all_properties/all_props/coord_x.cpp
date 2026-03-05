#include "coord_x.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace simio::analysis {

CoordXAnalyzer::CoordXAnalyzer(const CoordXConfig& cfg) : cfg_(cfg) {
    if (cfg_.nx <= 0) throw std::runtime_error("CoordXAnalyzer: nx must be > 0");
    if (cfg_.r_cw <= 0.0 || cfg_.r_aw <= 0.0 || cfg_.r_oo <= 0.0) {
        throw std::runtime_error("CoordXAnalyzer: cutoffs must be > 0");
    }
    x_centers_.assign(static_cast<size_t>(cfg_.nx), 0.0);
    for (int m = 0; m < MetricN; ++m) stats_[static_cast<size_t>(m)].init(cfg_.nx);
}

void CoordXAnalyzer::process_frame(const Topology& topo, const Frame& fr,
                                   const std::vector<MolState>& ms) {
    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Lz <= 0.0) throw std::runtime_error("CoordXAnalyzer: invalid box lengths");

    const bool whole_x = (cfg_.xmax <= cfg_.xmin);
    const double xminw = whole_x ? 0.0 : fr.pbc.wrap_pos(0, cfg_.xmin);
    const double xmaxw = whole_x ? 0.0 : fr.pbc.wrap_pos(0, cfg_.xmax);
    const double zminw = fr.pbc.wrap_pos(2, cfg_.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg_.zmax);

    const double xlen = whole_x ? Lx : interval_length_pbc(xminw, xmaxw, Lx);
    const double dx = xlen / static_cast<double>(cfg_.nx);
    if (xlen <= 0.0 || dx <= 0.0) throw std::runtime_error("CoordXAnalyzer: invalid x interval");

    if (!has_x_centers_) {
        has_x_centers_ = true;
        const double x0 = whole_x ? 0.0 : xminw;
        for (int i = 0; i < cfg_.nx; ++i) {
            x_centers_[static_cast<size_t>(i)] =
                fr.pbc.wrap_pos(0, x0 + (static_cast<double>(i) + 0.5) * dx);
        }
    }

    const int nx_cw = std::max(1, static_cast<int>(std::ceil((cfg_.r_cw + 0.05) / fr.grid.cx)));
    const int ny_cw = std::max(1, static_cast<int>(std::ceil((cfg_.r_cw + 0.05) / fr.grid.cy)));
    const int nz_cw = std::max(1, static_cast<int>(std::ceil((cfg_.r_cw + 0.05) / fr.grid.cz)));
    const int nx_aw = std::max(1, static_cast<int>(std::ceil((cfg_.r_aw + 0.05) / fr.grid.cx)));
    const int ny_aw = std::max(1, static_cast<int>(std::ceil((cfg_.r_aw + 0.05) / fr.grid.cy)));
    const int nz_aw = std::max(1, static_cast<int>(std::ceil((cfg_.r_aw + 0.05) / fr.grid.cz)));
    const int nx_oo = std::max(1, static_cast<int>(std::ceil((cfg_.r_oo + 0.05) / fr.grid.cx)));
    const int ny_oo = std::max(1, static_cast<int>(std::ceil((cfg_.r_oo + 0.05) / fr.grid.cy)));
    const int nz_oo = std::max(1, static_cast<int>(std::ceil((cfg_.r_oo + 0.05) / fr.grid.cz)));

    const double r2_cw = cfg_.r_cw * cfg_.r_cw;
    const double r2_aw = cfg_.r_aw * cfg_.r_aw;
    const double r2_oo = cfg_.r_oo * cfg_.r_oo;
    constexpr double kCosThetaHb = 0.8660254037844386;  // cos(30 deg)
    constexpr double kRhcl = 0.30;
    const double r2_hcl = kRhcl * kRhcl;
    constexpr double kNormEps = 1e-12;

    const size_t nmol = topo.mols.size();
    std::vector<Vec3d> pos(nmol, Vec3d{});
    std::vector<uint8_t> has_water_sites(nmol, 0);
    std::vector<Vec3d> water_h1u(nmol, Vec3d{});
    std::vector<Vec3d> water_h2u(nmol, Vec3d{});

    for (size_t mid = 0; mid < nmol; ++mid) {
        const MolSpan& m = topo.mols[mid];
        const MolCache& c = ms[mid].cache;
        if (m.natoms < 1) continue;

        if (m.type == MolType::Water && m.natoms >= 3) {
            if (c.flags & MolCache::HAS_REF) {
                pos[mid] = c.ref_wrapped;
            } else {
                Vec3d rO = fr.atoms.pos(m.first + 0);
                fr.pbc.wrap_pos3(rO);
                pos[mid] = rO;
            }

            if (c.flags & MolCache::HAS_SITES) {
                has_water_sites[mid] = 1;
                water_h1u[mid] = c.sites_u[1];
                water_h2u[mid] = c.sites_u[2];
            }
            continue;
        }

        if (c.flags & MolCache::HAS_KEY) {
            pos[mid] = c.key_wrapped;
            continue;
        }

        Vec3d r = fr.atoms.pos(m.first);
        fr.pbc.wrap_pos3(r);
        pos[mid] = r;
    }

    auto xbin_from = [&](const Vec3d& r, int& ix) -> bool {
        IntervalMap xm{};
        if (whole_x) {
            xm.inside = true;
            xm.u = r.v[0];
            xm.length = Lx;
        } else {
            xm = map_on_pbc_interval(r.v[0], xminw, xmaxw, Lx);
        }
        const IntervalMap zm = map_on_pbc_interval(r.v[2], zminw, zmaxw, Lz);
        if (!(xm.inside && zm.inside)) return false;
        ix = static_cast<int>(std::floor(xm.u / dx));
        if (ix < 0 || ix >= cfg_.nx) return false;
        return true;
    };

    std::array<std::vector<double>, MetricN> sum_count{};
    std::array<std::vector<int64_t>, MetricN> n_central{};
    for (int m = 0; m < MetricN; ++m) {
        sum_count[static_cast<size_t>(m)].assign(static_cast<size_t>(cfg_.nx), 0.0);
        n_central[static_cast<size_t>(m)].assign(static_cast<size_t>(cfg_.nx), 0);
    }

    auto in_type = [](MolType got, MolType want) { return got == want; };

    for (int mid = 0; mid < static_cast<int>(nmol); ++mid) {
        const MolSpan& m = topo.mols[static_cast<size_t>(mid)];
        const Vec3d center = pos[static_cast<size_t>(mid)];

        int ix = -1;
        if (!xbin_from(center, ix)) continue;

        if (m.type == MolType::Water) {
            Vec3d rH1u{};
            Vec3d rH2u{};
            if (has_water_sites[static_cast<size_t>(mid)]) {
                rH1u = water_h1u[static_cast<size_t>(mid)];
                rH2u = water_h2u[static_cast<size_t>(mid)];
            } else {
                const Vec3d rH1 = fr.atoms.pos(m.first + 1);
                const Vec3d rH2 = fr.atoms.pos(m.first + 2);
                rH1u = pbc_unwrap_to_ref(fr.pbc, center, rH1);
                rH2u = pbc_unwrap_to_ref(fr.pbc, center, rH2);
            }
            const Vec3d OH1 = rH1u - center;
            const Vec3d OH2 = rH2u - center;
            const double oh1_norm = norm3(OH1);
            const double oh2_norm = norm3(OH2);

            int ibc_count = 0;
            int iba_count = 0;
            int bww_count = 0;
            int hbww_don_count = 0;
            int hbww_acc_count = 0;
            int hbwcl_don_count = 0;

            bool h1_don_ww = false;
            bool h2_don_ww = false;
            bool h1_don_cl = false;
            bool h2_don_cl = false;

            fr.grid.for_candidates_box(center, nx_cw, ny_cw, nz_cw, [&](int cand_id) {
                if (!in_type(topo.mols[static_cast<size_t>(cand_id)].type, MolType::Cation)) return;
                const Vec3d dr = fr.pbc.min_image(pos[static_cast<size_t>(cand_id)] - center);
                const double d2 = dot3(dr, dr);
                if (d2 <= r2_cw) ++ibc_count;
            });

            fr.grid.for_candidates_box(center, nx_aw, ny_aw, nz_aw, [&](int cand_id) {
                if (!in_type(topo.mols[static_cast<size_t>(cand_id)].type, MolType::Anion)) return;
                const Vec3d dr = fr.pbc.min_image(pos[static_cast<size_t>(cand_id)] - center);
                const double d2 = dot3(dr, dr);
                if (d2 <= r2_aw) ++iba_count;

                if (d2 > r2_aw) return;
                const double ocl_norm = std::sqrt(d2);
                if (ocl_norm <= kNormEps) return;

                if (!h1_don_cl && oh1_norm > kNormEps) {
                    const double cosang = dot3(OH1, dr) / (oh1_norm * ocl_norm);
                    const Vec3d hcl = fr.pbc.min_image(pos[static_cast<size_t>(cand_id)] - rH1u);
                    if (cosang >= kCosThetaHb && dot3(hcl, hcl) <= r2_hcl) {
                        h1_don_cl = true;
                        ++hbwcl_don_count;
                    }
                }
                if (!h2_don_cl && oh2_norm > kNormEps) {
                    const double cosang = dot3(OH2, dr) / (oh2_norm * ocl_norm);
                    const Vec3d hcl = fr.pbc.min_image(pos[static_cast<size_t>(cand_id)] - rH2u);
                    if (cosang >= kCosThetaHb && dot3(hcl, hcl) <= r2_hcl) {
                        h2_don_cl = true;
                        ++hbwcl_don_count;
                    }
                }
            });

            fr.grid.for_candidates_box(center, nx_oo, ny_oo, nz_oo, [&](int cand_id) {
                if (!in_type(topo.mols[static_cast<size_t>(cand_id)].type, MolType::Water)) return;
                if (cand_id == mid) return;

                const Vec3d dr = fr.pbc.min_image(pos[static_cast<size_t>(cand_id)] - center);
                const double d2 = dot3(dr, dr);
                if (d2 <= r2_oo) ++bww_count;
                if (d2 > r2_oo) return;

                const double oo_norm = std::sqrt(d2);
                if (oo_norm <= kNormEps) return;

                if (!h1_don_ww && oh1_norm > kNormEps) {
                    const double cosang = dot3(OH1, dr) / (oh1_norm * oo_norm);
                    if (cosang >= kCosThetaHb) {
                        h1_don_ww = true;
                        ++hbww_don_count;
                    }
                }
                if (!h2_don_ww && oh2_norm > kNormEps) {
                    const double cosang = dot3(OH2, dr) / (oh2_norm * oo_norm);
                    if (cosang >= kCosThetaHb) {
                        h2_don_ww = true;
                        ++hbww_don_count;
                    }
                }

                if (hbww_acc_count >= 2) return;

                const MolSpan& wn = topo.mols[static_cast<size_t>(cand_id)];
                if (wn.natoms < 3) return;

                Vec3d nH1u{};
                Vec3d nH2u{};
                if (has_water_sites[static_cast<size_t>(cand_id)]) {
                    nH1u = water_h1u[static_cast<size_t>(cand_id)];
                    nH2u = water_h2u[static_cast<size_t>(cand_id)];
                } else {
                    const Vec3d nH1 = fr.atoms.pos(wn.first + 1);
                    const Vec3d nH2 = fr.atoms.pos(wn.first + 2);
                    nH1u = pbc_unwrap_to_ref(fr.pbc, pos[static_cast<size_t>(cand_id)], nH1);
                    nH2u = pbc_unwrap_to_ref(fr.pbc, pos[static_cast<size_t>(cand_id)], nH2);
                }

                const Vec3d OtoCenter = fr.pbc.min_image(center - pos[static_cast<size_t>(cand_id)]);
                const double oc_norm = norm3(OtoCenter);
                if (oc_norm <= kNormEps) return;

                bool donated_to_center = false;
                const Vec3d nOH1 = nH1u - pos[static_cast<size_t>(cand_id)];
                const Vec3d nOH2 = nH2u - pos[static_cast<size_t>(cand_id)];
                const double noh1_norm = norm3(nOH1);
                const double noh2_norm = norm3(nOH2);

                if (noh1_norm > kNormEps) {
                    const double cosang = dot3(nOH1, OtoCenter) / (noh1_norm * oc_norm);
                    if (cosang >= kCosThetaHb) donated_to_center = true;
                }
                if (!donated_to_center && noh2_norm > kNormEps) {
                    const double cosang = dot3(nOH2, OtoCenter) / (noh2_norm * oc_norm);
                    if (cosang >= kCosThetaHb) donated_to_center = true;
                }
                if (donated_to_center) ++hbww_acc_count;
            });

            const int hbww_tot_count = hbww_don_count + hbww_acc_count;

            sum_count[IBC][static_cast<size_t>(ix)] += static_cast<double>(ibc_count);
            sum_count[IBA][static_cast<size_t>(ix)] += static_cast<double>(iba_count);
            sum_count[BWW][static_cast<size_t>(ix)] += static_cast<double>(bww_count);
            sum_count[HBWW_DON][static_cast<size_t>(ix)] += static_cast<double>(hbww_don_count);
            sum_count[HBWW_ACC][static_cast<size_t>(ix)] += static_cast<double>(hbww_acc_count);
            sum_count[HBWW_TOT][static_cast<size_t>(ix)] += static_cast<double>(hbww_tot_count);
            sum_count[HBWCL_DON][static_cast<size_t>(ix)] += static_cast<double>(hbwcl_don_count);

            n_central[IBC][static_cast<size_t>(ix)] += 1;
            n_central[IBA][static_cast<size_t>(ix)] += 1;
            n_central[BWW][static_cast<size_t>(ix)] += 1;
            n_central[HBWW_DON][static_cast<size_t>(ix)] += 1;
            n_central[HBWW_ACC][static_cast<size_t>(ix)] += 1;
            n_central[HBWW_TOT][static_cast<size_t>(ix)] += 1;
            n_central[HBWCL_DON][static_cast<size_t>(ix)] += 1;

        } else if (m.type == MolType::Cation) {
            int bnc_count = 0;
            fr.grid.for_candidates_box(center, nx_cw, ny_cw, nz_cw, [&](int cand_id) {
                if (!in_type(topo.mols[static_cast<size_t>(cand_id)].type, MolType::Water)) return;
                const Vec3d dr = fr.pbc.min_image(pos[static_cast<size_t>(cand_id)] - center);
                if (dot3(dr, dr) <= r2_cw) ++bnc_count;
            });

            sum_count[BNC][static_cast<size_t>(ix)] += static_cast<double>(bnc_count);
            n_central[BNC][static_cast<size_t>(ix)] += 1;

        } else if (m.type == MolType::Anion) {
            int bna_count = 0;
            fr.grid.for_candidates_box(center, nx_aw, ny_aw, nz_aw, [&](int cand_id) {
                if (!in_type(topo.mols[static_cast<size_t>(cand_id)].type, MolType::Water)) return;
                const Vec3d dr = fr.pbc.min_image(pos[static_cast<size_t>(cand_id)] - center);
                if (dot3(dr, dr) <= r2_aw) ++bna_count;
            });

            sum_count[BNA][static_cast<size_t>(ix)] += static_cast<double>(bna_count);
            n_central[BNA][static_cast<size_t>(ix)] += 1;
        }
    }

    for (int m = 0; m < MetricN; ++m) {
        for (int i = 0; i < cfg_.nx; ++i) {
            const int64_t n = n_central[static_cast<size_t>(m)][static_cast<size_t>(i)];
            if (n <= 0) continue;
            const double v = sum_count[static_cast<size_t>(m)][static_cast<size_t>(i)] /
                             static_cast<double>(n);
            stats_[static_cast<size_t>(m)].add(i, v);
        }
    }

    nframes_ += 1;
}

void CoordXAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("CoordXAnalyzer: no frames processed");
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("CoordXAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "x_center_nm,IBC_mean,IBC_sem,IBA_mean,IBA_sem,BNC_mean,BNC_sem,BNA_mean,BNA_sem,BWW_mean,BWW_sem,"
           "HBww_don_mean,HBww_don_sem,HBww_acc_mean,HBww_acc_sem,HBww_tot_mean,HBww_tot_sem,"
           "HBwcl_don_mean,HBwcl_don_sem\n";

    for (size_t i = 0; i < x_centers_.size(); ++i) {
        const int ix = static_cast<int>(i);
        ofs << x_centers_[i] << "," << stats_[IBC].mean(ix) << "," << stats_[IBC].sem(ix) << ","
            << stats_[IBA].mean(ix) << "," << stats_[IBA].sem(ix) << "," << stats_[BNC].mean(ix)
            << "," << stats_[BNC].sem(ix) << "," << stats_[BNA].mean(ix) << ","
            << stats_[BNA].sem(ix) << "," << stats_[BWW].mean(ix) << "," << stats_[BWW].sem(ix)
            << "," << stats_[HBWW_DON].mean(ix) << "," << stats_[HBWW_DON].sem(ix) << ","
            << stats_[HBWW_ACC].mean(ix) << "," << stats_[HBWW_ACC].sem(ix) << ","
            << stats_[HBWW_TOT].mean(ix) << "," << stats_[HBWW_TOT].sem(ix) << ","
            << stats_[HBWCL_DON].mean(ix) << "," << stats_[HBWCL_DON].sem(ix) << "\n";
    }
}

}  // namespace simio::analysis
