#include "jump_msd.hpp"

#include "simio/analysis/intrinsics/channel_roi.hpp"
#include "simio/analysis/intrinsics/context.hpp"
#include "simio/analysis/intrinsics/in_channel_mask.hpp"
#include "simio/runtime/cache.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace simio::analysis {

JumpMsdAnalyzer::JumpMsdAnalyzer(const JumpMsdConfig& cfg) : cfg_(cfg) {
    if (cfg_.keep_frames <= 0) {
        throw std::runtime_error("JumpMsdAnalyzer: keep_frames must be > 0");
    }
    if (cfg_.bound_layer_nm < 0.0) {
        throw std::runtime_error("JumpMsdAnalyzer: bound_layer_nm must be >= 0");
    }

    const int nlag = cfg_.keep_frames + 1;
    n_samples_by_lag_.assign(static_cast<size_t>(nlag), 0);
    lag_time_ps_.init(nlag);

    msd_x_all_.init(nlag);
    msd_y_all_.init(nlag);
    msd_z_all_.init(nlag);
    msd_z_raw_all_.init(nlag);
    msd_r2_all_.init(nlag);

    max_jump_x_all_.init(nlag);

    for (int sid = 0; sid < 3; ++sid) {
        msd_x_sp_[static_cast<size_t>(sid)].init(nlag);
        msd_y_sp_[static_cast<size_t>(sid)].init(nlag);
        msd_z_sp_[static_cast<size_t>(sid)].init(nlag);
        msd_z_raw_sp_[static_cast<size_t>(sid)].init(nlag);
        msd_r2_sp_[static_cast<size_t>(sid)].init(nlag);
        max_jump_x_sp_[static_cast<size_t>(sid)].init(nlag);
        msd_x_channel_strict_sum_dx2_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag),
                                                                       0.0);
        msd_x_channel_strict_n_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0);
        msd_x_channel_loose_sum_dx2_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag),
                                                                      0.0);
        msd_x_channel_loose_n_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0);
        msd_z_channel_unwrapped_sum_dz2_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag),
                                                                          0.0);
        msd_z_channel_raw_sum_dz2_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0.0);
        msd_z_channel_n_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0);
        vacf_y_sum_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0.0);
        vacf_y_n_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0);
        vacf_x_channel_sum_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0.0);
        vacf_x_channel_n_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0);
        vacf_z_channel_raw_sum_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0.0);
        vacf_z_channel_raw_n_[static_cast<size_t>(sid)].assign(static_cast<size_t>(nlag), 0);
    }
    msd_x_channel_strict_lag_ps_sum_.assign(static_cast<size_t>(nlag), 0.0);
    msd_x_channel_strict_lag_ps_n_.assign(static_cast<size_t>(nlag), 0);
    msd_x_channel_loose_lag_ps_sum_.assign(static_cast<size_t>(nlag), 0.0);
    msd_x_channel_loose_lag_ps_n_.assign(static_cast<size_t>(nlag), 0);
    vacf_y_lag_ps_sum_.assign(static_cast<size_t>(nlag), 0.0);
    vacf_y_lag_ps_n_.assign(static_cast<size_t>(nlag), 0);
    vacf_x_channel_lag_ps_sum_.assign(static_cast<size_t>(nlag), 0.0);
    vacf_x_channel_lag_ps_n_.assign(static_cast<size_t>(nlag), 0);
    vacf_z_channel_raw_lag_ps_sum_.assign(static_cast<size_t>(nlag), 0.0);
    vacf_z_channel_raw_lag_ps_n_.assign(static_cast<size_t>(nlag), 0);
}

JumpMsdAnalyzer::JumpMsdAnalyzer(const JumpMsdConfig& cfg, simio::runtime::CacheStore& cache)
    : JumpMsdAnalyzer(cfg) {
    cache_ = &cache;
}

void JumpMsdAnalyzer::process_frame(const Topology& topo, const Frame& fr,
                                    const std::vector<MolState>& ms, int frame_idx) {
    if (nmol_ == 0) {
        nmol_ = topo.mols.size();
        species_idx_.assign(nmol_, -1);
        runlen_in_channel_.assign(nmol_, 0);
        for (size_t mid = 0; mid < nmol_; ++mid) {
            species_idx_[mid] = species_index_from_type(topo.mols[mid].type);
        }
    }

    if (topo.mols.size() != nmol_ || ms.size() != nmol_) {
        throw std::runtime_error("JumpMsdAnalyzer: topology/state size mismatch");
    }

    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Lz <= 0.0) throw std::runtime_error("JumpMsdAnalyzer: invalid box lengths");

    if (!has_roi_) {
        if (!cache_) {
            static thread_local simio::runtime::CacheStore fallback_cache;
            cache_ = &fallback_cache;
        }
        simio::analysis::intrinsics::IntrinsicContext ictx{*cache_};
        roi_ = simio::analysis::intrinsics::get_channel_roi_x(
            ictx, cfg_.x_channel_min, cfg_.x_channel_max, Lx);
        xlen_ = roi_.xlen;
        has_roi_ = true;
    }

    const double xlen = xlen_;
    const bool whole_x_channel = (xlen <= 0.0);

    const double zminw = fr.pbc.wrap_pos(2, cfg_.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg_.zmax);
    const double zlen = interval_length_pbc(zminw, zmaxw, Lz);
    if (zlen <= 0.0) throw std::runtime_error("JumpMsdAnalyzer: invalid z slab interval");

    const double bound_delta_eff = std::min(cfg_.bound_layer_nm, 0.5 * zlen);

    Snapshot cur{};
    cur.frame_idx = frame_idx;
    cur.step = fr.step;
    cur.time_ps = fr.time_ps;
    cur.key_cont.assign(nmol_, Vec3d{});
    cur.x_raw.assign(nmol_, 0.0);
    cur.z_raw.assign(nmol_, 0.0);
    cur.vx_raw.assign(nmol_, 0.0);
    cur.vy.assign(nmol_, 0.0);
    cur.vz_raw.assign(nmol_, 0.0);
    cur.valid.assign(nmol_, 0);
    cur.vel_valid.assign(nmol_, 0);
    cur.in_slab.assign(nmol_, 0);
    cur.in_channel.assign(nmol_, 0);
    cur.in_bound.assign(nmol_, 0);
    cur.in_core.assign(nmol_, 0);

    if (xw_tmp_.size() != nmol_) xw_tmp_.assign(nmol_, 0.0);
    for (size_t mid = 0; mid < nmol_; ++mid) {
        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) {
            xw_tmp_[mid] = 0.0;
            continue;
        }
        Vec3d keyw = c.key_wrapped;
        fr.pbc.wrap_pos3(keyw);
        xw_tmp_[mid] = keyw.v[0];
    }

    const std::int64_t frame_id = frame_counter_++;
    simio::analysis::intrinsics::IntrinsicContext ictx{*cache_};
    const auto mask = simio::analysis::intrinsics::get_in_channel_mask_x(
        ictx,
        frame_id,
        xw_tmp_.data(),
        xw_tmp_.size(),
        cfg_.x_channel_min,
        cfg_.x_channel_max,
        Lx);
    if (mask.in_channel.size() != nmol_) {
        throw std::runtime_error("JumpMsdAnalyzer: in_channel_mask size mismatch");
    }

    for (size_t mid = 0; mid < nmol_; ++mid) {
        const MolSpan& mol = topo.mols[mid];
        const MolState& st = ms[mid];
        bool is_valid = false;
        bool is_in_slab = false;
        bool is_in_channel = false;
        bool is_in_bound = false;
        bool is_in_core = false;

        if (mol.natoms > 0) {
            // Raw per-frame x from the trajectory.
            cur.x_raw[mid] = fr.atoms.x[static_cast<size_t>(mol.first)];
            // Raw per-frame z from the trajectory (no PBC correction/unwrapping).
            cur.z_raw[mid] = fr.atoms.z[static_cast<size_t>(mol.first)];
        }

        if ((st.cache.flags & MolCache::HAS_KEY) && st.time.has_prev) {
            cur.key_cont[mid] = st.time.key_cont;
            is_valid = true;

            Vec3d keyw = st.cache.key_wrapped;
            fr.pbc.wrap_pos3(keyw);
            const IntervalMap z_map = map_on_pbc_interval(keyw.v[2], zminw, zmaxw, Lz);
            is_in_slab = z_map.inside;
            const double x_u = roi_.map_x_to_channel(keyw.v[0]);
            const bool x_ok =
                whole_x_channel ? true
                                : (mask.in_channel[mid] != 0u && x_u >= 0.0 && x_u < xlen);
            is_in_channel = x_ok && is_in_slab;
            if (is_in_channel) {
                if (bound_delta_eff > 0.0) {
                    is_in_bound =
                        (z_map.u < bound_delta_eff) || (z_map.u > (zlen - bound_delta_eff));
                }
                is_in_core = !is_in_bound;
            }
        }

        cur.valid[mid] = is_valid ? 1 : 0;
        cur.in_slab[mid] = is_in_slab ? 1 : 0;
        cur.in_channel[mid] = is_in_channel ? 1 : 0;
        cur.in_bound[mid] = is_in_bound ? 1 : 0;
        cur.in_core[mid] = is_in_core ? 1 : 0;
        runlen_in_channel_[mid] = (is_valid && is_in_channel) ? (runlen_in_channel_[mid] + 1) : 0;
    }

    StateFrameRow state_row{};
    state_row.frame_idx = frame_idx;
    state_row.step = fr.step;
    state_row.time_ps = fr.time_ps;
    for (size_t mid = 0; mid < nmol_; ++mid) {
        if (!(cur.valid[mid] && cur.in_channel[mid])) continue;
        const int sid = species_idx_[mid];
        if (sid < 0 || sid >= 3) continue;
        const size_t s = static_cast<size_t>(sid);
        if (cur.in_bound[mid]) {
            state_row.n_bound[s] += 1;
        } else if (cur.in_core[mid]) {
            state_row.n_core[s] += 1;
        }
    }
    state_rows_.push_back(state_row);

    if (!history_.empty()) {
        const Snapshot& prev = history_.back();
        const double dt_step_ps = cur.time_ps - prev.time_ps;
        if (dt_step_ps > 0.0) {
            for (size_t mid = 0; mid < nmol_; ++mid) {
                if (!(cur.valid[mid] && prev.valid[mid])) continue;
                cur.vx_raw[mid] = (cur.x_raw[mid] - prev.x_raw[mid]) / dt_step_ps;
                cur.vy[mid] = (cur.key_cont[mid].v[1] - prev.key_cont[mid].v[1]) / dt_step_ps;
                cur.vz_raw[mid] = (cur.z_raw[mid] - prev.z_raw[mid]) / dt_step_ps;
                cur.vel_valid[mid] = 1;

                const int sid = species_idx_[mid];
                if (sid < 0 || sid >= 3) continue;
                const size_t s = static_cast<size_t>(sid);
                if (prev.in_core[mid] && cur.in_core[mid]) {
                    drift_vx_core_[s].add(cur.vx_raw[mid]);
                    drift_dt_ps_sum_ += dt_step_ps;
                    drift_dt_ps_n_ += 1;
                }
                if (prev.in_bound[mid] && cur.in_bound[mid]) {
                    drift_vx_bound_[s].add(cur.vx_raw[mid]);
                    drift_dt_ps_sum_ += dt_step_ps;
                    drift_dt_ps_n_ += 1;
                }
            }
        }
    }

    int64_t n_vacf_lag0_any = 0;
    for (size_t mid = 0; mid < nmol_; ++mid) {
        if (!(cur.vel_valid[mid] && cur.in_slab[mid])) continue;
        const int sid = species_idx_[mid];
        if (sid < 0 || sid >= 3) continue;
        const size_t s = static_cast<size_t>(sid);
        vacf_y_sum_[s][0] += cur.vy[mid] * cur.vy[mid];
        vacf_y_n_[s][0] += 1;
        n_vacf_lag0_any += 1;
    }
    if (n_vacf_lag0_any > 0) {
        vacf_y_lag_ps_n_[0] += n_vacf_lag0_any;
    }

    int64_t n_vacf_x_channel_lag0_any = 0;
    for (size_t mid = 0; mid < nmol_; ++mid) {
        if (!(cur.vel_valid[mid] && cur.in_channel[mid])) continue;
        const int sid = species_idx_[mid];
        if (sid < 0 || sid >= 3) continue;
        const size_t s = static_cast<size_t>(sid);
        vacf_x_channel_sum_[s][0] += cur.vx_raw[mid] * cur.vx_raw[mid];
        vacf_x_channel_n_[s][0] += 1;
        n_vacf_x_channel_lag0_any += 1;
    }
    if (n_vacf_x_channel_lag0_any > 0) {
        vacf_x_channel_lag_ps_n_[0] += n_vacf_x_channel_lag0_any;
    }

    int64_t n_vacf_z_channel_lag0_any = 0;
    for (size_t mid = 0; mid < nmol_; ++mid) {
        if (!(cur.vel_valid[mid] && cur.in_channel[mid])) continue;
        const int sid = species_idx_[mid];
        if (sid < 0 || sid >= 3) continue;
        const size_t s = static_cast<size_t>(sid);
        vacf_z_channel_raw_sum_[s][0] += cur.vz_raw[mid] * cur.vz_raw[mid];
        vacf_z_channel_raw_n_[s][0] += 1;
        n_vacf_z_channel_lag0_any += 1;
    }
    if (n_vacf_z_channel_lag0_any > 0) {
        vacf_z_channel_raw_lag_ps_n_[0] += n_vacf_z_channel_lag0_any;
    }

    for (const Snapshot& old : history_) {
        const int lag = cur.frame_idx - old.frame_idx;
        if (lag <= 0 || lag > cfg_.keep_frames) continue;
        const double lag_time_ps = cur.time_ps - old.time_ps;

        double sum_dx2_all = 0.0;
        double sum_dy2_all = 0.0;
        double sum_dz2_all = 0.0;
        double sum_dz2_raw_all = 0.0;
        int64_t n_all = 0;
        double max_jump_x_all = 0.0;
        bool has_all = false;

        std::array<double, 3> sum_dx2_sp{0.0, 0.0, 0.0};
        std::array<double, 3> sum_dy2_sp{0.0, 0.0, 0.0};
        std::array<double, 3> sum_dz2_sp{0.0, 0.0, 0.0};
        std::array<double, 3> sum_dz2_raw_sp{0.0, 0.0, 0.0};
        std::array<int64_t, 3> n_sp{0, 0, 0};
        std::array<double, 3> max_jump_x_sp{0.0, 0.0, 0.0};
        std::array<uint8_t, 3> has_sp{0, 0, 0};
        int64_t n_channel_strict_samples_any = 0;
        int64_t n_channel_loose_samples_any = 0;
        int64_t n_vacf_any = 0;
        int64_t n_vacf_x_channel_any = 0;
        int64_t n_vacf_z_channel_any = 0;

        for (size_t mid = 0; mid < nmol_; ++mid) {
            if (!(cur.valid[mid] && old.valid[mid])) continue;
            if (!old.in_slab[mid]) continue;

            const Vec3d dr = cur.key_cont[mid] - old.key_cont[mid];
            const double dx2 = dr.v[0] * dr.v[0];
            const double dy2 = dr.v[1] * dr.v[1];
            const double dz2 = dr.v[2] * dr.v[2];
            const double dz_raw = cur.z_raw[mid] - old.z_raw[mid];
            const double dz2_raw = dz_raw * dz_raw;
            const double dx_mi = fr.pbc.wrap_delta(0, dr.v[0]);
            const double jump_x = std::abs(dx_mi);

            sum_dx2_all += dx2;
            sum_dy2_all += dy2;
            sum_dz2_all += dz2;
            sum_dz2_raw_all += dz2_raw;
            n_all += 1;
            if (!has_all || jump_x > max_jump_x_all) {
                has_all = true;
                max_jump_x_all = jump_x;
            }

            const int sid = species_idx_[mid];
            if (sid >= 0 && sid < 3) {
                const size_t s = static_cast<size_t>(sid);
                sum_dx2_sp[s] += dx2;
                sum_dy2_sp[s] += dy2;
                sum_dz2_sp[s] += dz2;
                sum_dz2_raw_sp[s] += dz2_raw;
                n_sp[s] += 1;
                if (!has_sp[s] || jump_x > max_jump_x_sp[s]) {
                    has_sp[s] = 1;
                    max_jump_x_sp[s] = jump_x;
                }

                if (cur.vel_valid[mid] && old.vel_valid[mid]) {
                    vacf_y_sum_[s][static_cast<size_t>(lag)] += old.vy[mid] * cur.vy[mid];
                    vacf_y_n_[s][static_cast<size_t>(lag)] += 1;
                    n_vacf_any += 1;
                }

                // Loose channel MSD: origin must be in channel; endpoint can be outside.
                if (old.in_channel[mid]) {
                    const double dx_channel_loose = cur.key_cont[mid].v[0] - old.key_cont[mid].v[0];
                    msd_x_channel_loose_sum_dx2_[s][static_cast<size_t>(lag)] +=
                        dx_channel_loose * dx_channel_loose;
                    msd_x_channel_loose_n_[s][static_cast<size_t>(lag)] += 1;
                    n_channel_loose_samples_any += 1;
                }

                if (runlen_in_channel_[mid] >= lag + 1) {
                    if (cur.vel_valid[mid] && old.vel_valid[mid]) {
                        vacf_x_channel_sum_[s][static_cast<size_t>(lag)] +=
                            old.vx_raw[mid] * cur.vx_raw[mid];
                        vacf_x_channel_n_[s][static_cast<size_t>(lag)] += 1;
                        n_vacf_x_channel_any += 1;
                    }

                    if (cur.vel_valid[mid] && old.vel_valid[mid]) {
                        vacf_z_channel_raw_sum_[s][static_cast<size_t>(lag)] +=
                            old.vz_raw[mid] * cur.vz_raw[mid];
                        vacf_z_channel_raw_n_[s][static_cast<size_t>(lag)] += 1;
                        n_vacf_z_channel_any += 1;
                    }

                    const double dx_channel = cur.key_cont[mid].v[0] - old.key_cont[mid].v[0];
                    msd_x_channel_strict_sum_dx2_[s][static_cast<size_t>(lag)] +=
                        dx_channel * dx_channel;
                    msd_x_channel_strict_n_[s][static_cast<size_t>(lag)] += 1;

                    const double dz_unwrapped_channel =
                        cur.key_cont[mid].v[2] - old.key_cont[mid].v[2];
                    const double dz_raw_channel = cur.z_raw[mid] - old.z_raw[mid];
                    msd_z_channel_unwrapped_sum_dz2_[s][static_cast<size_t>(lag)] +=
                        dz_unwrapped_channel * dz_unwrapped_channel;
                    msd_z_channel_raw_sum_dz2_[s][static_cast<size_t>(lag)] +=
                        dz_raw_channel * dz_raw_channel;
                    msd_z_channel_n_[s][static_cast<size_t>(lag)] += 1;

                    if (lag == 1) {
                        const double abs_dz = std::abs(dz_raw_channel);
                        if (abs_dz > max_abs_dz_lag1_sp_[s]) max_abs_dz_lag1_sp_[s] = abs_dz;
                    }
                    n_channel_strict_samples_any += 1;
                }
            }
        }

        if (n_all > 0) {
            lag_time_ps_.add(lag, lag_time_ps);
            const double inv = 1.0 / static_cast<double>(n_all);
            msd_x_all_.add(lag, sum_dx2_all * inv);
            msd_y_all_.add(lag, sum_dy2_all * inv);
            msd_z_all_.add(lag, sum_dz2_all * inv);
            msd_z_raw_all_.add(lag, sum_dz2_raw_all * inv);
            msd_r2_all_.add(lag, (sum_dx2_all + sum_dy2_all + sum_dz2_all) * inv);
            max_jump_x_all_.add(lag, max_jump_x_all);
            n_samples_by_lag_[static_cast<size_t>(lag)] += n_all;
        }

        for (int sid = 0; sid < 3; ++sid) {
            const size_t s = static_cast<size_t>(sid);
            if (n_sp[s] <= 0) continue;

            const double inv = 1.0 / static_cast<double>(n_sp[s]);
            msd_x_sp_[s].add(lag, sum_dx2_sp[s] * inv);
            msd_y_sp_[s].add(lag, sum_dy2_sp[s] * inv);
            msd_z_sp_[s].add(lag, sum_dz2_sp[s] * inv);
            msd_z_raw_sp_[s].add(lag, sum_dz2_raw_sp[s] * inv);
            msd_r2_sp_[s].add(lag, (sum_dx2_sp[s] + sum_dy2_sp[s] + sum_dz2_sp[s]) * inv);
            max_jump_x_sp_[s].add(lag, max_jump_x_sp[s]);
        }

        if (n_channel_strict_samples_any > 0) {
            msd_x_channel_strict_lag_ps_sum_[static_cast<size_t>(lag)] +=
                lag_time_ps * static_cast<double>(n_channel_strict_samples_any);
            msd_x_channel_strict_lag_ps_n_[static_cast<size_t>(lag)] +=
                n_channel_strict_samples_any;
        }

        if (n_channel_loose_samples_any > 0) {
            msd_x_channel_loose_lag_ps_sum_[static_cast<size_t>(lag)] +=
                lag_time_ps * static_cast<double>(n_channel_loose_samples_any);
            msd_x_channel_loose_lag_ps_n_[static_cast<size_t>(lag)] +=
                n_channel_loose_samples_any;
        }

        if (n_vacf_any > 0) {
            vacf_y_lag_ps_sum_[static_cast<size_t>(lag)] +=
                lag_time_ps * static_cast<double>(n_vacf_any);
            vacf_y_lag_ps_n_[static_cast<size_t>(lag)] += n_vacf_any;
        }

        if (n_vacf_x_channel_any > 0) {
            vacf_x_channel_lag_ps_sum_[static_cast<size_t>(lag)] +=
                lag_time_ps * static_cast<double>(n_vacf_x_channel_any);
            vacf_x_channel_lag_ps_n_[static_cast<size_t>(lag)] += n_vacf_x_channel_any;
        }

        if (n_vacf_z_channel_any > 0) {
            vacf_z_channel_raw_lag_ps_sum_[static_cast<size_t>(lag)] +=
                lag_time_ps * static_cast<double>(n_vacf_z_channel_any);
            vacf_z_channel_raw_lag_ps_n_[static_cast<size_t>(lag)] += n_vacf_z_channel_any;
        }
    }

    history_.push_back(std::move(cur));
    const size_t max_history = static_cast<size_t>(cfg_.keep_frames + 1);
    while (history_.size() > max_history) history_.pop_front();

    nframes_ += 1;
}

void JumpMsdAnalyzer::write_channel_msd_x_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("JumpMsdAnalyzer: no frames processed");

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("JumpMsdAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "lag_frames,lag_ps_strict,lag_ps_loose,"
           "msd_x_channel_water_nm2_strict_mean,n_samples_water_strict,msd_x_channel_water_nm2_loose_mean,n_samples_water_loose,"
           "msd_x_channel_na_nm2_strict_mean,n_samples_na_strict,msd_x_channel_na_nm2_loose_mean,n_samples_na_loose,"
           "msd_x_channel_cl_nm2_strict_mean,n_samples_cl_strict,msd_x_channel_cl_nm2_loose_mean,n_samples_cl_loose\n";

    for (int lag = 1; lag <= cfg_.keep_frames; ++lag) {
        const size_t i = static_cast<size_t>(lag);
        const int64_t nlag_strict = msd_x_channel_strict_lag_ps_n_[i];
        const int64_t nlag_loose = msd_x_channel_loose_lag_ps_n_[i];
        const double lag_ps_strict = (nlag_strict > 0)
                                         ? (msd_x_channel_strict_lag_ps_sum_[i] /
                                            static_cast<double>(nlag_strict))
                                         : 0.0;
        const double lag_ps_loose = (nlag_loose > 0)
                                        ? (msd_x_channel_loose_lag_ps_sum_[i] /
                                           static_cast<double>(nlag_loose))
                                        : 0.0;

        const int64_t nw_strict = msd_x_channel_strict_n_[0][i];
        const int64_t nna_strict = msd_x_channel_strict_n_[1][i];
        const int64_t ncl_strict = msd_x_channel_strict_n_[2][i];
        const int64_t nw_loose = msd_x_channel_loose_n_[0][i];
        const int64_t nna_loose = msd_x_channel_loose_n_[1][i];
        const int64_t ncl_loose = msd_x_channel_loose_n_[2][i];

        const double msd_w_strict = (nw_strict > 0)
                                        ? (msd_x_channel_strict_sum_dx2_[0][i] /
                                           static_cast<double>(nw_strict))
                                        : 0.0;
        const double msd_na_strict = (nna_strict > 0)
                                         ? (msd_x_channel_strict_sum_dx2_[1][i] /
                                            static_cast<double>(nna_strict))
                                         : 0.0;
        const double msd_cl_strict = (ncl_strict > 0)
                                         ? (msd_x_channel_strict_sum_dx2_[2][i] /
                                            static_cast<double>(ncl_strict))
                                         : 0.0;
        const double msd_w_loose = (nw_loose > 0)
                                       ? (msd_x_channel_loose_sum_dx2_[0][i] /
                                          static_cast<double>(nw_loose))
                                       : 0.0;
        const double msd_na_loose = (nna_loose > 0)
                                        ? (msd_x_channel_loose_sum_dx2_[1][i] /
                                           static_cast<double>(nna_loose))
                                        : 0.0;
        const double msd_cl_loose = (ncl_loose > 0)
                                        ? (msd_x_channel_loose_sum_dx2_[2][i] /
                                           static_cast<double>(ncl_loose))
                                        : 0.0;

        ofs << lag << "," << lag_ps_strict << "," << lag_ps_loose << "," << msd_w_strict << ","
            << nw_strict << "," << msd_w_loose << "," << nw_loose << "," << msd_na_strict << ","
            << nna_strict << "," << msd_na_loose << "," << nna_loose << "," << msd_cl_strict
            << "," << ncl_strict << "," << msd_cl_loose << "," << ncl_loose << "\n";
    }
}

void JumpMsdAnalyzer::write_channel_msd_z_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("JumpMsdAnalyzer: no frames processed");

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("JumpMsdAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "lag_frames,lag_ps,"
           "msd_z_channel_unwrapped_water_nm2_mean,msd_z_channel_raw_water_nm2_mean,n_samples_water,"
           "msd_z_channel_unwrapped_na_nm2_mean,msd_z_channel_raw_na_nm2_mean,n_samples_na,"
           "msd_z_channel_unwrapped_cl_nm2_mean,msd_z_channel_raw_cl_nm2_mean,n_samples_cl,"
           "max_abs_dz_lag1_water_nm,max_abs_dz_lag1_na_nm,max_abs_dz_lag1_cl_nm\n";

    for (int lag = 1; lag <= cfg_.keep_frames; ++lag) {
        const size_t i = static_cast<size_t>(lag);
        const int64_t nlag = msd_x_channel_strict_lag_ps_n_[i];
        const double lag_ps =
            (nlag > 0) ? (msd_x_channel_strict_lag_ps_sum_[i] / static_cast<double>(nlag)) : 0.0;

        const int64_t nw = msd_z_channel_n_[0][i];
        const int64_t nna = msd_z_channel_n_[1][i];
        const int64_t ncl = msd_z_channel_n_[2][i];

        const double msd_w_u = (nw > 0)
                                   ? (msd_z_channel_unwrapped_sum_dz2_[0][i] /
                                      static_cast<double>(nw))
                                   : 0.0;
        const double msd_w_r =
            (nw > 0) ? (msd_z_channel_raw_sum_dz2_[0][i] / static_cast<double>(nw)) : 0.0;

        const double msd_na_u = (nna > 0)
                                    ? (msd_z_channel_unwrapped_sum_dz2_[1][i] /
                                       static_cast<double>(nna))
                                    : 0.0;
        const double msd_na_r =
            (nna > 0) ? (msd_z_channel_raw_sum_dz2_[1][i] / static_cast<double>(nna)) : 0.0;

        const double msd_cl_u = (ncl > 0)
                                    ? (msd_z_channel_unwrapped_sum_dz2_[2][i] /
                                       static_cast<double>(ncl))
                                    : 0.0;
        const double msd_cl_r =
            (ncl > 0) ? (msd_z_channel_raw_sum_dz2_[2][i] / static_cast<double>(ncl)) : 0.0;

        ofs << lag << "," << lag_ps << "," << msd_w_u << "," << msd_w_r << "," << nw << ","
            << msd_na_u << "," << msd_na_r << "," << nna << "," << msd_cl_u << "," << msd_cl_r
            << "," << ncl << "," << max_abs_dz_lag1_sp_[0] << "," << max_abs_dz_lag1_sp_[1]
            << "," << max_abs_dz_lag1_sp_[2] << "\n";
    }
}

void JumpMsdAnalyzer::write_state_channel_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("JumpMsdAnalyzer: no frames processed");

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("JumpMsdAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "frame_idx,step,time_ps,"
           "n_water_bound,n_water_core,n_na_bound,n_na_core,n_cl_bound,n_cl_core\n";
    for (const auto& r : state_rows_) {
        ofs << r.frame_idx << "," << r.step << "," << r.time_ps << "," << r.n_bound[0] << ","
            << r.n_core[0] << "," << r.n_bound[1] << "," << r.n_core[1] << "," << r.n_bound[2]
            << "," << r.n_core[2] << "\n";
    }
}

void JumpMsdAnalyzer::write_drift_channel_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("JumpMsdAnalyzer: no frames processed");

    auto mean = [](const DriftAccum& a) -> double {
        if (a.n <= 0) return 0.0;
        return a.sum / static_cast<double>(a.n);
    };
    auto sem = [](const DriftAccum& a) -> double {
        if (a.n <= 1) return 0.0;
        const double m = a.sum / static_cast<double>(a.n);
        double var =
            (a.sumsq - static_cast<double>(a.n) * m * m) / static_cast<double>(a.n - 1);
        if (var < 0.0) var = 0.0;
        return std::sqrt(var / static_cast<double>(a.n));
    };

    const double dt_ps_mean =
        (drift_dt_ps_n_ > 0) ? (drift_dt_ps_sum_ / static_cast<double>(drift_dt_ps_n_)) : 0.0;

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("JumpMsdAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "dt_ps_mean,"
           "vx_water_core_nm_per_ps_mean,vx_water_core_nm_per_ps_sem,n_samples_water_core,"
           "vx_na_core_nm_per_ps_mean,vx_na_core_nm_per_ps_sem,n_samples_na_core,"
           "vx_cl_core_nm_per_ps_mean,vx_cl_core_nm_per_ps_sem,n_samples_cl_core,"
           "vx_water_bound_nm_per_ps_mean,vx_water_bound_nm_per_ps_sem,n_samples_water_bound,"
           "vx_na_bound_nm_per_ps_mean,vx_na_bound_nm_per_ps_sem,n_samples_na_bound,"
           "vx_cl_bound_nm_per_ps_mean,vx_cl_bound_nm_per_ps_sem,n_samples_cl_bound\n";

    ofs << dt_ps_mean << "," << mean(drift_vx_core_[0]) << "," << sem(drift_vx_core_[0]) << ","
        << drift_vx_core_[0].n << "," << mean(drift_vx_core_[1]) << ","
        << sem(drift_vx_core_[1]) << "," << drift_vx_core_[1].n << ","
        << mean(drift_vx_core_[2]) << "," << sem(drift_vx_core_[2]) << ","
        << drift_vx_core_[2].n << "," << mean(drift_vx_bound_[0]) << ","
        << sem(drift_vx_bound_[0]) << "," << drift_vx_bound_[0].n << ","
        << mean(drift_vx_bound_[1]) << "," << sem(drift_vx_bound_[1]) << ","
        << drift_vx_bound_[1].n << "," << mean(drift_vx_bound_[2]) << ","
        << sem(drift_vx_bound_[2]) << "," << drift_vx_bound_[2].n << "\n";
}

void JumpMsdAnalyzer::write_vacf_y_csv(const std::string& path,
                                       std::array<double, 3>* plateau_last10_nm2_per_ps) const {
    if (nframes_ <= 0) throw std::runtime_error("JumpMsdAnalyzer: no frames processed");

    const int nlag = cfg_.keep_frames + 1;
    std::vector<double> lag_ps(static_cast<size_t>(nlag), 0.0);
    for (int lag = 1; lag < nlag; ++lag) {
        const size_t i = static_cast<size_t>(lag);
        const int64_t n = vacf_y_lag_ps_n_[i];
        if (n > 0) lag_ps[i] = vacf_y_lag_ps_sum_[i] / static_cast<double>(n);
    }

    std::array<std::vector<double>, 3> vacf_mean{};
    std::array<std::vector<double>, 3> dy_from_vacf{};
    std::array<double, 3> plateau_nm2_per_ps{0.0, 0.0, 0.0};

    for (int sid = 0; sid < 3; ++sid) {
        const size_t s = static_cast<size_t>(sid);
        vacf_mean[s].assign(static_cast<size_t>(nlag), 0.0);
        dy_from_vacf[s].assign(static_cast<size_t>(nlag), 0.0);

        for (int lag = 0; lag < nlag; ++lag) {
            const size_t i = static_cast<size_t>(lag);
            const int64_t n = vacf_y_n_[s][i];
            if (n > 0) vacf_mean[s][i] = vacf_y_sum_[s][i] / static_cast<double>(n);
        }

        for (int lag = 1; lag < nlag; ++lag) {
            const size_t i = static_cast<size_t>(lag);
            const size_t ip = static_cast<size_t>(lag - 1);
            double dt = lag_ps[i] - lag_ps[ip];
            if (dt < 0.0) dt = 0.0;
            dy_from_vacf[s][i] =
                dy_from_vacf[s][ip] + 0.5 * (vacf_mean[s][ip] + vacf_mean[s][i]) * dt;
        }

        const int tail = (nlag < 10) ? nlag : 10;
        double sum_tail = 0.0;
        for (int k = 0; k < tail; ++k) {
            sum_tail += dy_from_vacf[s][static_cast<size_t>(nlag - tail + k)];
        }
        if (tail > 0) plateau_nm2_per_ps[s] = sum_tail / static_cast<double>(tail);
    }

    if (plateau_last10_nm2_per_ps) {
        *plateau_last10_nm2_per_ps = plateau_nm2_per_ps;
    }

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("JumpMsdAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "lag_frames,lag_ps,"
           "vacf_y_water_nm2_per_ps2,D_y_from_vacf_water_nm2_per_ps,n_samples_water,"
           "vacf_y_na_nm2_per_ps2,D_y_from_vacf_na_nm2_per_ps,n_samples_na,"
           "vacf_y_cl_nm2_per_ps2,D_y_from_vacf_cl_nm2_per_ps,n_samples_cl\n";

    for (int lag = 0; lag < nlag; ++lag) {
        const size_t i = static_cast<size_t>(lag);
        ofs << lag << "," << lag_ps[i] << "," << vacf_mean[0][i] << "," << dy_from_vacf[0][i]
            << "," << vacf_y_n_[0][i] << "," << vacf_mean[1][i] << "," << dy_from_vacf[1][i]
            << "," << vacf_y_n_[1][i] << "," << vacf_mean[2][i] << "," << dy_from_vacf[2][i]
            << "," << vacf_y_n_[2][i] << "\n";
    }
}

void JumpMsdAnalyzer::write_vacf_x_channel_csv(
    const std::string& path, std::array<double, 3>* plateau_last10_nm2_per_ps) const {
    if (nframes_ <= 0) throw std::runtime_error("JumpMsdAnalyzer: no frames processed");

    const int nlag = cfg_.keep_frames + 1;
    std::vector<double> lag_ps(static_cast<size_t>(nlag), 0.0);
    for (int lag = 1; lag < nlag; ++lag) {
        const size_t i = static_cast<size_t>(lag);
        const int64_t n = vacf_x_channel_lag_ps_n_[i];
        if (n > 0) lag_ps[i] = vacf_x_channel_lag_ps_sum_[i] / static_cast<double>(n);
    }

    std::array<std::vector<double>, 3> vacf_mean{};
    std::array<std::vector<double>, 3> int_from_vacf{};
    std::array<double, 3> plateau_nm2_per_ps{0.0, 0.0, 0.0};

    for (int sid = 0; sid < 3; ++sid) {
        const size_t s = static_cast<size_t>(sid);
        vacf_mean[s].assign(static_cast<size_t>(nlag), 0.0);
        int_from_vacf[s].assign(static_cast<size_t>(nlag), 0.0);

        for (int lag = 0; lag < nlag; ++lag) {
            const size_t i = static_cast<size_t>(lag);
            const int64_t n = vacf_x_channel_n_[s][i];
            if (n > 0) vacf_mean[s][i] = vacf_x_channel_sum_[s][i] / static_cast<double>(n);
        }

        for (int lag = 1; lag < nlag; ++lag) {
            const size_t i = static_cast<size_t>(lag);
            const size_t ip = static_cast<size_t>(lag - 1);
            double dt = lag_ps[i] - lag_ps[ip];
            if (dt < 0.0) dt = 0.0;
            int_from_vacf[s][i] =
                int_from_vacf[s][ip] + 0.5 * (vacf_mean[s][ip] + vacf_mean[s][i]) * dt;
        }

        const int tail = (nlag < 10) ? nlag : 10;
        double sum_tail = 0.0;
        for (int k = 0; k < tail; ++k) {
            sum_tail += int_from_vacf[s][static_cast<size_t>(nlag - tail + k)];
        }
        if (tail > 0) plateau_nm2_per_ps[s] = sum_tail / static_cast<double>(tail);
    }

    if (plateau_last10_nm2_per_ps) {
        *plateau_last10_nm2_per_ps = plateau_nm2_per_ps;
    }

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("JumpMsdAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "lag_frames,lag_ps,"
           "vacf_x_channel_water_nm2_per_ps2,I_x_from_vacf_channel_water_nm2_per_ps,n_samples_water,"
           "vacf_x_channel_na_nm2_per_ps2,I_x_from_vacf_channel_na_nm2_per_ps,n_samples_na,"
           "vacf_x_channel_cl_nm2_per_ps2,I_x_from_vacf_channel_cl_nm2_per_ps,n_samples_cl\n";

    for (int lag = 0; lag < nlag; ++lag) {
        const size_t i = static_cast<size_t>(lag);
        ofs << lag << "," << lag_ps[i] << "," << vacf_mean[0][i] << "," << int_from_vacf[0][i]
            << "," << vacf_x_channel_n_[0][i] << "," << vacf_mean[1][i] << ","
            << int_from_vacf[1][i] << "," << vacf_x_channel_n_[1][i] << ","
            << vacf_mean[2][i] << "," << int_from_vacf[2][i] << ","
            << vacf_x_channel_n_[2][i] << "\n";
    }
}

void JumpMsdAnalyzer::write_vacf_z_channel_raw_csv(
    const std::string& path, std::array<double, 3>* plateau_last10_nm2_per_ps) const {
    if (nframes_ <= 0) throw std::runtime_error("JumpMsdAnalyzer: no frames processed");

    const int nlag = cfg_.keep_frames + 1;
    std::vector<double> lag_ps(static_cast<size_t>(nlag), 0.0);
    for (int lag = 1; lag < nlag; ++lag) {
        const size_t i = static_cast<size_t>(lag);
        const int64_t n = vacf_z_channel_raw_lag_ps_n_[i];
        if (n > 0) lag_ps[i] = vacf_z_channel_raw_lag_ps_sum_[i] / static_cast<double>(n);
    }

    std::array<std::vector<double>, 3> vacf_mean{};
    std::array<std::vector<double>, 3> int_from_vacf{};
    std::array<double, 3> plateau_nm2_per_ps{0.0, 0.0, 0.0};

    for (int sid = 0; sid < 3; ++sid) {
        const size_t s = static_cast<size_t>(sid);
        vacf_mean[s].assign(static_cast<size_t>(nlag), 0.0);
        int_from_vacf[s].assign(static_cast<size_t>(nlag), 0.0);

        for (int lag = 0; lag < nlag; ++lag) {
            const size_t i = static_cast<size_t>(lag);
            const int64_t n = vacf_z_channel_raw_n_[s][i];
            if (n > 0) vacf_mean[s][i] = vacf_z_channel_raw_sum_[s][i] / static_cast<double>(n);
        }

        for (int lag = 1; lag < nlag; ++lag) {
            const size_t i = static_cast<size_t>(lag);
            const size_t ip = static_cast<size_t>(lag - 1);
            double dt = lag_ps[i] - lag_ps[ip];
            if (dt < 0.0) dt = 0.0;
            int_from_vacf[s][i] =
                int_from_vacf[s][ip] + 0.5 * (vacf_mean[s][ip] + vacf_mean[s][i]) * dt;
        }

        const int tail = (nlag < 10) ? nlag : 10;
        double sum_tail = 0.0;
        for (int k = 0; k < tail; ++k) {
            sum_tail += int_from_vacf[s][static_cast<size_t>(nlag - tail + k)];
        }
        if (tail > 0) plateau_nm2_per_ps[s] = sum_tail / static_cast<double>(tail);
    }

    if (plateau_last10_nm2_per_ps) {
        *plateau_last10_nm2_per_ps = plateau_nm2_per_ps;
    }

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("JumpMsdAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "lag_frames,lag_ps,"
           "vacf_z_channel_raw_water_nm2_per_ps2,I_z_from_vacf_channel_raw_water_nm2_per_ps,n_samples_water,"
           "vacf_z_channel_raw_na_nm2_per_ps2,I_z_from_vacf_channel_raw_na_nm2_per_ps,n_samples_na,"
           "vacf_z_channel_raw_cl_nm2_per_ps2,I_z_from_vacf_channel_raw_cl_nm2_per_ps,n_samples_cl\n";

    for (int lag = 0; lag < nlag; ++lag) {
        const size_t i = static_cast<size_t>(lag);
        ofs << lag << "," << lag_ps[i] << "," << vacf_mean[0][i] << "," << int_from_vacf[0][i]
            << "," << vacf_z_channel_raw_n_[0][i] << "," << vacf_mean[1][i] << ","
            << int_from_vacf[1][i] << "," << vacf_z_channel_raw_n_[1][i] << ","
            << vacf_mean[2][i] << "," << int_from_vacf[2][i] << ","
            << vacf_z_channel_raw_n_[2][i] << "\n";
    }
}

void JumpMsdAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("JumpMsdAnalyzer: no frames processed");

    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("JumpMsdAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "lag_frames,lag_time_ps_mean,lag_time_ps_sem,n_pairs,n_samples_all,"
           "msd_x_nm2_mean,msd_x_nm2_sem,msd_y_nm2_mean,msd_y_nm2_sem,msd_z_nm2_mean,msd_z_nm2_sem,msd_z_raw_nm2_mean,msd_z_raw_nm2_sem,msd_r2_nm2_mean,msd_r2_nm2_sem,"
           "msd_x_water_nm2_mean,msd_x_water_nm2_sem,msd_y_water_nm2_mean,msd_y_water_nm2_sem,msd_z_water_nm2_mean,msd_z_water_nm2_sem,msd_z_raw_water_nm2_mean,msd_z_raw_water_nm2_sem,msd_r2_water_nm2_mean,msd_r2_water_nm2_sem,"
           "msd_x_na_nm2_mean,msd_x_na_nm2_sem,msd_y_na_nm2_mean,msd_y_na_nm2_sem,msd_z_na_nm2_mean,msd_z_na_nm2_sem,msd_z_raw_na_nm2_mean,msd_z_raw_na_nm2_sem,msd_r2_na_nm2_mean,msd_r2_na_nm2_sem,"
           "msd_x_cl_nm2_mean,msd_x_cl_nm2_sem,msd_y_cl_nm2_mean,msd_y_cl_nm2_sem,msd_z_cl_nm2_mean,msd_z_cl_nm2_sem,msd_z_raw_cl_nm2_mean,msd_z_raw_cl_nm2_sem,msd_r2_cl_nm2_mean,msd_r2_cl_nm2_sem,"
           "max_jump_x_nm_mean,max_jump_x_nm_sem,max_jump_x_water_nm_mean,max_jump_x_water_nm_sem,max_jump_x_na_nm_mean,max_jump_x_na_nm_sem,max_jump_x_cl_nm_mean,max_jump_x_cl_nm_sem\n";

    for (int lag = 1; lag <= cfg_.keep_frames; ++lag) {
        const size_t i = static_cast<size_t>(lag);
        const int64_t n_pairs = lag_time_ps_.n_used[i];

        ofs << lag << "," << lag_time_ps_.mean(lag) << "," << lag_time_ps_.sem(lag) << ","
            << n_pairs << "," << n_samples_by_lag_[i] << ","
            << msd_x_all_.mean(lag) << "," << msd_x_all_.sem(lag) << ","
            << msd_y_all_.mean(lag) << "," << msd_y_all_.sem(lag) << ","
            << msd_z_all_.mean(lag) << "," << msd_z_all_.sem(lag) << ","
            << msd_z_raw_all_.mean(lag) << "," << msd_z_raw_all_.sem(lag) << ","
            << msd_r2_all_.mean(lag) << "," << msd_r2_all_.sem(lag) << ","
            << msd_x_sp_[0].mean(lag) << "," << msd_x_sp_[0].sem(lag) << ","
            << msd_y_sp_[0].mean(lag) << "," << msd_y_sp_[0].sem(lag) << ","
            << msd_z_sp_[0].mean(lag) << "," << msd_z_sp_[0].sem(lag) << ","
            << msd_z_raw_sp_[0].mean(lag) << "," << msd_z_raw_sp_[0].sem(lag) << ","
            << msd_r2_sp_[0].mean(lag) << "," << msd_r2_sp_[0].sem(lag) << ","
            << msd_x_sp_[1].mean(lag) << "," << msd_x_sp_[1].sem(lag) << ","
            << msd_y_sp_[1].mean(lag) << "," << msd_y_sp_[1].sem(lag) << ","
            << msd_z_sp_[1].mean(lag) << "," << msd_z_sp_[1].sem(lag) << ","
            << msd_z_raw_sp_[1].mean(lag) << "," << msd_z_raw_sp_[1].sem(lag) << ","
            << msd_r2_sp_[1].mean(lag) << "," << msd_r2_sp_[1].sem(lag) << ","
            << msd_x_sp_[2].mean(lag) << "," << msd_x_sp_[2].sem(lag) << ","
            << msd_y_sp_[2].mean(lag) << "," << msd_y_sp_[2].sem(lag) << ","
            << msd_z_sp_[2].mean(lag) << "," << msd_z_sp_[2].sem(lag) << ","
            << msd_z_raw_sp_[2].mean(lag) << "," << msd_z_raw_sp_[2].sem(lag) << ","
            << msd_r2_sp_[2].mean(lag) << "," << msd_r2_sp_[2].sem(lag) << ","
            << max_jump_x_all_.mean(lag) << "," << max_jump_x_all_.sem(lag) << ","
            << max_jump_x_sp_[0].mean(lag) << "," << max_jump_x_sp_[0].sem(lag) << ","
            << max_jump_x_sp_[1].mean(lag) << "," << max_jump_x_sp_[1].sem(lag) << ","
            << max_jump_x_sp_[2].mean(lag) << "," << max_jump_x_sp_[2].sem(lag) << "\n";
    }
}

}  // namespace simio::analysis
