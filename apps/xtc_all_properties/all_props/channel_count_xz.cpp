#include "channel_count_xz.hpp"

#include "simio/analysis/intrinsics/context.hpp"
#include "simio/analysis/intrinsics/in_channel_mask.hpp"
#include "simio/runtime/cache.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace simio::analysis {

ChannelCountXZAnalyzer::ChannelCountXZAnalyzer(const ChannelCountXZConfig& cfg) : cfg_(cfg) {}

ChannelCountXZAnalyzer::ChannelCountXZAnalyzer(const ChannelCountXZConfig& cfg,
                                               simio::runtime::CacheStore& cache)
    : ChannelCountXZAnalyzer(cfg) {
    cache_ = &cache;
}

void ChannelCountXZAnalyzer::process_frame(const Topology& topo, const Frame& fr,
                                           const std::vector<MolState>& ms, int frame_idx) {
    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Lz <= 0.0) {
        throw std::runtime_error("ChannelCountXZAnalyzer: invalid box lengths");
    }

    if (!cache_) {
        static thread_local simio::runtime::CacheStore fallback_cache;
        cache_ = &fallback_cache;
    }

    const double zminw = fr.pbc.wrap_pos(2, cfg_.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg_.zmax);

    FrameRow row{};
    row.frame_idx = frame_idx;
    row.step = fr.step;
    row.time_ps = fr.time_ps;

    const size_t nmol = topo.mols.size();
    if (xw_tmp_.size() != nmol) xw_tmp_.assign(nmol, 0.0);
    for (size_t mid = 0; mid < nmol; ++mid) {
        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) {
            xw_tmp_[mid] = 0.0;
            continue;
        }
        Vec3d key = c.key_wrapped;
        fr.pbc.wrap_pos3(key);
        xw_tmp_[mid] = key.v[0];
    }

    const std::int64_t frame_id = frame_counter_++;
    simio::analysis::intrinsics::IntrinsicContext ictx{*cache_};
    const auto mask = simio::analysis::intrinsics::get_in_channel_mask_x(
        ictx,
        frame_id,
        xw_tmp_.data(),
        xw_tmp_.size(),
        cfg_.xmin,
        cfg_.xmax,
        Lx);
    if (mask.in_channel.size() != nmol) {
        throw std::runtime_error("ChannelCountXZAnalyzer: in_channel_mask size mismatch");
    }

    for (size_t mid = 0; mid < nmol; ++mid) {
        const int sid = species_index_from_type(topo.mols[mid].type);
        if (sid < 0 || sid >= 3) continue;

        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) continue;

        Vec3d key = c.key_wrapped;
        fr.pbc.wrap_pos3(key);
        if (mask.in_channel[mid] == 0u) continue;
        if (!in_range_pbc(key.v[2], zminw, zmaxw)) continue;

        row.count[static_cast<size_t>(sid)] += 1;
    }

    rows_.push_back(row);
    nframes_ += 1;
}

void ChannelCountXZAnalyzer::write_csv(const std::string& path) const {
    if (nframes_ <= 0) throw std::runtime_error("ChannelCountXZAnalyzer: no frames processed");
    std::ofstream ofs(path);
    if (!ofs) throw std::runtime_error("ChannelCountXZAnalyzer: failed to open output CSV: " + path);

    ofs << std::setprecision(12);
    ofs << "frame_idx,step,time_ps,sol,na,cl\n";
    for (const FrameRow& r : rows_) {
        ofs << r.frame_idx << "," << r.step << "," << r.time_ps << "," << r.count[0] << ","
            << r.count[1] << "," << r.count[2] << "\n";
    }
}

}  // namespace simio::analysis
