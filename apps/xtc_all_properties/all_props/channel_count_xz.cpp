#include "channel_count_xz.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace simio::analysis {

ChannelCountXZAnalyzer::ChannelCountXZAnalyzer(const ChannelCountXZConfig& cfg) : cfg_(cfg) {}

void ChannelCountXZAnalyzer::process_frame(const Topology& topo, const Frame& fr,
                                           const std::vector<MolState>& ms, int frame_idx) {
    const double Lx = fr.pbc.L[0];
    const double Lz = fr.pbc.L[2];
    if (Lx <= 0.0 || Lz <= 0.0) {
        throw std::runtime_error("ChannelCountXZAnalyzer: invalid box lengths");
    }

    const double xminw = fr.pbc.wrap_pos(0, cfg_.xmin);
    const double xmaxw = fr.pbc.wrap_pos(0, cfg_.xmax);
    const double zminw = fr.pbc.wrap_pos(2, cfg_.zmin);
    const double zmaxw = fr.pbc.wrap_pos(2, cfg_.zmax);

    FrameRow row{};
    row.frame_idx = frame_idx;
    row.step = fr.step;
    row.time_ps = fr.time_ps;

    for (size_t mid = 0; mid < topo.mols.size(); ++mid) {
        const int sid = species_index_from_type(topo.mols[mid].type);
        if (sid < 0 || sid >= 3) continue;

        const MolCache& c = ms[mid].cache;
        if (!(c.flags & MolCache::HAS_KEY)) continue;

        Vec3d key = c.key_wrapped;
        fr.pbc.wrap_pos3(key);
        if (!in_range_pbc(key.v[0], xminw, xmaxw)) continue;
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
