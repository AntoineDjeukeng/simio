#include "simio/properties/DipoleXZ.hpp"

#include <utility>

namespace simio::properties {

DipoleXZProperty::DipoleXZProperty(const DipoleXZConfig& cfg) : cfg_(cfg) {}

const char* DipoleXZProperty::name() const { return "DipoleXZProperty"; }

uint64_t DipoleXZProperty::requires() const {
    return layered::field_mask(layered::Field::Layer0Bins) |
           layered::field_mask(layered::Field::Dipole) |
           layered::field_mask(layered::Field::InChannel);
}

void DipoleXZProperty::compute_frame(const layered::FrameContext& ctx) {
    const size_t nmol = ctx.topo.mols.size();
    if (prev_in_channel_.size() != nmol) {
        prev_in_channel_.assign(nmol, 0);
    }

    DipoleXZFrame out;
    out.nx = ctx.layer0.nx;
    out.nz = ctx.layer0.nz;
    out.counts.assign((size_t)(out.nx * out.nz), 0);
    out.sum_mu.assign((size_t)(out.nx * out.nz), 0.0);

    double total_mu = 0.0;
    for (size_t mid = 0; mid < nmol; ++mid) {
        const bool in_now = (ctx.layer1.is_in_channel[mid] != 0);
        const bool in_prev = (prev_in_channel_[mid] != 0);
        const bool include = include_molecule(cfg_.selection, in_now, in_prev);

        if (include && ctx.layer1.has_dipole[mid]) {
            const int bin = ctx.layer0.bin_id[mid];
            if (bin >= 0 && bin < out.nx * out.nz) {
                const double mu = norm3(ctx.layer1.dipole[mid]);
                ++out.counts[(size_t)bin];
                out.sum_mu[(size_t)bin] += mu;
                ++out.total;
                total_mu += mu;
            }
        }

        prev_in_channel_[mid] = in_now ? 1 : 0;
    }

    if (out.total > 0) out.mean_mu = total_mu / (double)out.total;
    frames_.push_back(std::move(out));
}

void DipoleXZProperty::finalize() {}

const std::vector<DipoleXZFrame>& DipoleXZProperty::frames() const { return frames_; }

}  // namespace simio::properties
