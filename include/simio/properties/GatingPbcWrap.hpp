#pragma once

#include <vector>

#include "simio/layered_pipeline.hpp"
#include "simio/properties/GatingTypes.hpp"

namespace simio::properties {

struct GatingPbcWrapConfig {
    MoleculeSelection selection = MoleculeSelection::InChannelNow;
};

class GatingPbcWrapProperty final : public layered::PropertyKernel {
  public:
    explicit GatingPbcWrapProperty(const GatingPbcWrapConfig& cfg = {});

    const char* name() const override;
    uint64_t requires() const override;
    void compute_frame(const layered::FrameContext& ctx) override;
    void finalize() override;

    const std::vector<GatingFrameCount>& frames() const;
    const GatingFrameCount& cumulative() const;

  private:
    GatingPbcWrapConfig cfg_{};
    std::vector<Vec3d> prev_key_time_wrapped_{};
    std::vector<uint8_t> has_prev_{};
    std::vector<uint8_t> prev_in_channel_{};
    std::vector<GatingFrameCount> frames_{};
    GatingFrameCount cumulative_{};
};

}  // namespace simio::properties
