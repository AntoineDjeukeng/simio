#pragma once

#include <cstdint>
#include <vector>

#include "simio/layered_pipeline.hpp"
#include "simio/properties/GatingTypes.hpp"

namespace simio::properties {

struct DipoleXZConfig {
    MoleculeSelection selection = MoleculeSelection::All;
};

struct DipoleXZFrame {
    int nx = 0;
    int nz = 0;
    int64_t total = 0;
    double mean_mu = 0.0;
    std::vector<int64_t> counts;
    std::vector<double> sum_mu;
};

class DipoleXZProperty final : public layered::PropertyKernel {
  public:
    explicit DipoleXZProperty(const DipoleXZConfig& cfg = {});

    const char* name() const override;
    uint64_t requires() const override;
    void compute_frame(const layered::FrameContext& ctx) override;
    void finalize() override;

    const std::vector<DipoleXZFrame>& frames() const;

  private:
    DipoleXZConfig cfg_{};
    std::vector<uint8_t> prev_in_channel_{};
    std::vector<DipoleXZFrame> frames_{};
};

}  // namespace simio::properties

