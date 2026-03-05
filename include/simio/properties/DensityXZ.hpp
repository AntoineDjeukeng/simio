#pragma once

#include <cstdint>
#include <vector>

#include "simio/layered_pipeline.hpp"
#include "simio/properties/GatingTypes.hpp"

namespace simio::properties {

enum class DensityXMode : uint8_t {
    ChannelWindow = 0,
    FullBox = 1,
};

struct DensityRoi {
    double x_min = 0.0;
    double x_max = 0.0;
    double z_min = 0.0;
    double z_max = 0.0;
    DensityXMode x_mode = DensityXMode::ChannelWindow;
};

struct DensityXZConfig {
    MoleculeSelection selection = MoleculeSelection::All;
    DensityRoi roi{};
    int nx = 32;
    int nz = 16;
    bool normalize_number_density = true;
};

struct DensityXZFrame {
    int nx = 0;
    int nz = 0;
    double dx = 0.0;
    double dz = 0.0;
    double bin_volume = 0.0;
    int64_t selected_total = 0;
    int64_t binned_total = 0;
    int64_t bin_oob_count = 0;
    int64_t selected_water = 0;
    int64_t selected_cation = 0;
    int64_t selected_anion = 0;
    int64_t selected_other = 0;
    std::vector<int64_t> counts;
    std::vector<double> rho;
};

class DensityXZProperty final : public layered::PropertyKernel {
  public:
    explicit DensityXZProperty(const DensityXZConfig& cfg = {});

    const char* name() const override;
    uint64_t requires() const override;
    void compute_frame(const layered::FrameContext& ctx) override;
    void finalize() override;

    const std::vector<DensityXZFrame>& frames() const;

  private:
    DensityXZConfig cfg_{};
    std::vector<uint8_t> prev_in_channel_{};
    std::vector<DensityXZFrame> frames_{};
};

}  // namespace simio::properties
