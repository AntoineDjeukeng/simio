#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

namespace simio::runtime { class CacheStore; }

namespace simio::analysis::intrinsics {

struct IntrinsicContext;

struct XMask {
  std::int64_t frame_index{-1};
  std::vector<std::uint8_t> in_channel; // size n
};

XMask get_in_channel_mask_x(
    IntrinsicContext& ctx,
    std::int64_t frame_index,
    const double* x_wrapped, std::size_t n,
    double xmin, double xmax, double Lx);

XMask get_or_build_in_channel_mask_x(
    simio::runtime::CacheStore& cache,
    std::int64_t frame_index,
    const double* x_wrapped, std::size_t n,
    double xmin, double xmax, double Lx);

} // namespace simio::analysis::intrinsics
