#include "simio/analysis/intrinsics/in_channel_mask.hpp"
#include "simio/analysis/intrinsics/context.hpp"
#include "simio/analysis/intrinsics/channel_roi.hpp"
#include "simio/runtime/cache.hpp"

#include <cstdint>
#include <cstring>

namespace simio::analysis::intrinsics {

static std::uint64_t u64_of_double(double x) {
  std::uint64_t u;
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  std::memcpy(&u, &x, sizeof(u));
  return u;
}

static std::uint64_t hash_mask(double xmin, double xmax, double Lx, std::size_t n) {
  std::uint64_t h = 1469598103934665603ULL;
  auto mix = [&](std::uint64_t v) { h ^= v; h *= 1099511628211ULL; };
  mix(u64_of_double(xmin));
  mix(u64_of_double(xmax));
  mix(u64_of_double(Lx));
  mix((std::uint64_t)n);
  return h;
}

XMask get_in_channel_mask_x(
    IntrinsicContext& ctx,
    std::int64_t frame_index,
    const double* x_wrapped, std::size_t n,
    double xmin, double xmax, double Lx)
{
  using namespace simio::runtime;
  auto& cache = ctx.cache;

  CacheKey k{
    "intrinsic.in_channel_mask_x",
    CacheScope::Frame,
    hash_mask(xmin, xmax, Lx, n),
    frame_index
  };

  if (const Blob* b = cache.get(k)) {
    if (b->bytes.size() == n) {
      XMask out;
      out.frame_index = frame_index;
      out.in_channel.assign(b->bytes.begin(), b->bytes.end());
      return out;
    }
  }

  // Build using canonical ROI
  ChannelRoiX roi = get_channel_roi_x(ctx, xmin, xmax, Lx);

  XMask out;
  out.frame_index = frame_index;
  out.in_channel.resize(n);

  for (std::size_t i = 0; i < n; ++i) {
    const double xw = x_wrapped[i]; // assumed wrapped to [0,Lx)
    out.in_channel[i] = (std::uint8_t)(roi.contains_x(xw) ? 1 : 0);
  }

  Blob blob;
  blob.bytes.assign(out.in_channel.begin(), out.in_channel.end());
  (void)cache.put_strict(k, std::move(blob));
  return out;
}

XMask get_or_build_in_channel_mask_x(
    simio::runtime::CacheStore& cache,
    std::int64_t frame_index,
    const double* x_wrapped, std::size_t n,
    double xmin, double xmax, double Lx)
{
  IntrinsicContext ctx{cache};
  return get_in_channel_mask_x(ctx, frame_index, x_wrapped, n, xmin, xmax, Lx);
}

} // namespace simio::analysis::intrinsics
