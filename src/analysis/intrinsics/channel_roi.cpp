#include "simio/analysis/intrinsics/channel_roi.hpp"
#include "simio/analysis/intrinsics/context.hpp"
#include "simio/runtime/cache.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>

namespace simio::analysis::intrinsics {

static double wrap_0_L(double x, double L) {
  if (L <= 0.0) return x;
  x = std::fmod(x, L);
  if (x < 0.0) x += L;
  return x;
}

bool ChannelRoiX::contains_x(double x) const {
  // Match existing PBC interval semantics: if wraps, interval is [xmin_w,Lx) U [0,xmax_w]
  if (!wraps) return (x >= xmin_w && x <= xmax_w);
  return (x >= xmin_w || x <= xmax_w);
}

double ChannelRoiX::map_x_to_channel(double x) const {
  // Map x along the interval starting at xmin_w.
  // For non-wrapping: s = x - xmin_w
  // For wrapping: if x >= xmin_w => s = x - xmin_w else s = (Lx - xmin_w) + x
  if (!wraps) return x - xmin_w;
  if (x >= xmin_w) return x - xmin_w;
  return (Lx - xmin_w) + x;
}

static std::uint64_t u64_of_double(double x) {
  std::uint64_t u;
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  std::memcpy(&u, &x, sizeof(u));
  return u;
}

static std::uint64_t hash_roi(double xmin, double xmax, double Lx) {
  std::uint64_t h = 1469598103934665603ULL;
  auto mix = [&](std::uint64_t v) { h ^= v; h *= 1099511628211ULL; };
  mix(u64_of_double(xmin));
  mix(u64_of_double(xmax));
  mix(u64_of_double(Lx));
  return h;
}

ChannelRoiX get_channel_roi_x(IntrinsicContext& ctx, double xmin, double xmax, double Lx) {
  using namespace simio::runtime;
  auto& cache = ctx.cache;

  CacheKey k{"intrinsic.channel_roi_x", CacheScope::Global, hash_roi(xmin, xmax, Lx), -1};

  if (const Blob* b = cache.get(k)) {
    if (b->bytes.size() == sizeof(ChannelRoiX)) {
      ChannelRoiX roi;
      std::memcpy(&roi, b->bytes.data(), sizeof(ChannelRoiX));
      return roi;
    }
  }

  ChannelRoiX roi;
  roi.xmin = xmin; roi.xmax = xmax; roi.Lx = Lx;
  roi.xmin_w = wrap_0_L(xmin, Lx);
  roi.xmax_w = wrap_0_L(xmax, Lx);
  roi.wraps = (roi.xmin_w > roi.xmax_w);

  if (!roi.wraps) roi.xlen = roi.xmax_w - roi.xmin_w;
  else roi.xlen = (Lx - roi.xmin_w) + roi.xmax_w;

  Blob out;
  out.bytes.resize(sizeof(ChannelRoiX));
  std::memcpy(out.bytes.data(), &roi, sizeof(ChannelRoiX));
  (void)cache.put_strict(k, std::move(out));
  return roi;
}

ChannelRoiX get_or_build_channel_roi_x(simio::runtime::CacheStore& cache, double xmin, double xmax, double Lx) {
  IntrinsicContext ctx{cache};
  return get_channel_roi_x(ctx, xmin, xmax, Lx);
}

} // namespace simio::analysis::intrinsics
