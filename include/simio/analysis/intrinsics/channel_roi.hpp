#pragma once
#include <cstdint>

namespace simio::runtime { class CacheStore; }

namespace simio::analysis::intrinsics {

struct IntrinsicContext;

struct ChannelRoiX {
  double xmin{0.0};
  double xmax{0.0};
  double Lx{0.0};

  // normalized bounds in [0,Lx)
  double xmin_w{0.0};
  double xmax_w{0.0};
  bool wraps{false};
  double xlen{0.0}; // length along PBC interval

  // x is assumed in [0,Lx) (caller should wrap if needed)
  bool contains_x(double x) const;
  // maps x into s in [0, xlen) along the interval
  double map_x_to_channel(double x) const;
};

ChannelRoiX get_channel_roi_x(IntrinsicContext& ctx, double xmin, double xmax, double Lx);
ChannelRoiX get_or_build_channel_roi_x(simio::runtime::CacheStore& cache, double xmin, double xmax, double Lx);

} // namespace simio::analysis::intrinsics
