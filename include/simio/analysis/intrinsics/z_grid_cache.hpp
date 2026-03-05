#pragma once
#include <cstdint>
#include <vector>

namespace simio::runtime { class CacheStore; }

namespace simio::analysis::intrinsics {

struct IntrinsicContext;

// IMPORTANT: centers_rel are relative to zmin (i.e., (i+0.5)*dz) to match existing CSV semantics.
struct ZGrid {
  double zmin{0.0};
  double zmax{0.0};
  int nz{0};
  double dz{0.0};
  std::vector<double> centers_rel; // size nz, (i+0.5)*dz
};

ZGrid get_or_build_z_grid(simio::runtime::CacheStore& cache, double zmin, double zmax, int nz);
ZGrid get_z_grid(IntrinsicContext& ctx, double zmin, double zmax, int nz);

} // namespace simio::analysis::intrinsics
