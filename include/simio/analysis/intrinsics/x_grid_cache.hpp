#pragma once
#include <cstdint>
#include <vector>

namespace simio::runtime { class CacheStore; }

namespace simio::analysis::intrinsics {

// X grid used by 1D x-binned properties.
// IMPORTANT: centers_rel are relative to xmin (i.e., (i+0.5)*dx), matching current outputs.
struct XGrid {
  double xmin{0.0};
  double xmax{0.0};
  int nx{0};
  double dx{0.0};
  std::vector<double> centers_rel; // size nx, (i+0.5)*dx
};

// Build or fetch from CacheStore under key "intrinsic.x_grid".
// params_hash depends on xmin, xmax, nx.
XGrid get_or_build_x_grid(simio::runtime::CacheStore& cache, double xmin, double xmax, int nx);

} // namespace simio::analysis::intrinsics
