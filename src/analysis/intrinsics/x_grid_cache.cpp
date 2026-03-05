#include "simio/analysis/intrinsics/x_grid_cache.hpp"
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

static std::uint64_t hash_xgrid(double xmin, double xmax, int nx) {
  // FNV-1a style mix
  std::uint64_t h = 1469598103934665603ULL;
  auto mix = [&](std::uint64_t v) {
    h ^= v;
    h *= 1099511628211ULL;
  };
  mix(u64_of_double(xmin));
  mix(u64_of_double(xmax));
  mix((std::uint64_t)(std::uint32_t)nx);
  return h;
}

XGrid get_or_build_x_grid(simio::runtime::CacheStore& cache, double xmin, double xmax, int nx) {
  using namespace simio::runtime;

  CacheKey k{"intrinsic.x_grid", CacheScope::Global, hash_xgrid(xmin, xmax, nx), -1};

  if (const Blob* b = cache.get(k)) {
    const std::size_t need = sizeof(double)*4 + sizeof(std::int32_t) + sizeof(double)*(std::size_t)nx;
    if (b->bytes.size() == need) {
      XGrid g;
      const std::uint8_t* p = b->bytes.data();
      std::memcpy(&g.xmin, p, sizeof(double)); p += sizeof(double);
      std::memcpy(&g.xmax, p, sizeof(double)); p += sizeof(double);
      std::memcpy(&g.dx,   p, sizeof(double)); p += sizeof(double);
      double nx_as_double = 0.0;
      std::memcpy(&nx_as_double, p, sizeof(double)); p += sizeof(double); // reserved for future; keep layout stable
      std::int32_t nx_i = 0;
      std::memcpy(&nx_i, p, sizeof(nx_i)); p += sizeof(nx_i);
      if (nx_i == nx) {
        g.nx = nx;
        g.centers_rel.resize((std::size_t)nx);
        std::memcpy(g.centers_rel.data(), p, sizeof(double)*(std::size_t)nx);
        return g;
      }
    }
    // If blob mismatch, fallthrough to rebuild.
  }

  // Build matching current behavior:
  // Lx = xmax - xmin
  // dx = Lx / nx
  // centers_rel[i] = (i+0.5)*dx  (RELATIVE to xmin)
  XGrid g;
  g.xmin = xmin; g.xmax = xmax; g.nx = nx;
  const double Lx = (xmax - xmin);
  g.dx = Lx / (double)nx;
  g.centers_rel.resize((std::size_t)nx);
  for (int i = 0; i < nx; ++i) {
    g.centers_rel[(std::size_t)i] = ((double)i + 0.5) * g.dx;
  }

  Blob out;
  out.bytes.resize(sizeof(double)*4 + sizeof(std::int32_t) + sizeof(double)*(std::size_t)nx);
  std::uint8_t* q = out.bytes.data();
  std::memcpy(q, &g.xmin, sizeof(double)); q += sizeof(double);
  std::memcpy(q, &g.xmax, sizeof(double)); q += sizeof(double);
  std::memcpy(q, &g.dx,   sizeof(double)); q += sizeof(double);
  double reserved = 0.0;
  std::memcpy(q, &reserved, sizeof(double)); q += sizeof(double);
  std::int32_t nx_i = (std::int32_t)nx;
  std::memcpy(q, &nx_i, sizeof(nx_i)); q += sizeof(nx_i);
  std::memcpy(q, g.centers_rel.data(), sizeof(double)*(std::size_t)nx);

  // Use strict mode if enabled by env var in main (we don't read env here).
  // The caller controls strictness by choosing put vs put_strict; here we default to put_strict always.
  // Duplicate here indicates a programming error (same run key built twice).
  (void)cache.put_strict(k, std::move(out));

  return g;
}

} // namespace simio::analysis::intrinsics
