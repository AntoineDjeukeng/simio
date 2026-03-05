#include "simio/analysis/intrinsics/z_grid_cache.hpp"
#include "simio/analysis/intrinsics/context.hpp"
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

static std::uint64_t hash_zgrid(double zmin, double zmax, int nz) {
  std::uint64_t h = 1469598103934665603ULL;
  auto mix = [&](std::uint64_t v) {
    h ^= v;
    h *= 1099511628211ULL;
  };
  mix(u64_of_double(zmin));
  mix(u64_of_double(zmax));
  mix((std::uint64_t)(std::uint32_t)nz);
  return h;
}

ZGrid get_z_grid(IntrinsicContext& ctx, double zmin, double zmax, int nz) {
  using namespace simio::runtime;
  auto& cache = ctx.cache;

  CacheKey k{"intrinsic.z_grid", CacheScope::Global, hash_zgrid(zmin, zmax, nz), -1};

  if (const Blob* b = cache.get(k)) {
    const std::size_t need = sizeof(double)*4 + sizeof(std::int32_t) + sizeof(double)*(std::size_t)nz;
    if (b->bytes.size() == need) {
      ZGrid g;
      const std::uint8_t* p = b->bytes.data();
      std::memcpy(&g.zmin, p, sizeof(double)); p += sizeof(double);
      std::memcpy(&g.zmax, p, sizeof(double)); p += sizeof(double);
      std::memcpy(&g.dz,   p, sizeof(double)); p += sizeof(double);
      double reserved = 0.0;
      std::memcpy(&reserved, p, sizeof(double)); p += sizeof(double);
      std::int32_t nz_i = 0;
      std::memcpy(&nz_i, p, sizeof(nz_i)); p += sizeof(nz_i);
      if (nz_i == nz) {
        g.nz = nz;
        g.centers_rel.resize((std::size_t)nz);
        std::memcpy(g.centers_rel.data(), p, sizeof(double)*(std::size_t)nz);
        return g;
      }
    }
  }

  ZGrid g;
  g.zmin = zmin; g.zmax = zmax; g.nz = nz;
  const double Lz = (zmax - zmin);
  g.dz = Lz / (double)nz;
  g.centers_rel.resize((std::size_t)nz);
  for (int i = 0; i < nz; ++i) {
    g.centers_rel[(std::size_t)i] = ((double)i + 0.5) * g.dz;
  }

  Blob out;
  out.bytes.resize(sizeof(double)*4 + sizeof(std::int32_t) + sizeof(double)*(std::size_t)nz);
  std::uint8_t* q = out.bytes.data();
  std::memcpy(q, &g.zmin, sizeof(double)); q += sizeof(double);
  std::memcpy(q, &g.zmax, sizeof(double)); q += sizeof(double);
  std::memcpy(q, &g.dz,   sizeof(double)); q += sizeof(double);
  double reserved = 0.0;
  std::memcpy(q, &reserved, sizeof(double)); q += sizeof(double);
  std::int32_t nz_i = (std::int32_t)nz;
  std::memcpy(q, &nz_i, sizeof(nz_i)); q += sizeof(nz_i);
  std::memcpy(q, g.centers_rel.data(), sizeof(double)*(std::size_t)nz);

  (void)cache.put_strict(k, std::move(out));
  return g;
}

ZGrid get_or_build_z_grid(simio::runtime::CacheStore& cache, double zmin, double zmax, int nz) {
  IntrinsicContext ctx{cache};
  return get_z_grid(ctx, zmin, zmax, nz);
}

} // namespace simio::analysis::intrinsics
