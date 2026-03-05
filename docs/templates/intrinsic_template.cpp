#include "simio/analysis/intrinsics/context.hpp"
#include "simio/runtime/cache.hpp"
#include <cstring>

#include "intrinsic_template.hpp"

namespace simio::analysis::intrinsics {

MyIntrinsic get_my_intrinsic(IntrinsicContext& ctx, int param)
{
  using namespace simio::runtime;

  auto& cache = ctx.cache;

  CacheKey key{
    "intrinsic.my_intrinsic",
    CacheScope::Global,
    (std::uint64_t)param,
    -1
  };

  if (const Blob* b = cache.get(key)) {
    MyIntrinsic out;
    std::memcpy(&out.value, b->bytes.data(), sizeof(int));
    return out;
  }

  MyIntrinsic out;
  out.value = param;

  Blob blob;
  blob.bytes.resize(sizeof(int));
  std::memcpy(blob.bytes.data(), &out.value, sizeof(int));

  cache.put_strict(key, std::move(blob));
  return out;
}

}
