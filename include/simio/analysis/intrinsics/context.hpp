#pragma once

namespace simio::runtime { class CacheStore; }

namespace simio::analysis::intrinsics {

struct IntrinsicContext {
  simio::runtime::CacheStore& cache;
};

} // namespace simio::analysis::intrinsics
