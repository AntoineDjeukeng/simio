#pragma once
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace simio::runtime {

// Scope of cached results
enum class CacheScope : uint8_t {
  Frame,   // per-frame cached values
  Window,  // per-window cached values (e.g., VACF accumulators)
  Global   // topology constants, etc.
};

// Minimal cache key: node_id + scope + params_hash + frame/window id
struct CacheKey {
  std::string node_id;
  CacheScope scope{CacheScope::Frame};
  std::uint64_t params_hash{0};
  std::int64_t instance_id{0}; // frame index or window id; -1 for non-instance scoped

  bool operator==(const CacheKey& o) const noexcept {
    return node_id == o.node_id &&
           scope == o.scope &&
           params_hash == o.params_hash &&
           instance_id == o.instance_id;
  }
};

struct CacheKeyHash {
  std::size_t operator()(const CacheKey& k) const noexcept {
    // Simple but stable combination
    std::size_t h = std::hash<std::string>{}(k.node_id);
    h ^= (std::size_t)k.scope + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    h ^= (std::size_t)k.params_hash + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    h ^= (std::size_t)k.instance_id + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
    return h;
  }
};

// Untyped payload for now. We'll replace with typed handles later.
struct Blob {
  std::vector<std::uint8_t> bytes;
};

class CacheStore {
public:
  struct Stats {
    std::uint64_t gets{0};
    std::uint64_t hits{0};
    std::uint64_t puts{0};
    std::uint64_t duplicate_puts{0};
  };

  const Blob* get(const CacheKey& key) {
    ++stats_.gets;
    auto it = map_.find(key);
    if (it == map_.end()) return nullptr;
    ++stats_.hits;
    return &it->second;
  }

  void put(const CacheKey& key, Blob value) {
    ++stats_.puts;
    map_[key] = std::move(value);
  }

  // Strict put: returns false if key already exists (detect duplicate computation)
  bool put_strict(const CacheKey& key, Blob value) {
    ++stats_.puts;
    auto [it, inserted] = map_.emplace(key, std::move(value));
    if (!inserted) {
      ++stats_.duplicate_puts;
      return false;
    }
    return true;
  }

  void clear() { map_.clear(); stats_ = {}; }

  const Stats& stats() const { return stats_; }

private:
  std::unordered_map<CacheKey, Blob, CacheKeyHash> map_;
  Stats stats_{};
};

} // namespace simio::runtime
