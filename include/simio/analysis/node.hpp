#pragma once
#include <string>
#include <vector>

namespace simio::runtime { class CacheStore; }

namespace simio::analysis {

struct NodeDesc {
  std::string id;                 // stable node id (e.g., "density_x")
  std::vector<std::string> deps;  // dependency node ids
};

class INode {
public:
  virtual ~INode() = default;
  virtual NodeDesc desc() const = 0;

  // Note: We will introduce typed contexts later.
  virtual void compute(simio::runtime::CacheStore& cache) = 0;
};

} // namespace simio::analysis
