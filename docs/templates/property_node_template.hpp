#pragma once
#include "simio/analysis/node.hpp"
#include "simio/analysis/intrinsics/context.hpp"

namespace simio::analysis {

class MyPropertyNode : public INode {
public:

  NodeDesc desc() const override
  {
    return NodeDesc{
      "property.my_property",
      {"intrinsic.x_grid"}
    };
  }

  void compute(simio::runtime::CacheStore& cache) override
  {
    simio::analysis::intrinsics::IntrinsicContext ctx{cache};

    // fetch intrinsics
    // compute property
    // store outputs
  }
};

}
