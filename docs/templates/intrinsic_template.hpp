#pragma once

namespace simio::analysis::intrinsics {

struct IntrinsicContext;

struct MyIntrinsic {
  int value{0};
};

MyIntrinsic get_my_intrinsic(IntrinsicContext& ctx, int param);

}
