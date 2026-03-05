#include "simio/runtime/cache.hpp"
#include <cassert>
#include <cstdint>
#include <iostream>

int main() {
  using namespace simio::runtime;

  CacheStore c;
  CacheKey k1{"nodeA", CacheScope::Frame, 123, 0};

  // miss
  assert(c.get(k1) == nullptr);
  assert(c.stats().gets == 1);
  assert(c.stats().hits == 0);

  // put + hit
  Blob b;
  b.bytes = {1,2,3,4};
  c.put(k1, b);
  assert(c.stats().puts == 1);

  auto* out = c.get(k1);
  assert(out != nullptr);
  assert(out->bytes.size() == 4);
  assert(out->bytes[0] == 1);

  std::cout << "OK\n";
  return 0;
}
