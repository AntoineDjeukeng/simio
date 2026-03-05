#include "simio/runtime/cache.hpp"
#include <cassert>
#include <iostream>

int main() {
  using namespace simio::runtime;

  CacheStore c;
  CacheKey k{"nodeA", CacheScope::Frame, 1, 0};

  Blob b1; b1.bytes = {1};
  Blob b2; b2.bytes = {2};

  bool ok1 = c.put_strict(k, std::move(b1));
  bool ok2 = c.put_strict(k, std::move(b2));

  assert(ok1 == true);
  assert(ok2 == false);
  assert(c.stats().duplicate_puts == 1);

  auto* out = c.get(k);
  assert(out && out->bytes.size() == 1 && out->bytes[0] == 1);

  std::cout << "OK\n";
  return 0;
}
