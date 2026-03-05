# Intrinsics

Intrinsics are reusable intermediate computations shared by multiple properties (e.g., x-binning grid, molecule COM, region flags).
They are *not* output properties; they exist to eliminate duplicated work and make dependencies explicit.

## Naming
- Cache key node id MUST be prefixed with: `intrinsic.`
  - Example: `intrinsic.x_grid`

## Scopes
Use `simio::runtime::CacheScope`:
- `Global`: depends only on run configuration/topology; computed once per run.
- `Frame`: depends on current frame; computed once per frame.
- `Window`: depends on a time window; computed once per window (ring-buffer sized).

## Cache key rules
A cache key is:
- `node_id`: stable intrinsic name (e.g., `intrinsic.x_grid`)
- `scope`: Global/Frame/Window
- `params_hash`: hashes the intrinsic parameters (binning bounds, cutoffs, selection, etc.)
- `instance_id`: frame index / window id; `-1` for Global

### Compute-once invariant
Intrinsics MUST be stored using `CacheStore::put_strict()` so duplicates are detected.

## API pattern
Preferred access pattern is context-based:
```cpp
simio::analysis::intrinsics::IntrinsicContext ictx{cache};
auto grid = simio::analysis::intrinsics::get_x_grid(ictx, xmin, xmax, nx);
```
