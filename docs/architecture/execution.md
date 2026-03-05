# Execution architecture

This repository follows a **compute-once** model: reusable computations are implemented as **intrinsics** and cached; analyzers/properties must consume intrinsics rather than recomputing equivalent logic.

## Entry point and configuration

All runs are configured through `simio::runtime::RunConfig` (CLI or JSON), currently used by:

- `apps/xtc_all_properties/main.cpp`

Rule: parsing/validation belongs in `RunConfig`; analyzers receive typed config structs derived from it.

## Cache model

Cache is provided by `simio::runtime::CacheStore`.

A cached value is addressed by:

- `node_id` (stable string, e.g. `intrinsic.x_grid`)
- `scope` (`Global`, `Frame`, `Window`)
- `params_hash` (hash of all parameters that affect the result)
- `instance_id` (frame index/window id; `-1` for non-instance scoped)

Strict mode (`CacheStore::put_strict`) is used for intrinsics to detect duplicate computation.

Stats are available via `[simio-cache]` log lines when `SIMIO_CACHE_STATS=1`.

## Intrinsics

Intrinsics live under:

- `include/simio/analysis/intrinsics/`
- `src/analysis/intrinsics/`

They are accessed via `IntrinsicContext`, which wraps a shared `CacheStore`.

Current intrinsics:

### Global
- `x_grid` : cached x bin centers/dx (relative centers)
- `z_grid` : cached z bin centers/dz (relative centers)
- `channel_roi_x` : cached PBC interval helpers (`contains_x`, `map_x_to_channel`)

### Frame
- `in_channel_mask_x` : cached per-frame membership mask for x-channel inclusion
  - intended to be reused by multiple analyzers in the same run/frame

## Analyzer/property rules

### Compute-once invariant
If an analyzer needs a reusable quantity (bin grids, ROI mapping, membership tests), it must be implemented as an intrinsic and retrieved through the cache.

Analyzers must **not** implement ad-hoc versions of the same computation.

### Hot-loop performance
- No per-frame dynamic allocations inside hot loops.
- Use preallocated member buffers (e.g., `xw_tmp_`).
- Intrinsic retrieval should happen once per frame at most (and cached).

## How to add a new intrinsic

1. Choose scope: `Global`, `Frame`, or `Window`.
2. Define `node_id` and `params_hash` (include all parameters).
3. Use `CacheStore::get` then build if missing.
4. Store with `put_strict` to enforce compute-once.
5. Add a unit test in `tests/unit/` for:
   - cache hit behavior
   - strict duplicate detection
6. Wire the new source into the build (Makefile/CMake if needed).

Templates: `docs/templates/intrinsic_template.*`

## How to add a new analyzer/property

1. Add a config struct and validate inputs via `RunConfig`.
2. Identify required intrinsics and consume them via `IntrinsicContext`.
3. Ensure no duplicated computation (prefer intrinsics).
4. Add to `apps/xtc_all_properties` wiring.
5. Validate with:
   - `scripts/run_smoke.sh`
   - `scripts/compare_outputs.py` against baseline
6. Commit in a single focused commit with a clear message.

Templates: `docs/templates/property_node_template.hpp`
