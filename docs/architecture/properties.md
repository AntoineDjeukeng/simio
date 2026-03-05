# Properties

A property is an analysis output (CSV, profiles, fluxes, MSD/VACF, etc.) computed from frames, topology, and intrinsics.

## Rules

1. Properties must not recompute intrinsics locally.
2. Hot loops must avoid cache lookups. Fetch intrinsics once (constructor/init), store results locally.
3. If multiple properties share an intermediate computation, it must become an intrinsic.

## Adding a new property (current analyzer style)

1. Implement analyzer in `apps/xtc_all_properties/all_props/`
2. Fetch intrinsics in constructor using `IntrinsicContext`
3. Store intrinsic results (dx, centers, etc.)
4. Use stored values in `process_frame()`
5. Write output in `write_csv()`

## Future direction (node graph)

Long-term properties will derive from `INode`:

- declare `NodeDesc` (id + dependencies)
- implement `compute(CacheStore&)`
- scheduler executes dependencies automatically
