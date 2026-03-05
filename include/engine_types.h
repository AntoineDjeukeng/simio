#ifndef ENGINE_TYPES_H
#define ENGINE_TYPES_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum PropertyScope {
    PROPERTY_SCOPE_SYSTEM = 0,
    PROPERTY_SCOPE_FRAME = 1,
    PROPERTY_SCOPE_WINDOW = 2
} PropertyScope;

typedef enum OutputKind {
    OUTPUT_KIND_SCALAR_F64 = 0,
    OUTPUT_KIND_ARRAY_F64 = 1,
    OUTPUT_KIND_ARRAY_F32 = 2,
    OUTPUT_KIND_ARRAY_U32 = 3,
    OUTPUT_KIND_OPAQUE = 4
} OutputKind;

typedef struct MoleculeIndex {
    uint32_t ntypes;
    uint32_t nmol;

    /* Immutable arrays from SYSTEM pool */
    const uint32_t *type_counts;      /* [ntypes] */
    const uint32_t **type_mol_ids;    /* [ntypes][type_counts[t]] */
} MoleculeIndex;

typedef struct SystemTopology {
    uint32_t natoms;
    uint32_t nmol;
    uint32_t ntypes;

    /* Immutable molecule->atom map: atom range of mol m is [off[m], off[m+1]) */
    const uint32_t *mol_atom_offsets; /* [nmol+1] */
    const uint16_t *mol_type;         /* [nmol] */

    const MoleculeIndex *mol_index;

    const double *constants;
    uint32_t nconstants;
} SystemTopology;

typedef struct CellList {
    uint32_t nx;
    uint32_t ny;
    uint32_t nz;
    uint32_t ncells;

    float cell_size[3];
    float inv_cell_size[3];
    uint8_t pbc[3];

    uint32_t max_items;
    uint32_t nitems;

    uint32_t *cell_offsets;   /* [ncells+1] */
    uint32_t *cell_ids;       /* [max_items] */
    uint32_t *mol_cell;       /* [nmol] optional but allocated in this skeleton */

    uint8_t ready;
    uint64_t built_frame_index;
} CellList;

typedef struct FrameSlot {
    float *positions_xyz;     /* [natoms*3], owned by FRAME pool */
    float box[3];
    double time_ps;
    uint64_t step;
    uint64_t frame_index;
    uint8_t occupied;

    float *mol_com_xyz;       /* [nmol*3], primitive cache */
    uint8_t mol_com_ready;
    uint64_t mol_com_frame_index;

    CellList cell_list;       /* primitive cache */
} FrameSlot;

typedef struct FrameRing {
    uint32_t k_slots;
    uint64_t frames_seen;
    FrameSlot *slots;         /* [k_slots] */
} FrameRing;

typedef struct Workspace {
    uint32_t nthreads;
    uint32_t ncells;
    uint32_t nbins;

    uint32_t *counts_tls;     /* [nthreads*ncells] */
    uint32_t *fill_tls;       /* [nthreads*ncells] */
    uint32_t *cell_totals;    /* [ncells] */

    double *bins_tls;         /* [nthreads*nbins] */
    double *bins_reduced;     /* [nbins] */
} Workspace;

typedef struct OutputSpec {
    const char *label;
    PropertyScope scope;
    OutputKind kind;
    uint64_t required_bytes;
} OutputSpec;

typedef struct CacheKey {
    uint32_t node_id;
    PropertyScope scope;
    uint64_t params_hash;

    /* Scope index fields */
    uint64_t frame_index;
    uint64_t window_start;
    uint64_t window_end;
} CacheKey;

typedef struct PropertyOutput {
    void *data;
    uint64_t bytes;
    uint64_t capacity_bytes;
} PropertyOutput;

struct EngineContext;
typedef int (*PropertyComputeFn)(struct EngineContext *ctx, PropertyOutput *out);

typedef struct PropertyNode {
    uint32_t node_id;
    const char *name;
    PropertyScope scope;
    uint64_t params_hash;

    const uint32_t *deps;     /* node ids */
    uint32_t ndeps;

    OutputSpec output;
    PropertyComputeFn compute;
    void *user_data;
} PropertyNode;

typedef struct PropertyCacheEntry {
    CacheKey key;
    PropertyOutput out;
    uint8_t valid;
} PropertyCacheEntry;

typedef struct PropertyCachePool {
    uint32_t max_nodes;
    uint32_t ring_size;
    uint32_t max_output_bytes;

    PropertyCacheEntry *system_entries;   /* [max_nodes] */
    PropertyCacheEntry *frame_entries;    /* [max_nodes*ring_size] */
    PropertyCacheEntry *window_entries;   /* [max_nodes*ring_size] */

    uint8_t *system_storage;              /* [max_nodes*max_output_bytes] */
    uint8_t *frame_storage;               /* [max_nodes*ring_size*max_output_bytes] */
    uint8_t *window_storage;              /* [max_nodes*ring_size*max_output_bytes] */
} PropertyCachePool;

typedef struct SchedulerConfig {
    uint32_t ring_size;
    uint32_t nthreads;

    uint32_t max_nodes;
    uint32_t max_output_bytes;

    uint32_t max_cells;
    uint32_t max_cell_items;
    uint32_t max_bins;

    uint32_t cell_dims[3];
    float cell_size[3];
    float default_cutoff;
    uint8_t pbc[3];
} SchedulerConfig;

typedef struct FrameInput {
    const float *positions_xyz; /* [natoms*3] */
    float box[3];
    double time_ps;
    uint64_t step;
} FrameInput;

typedef struct DistancePolicy {
    uint8_t periodic[3];
} DistancePolicy;

typedef void (*NeighborVisitorFn)(
    uint32_t ref_mol,
    uint32_t nbr_mol,
    float distance2,
    void *user_data
);

typedef struct EngineContext {
    const SystemTopology *system;

    FrameRing frame_ring;
    Workspace workspace;
    PropertyCachePool cache;

    PropertyNode *nodes;       /* enabled DAG nodes */
    uint32_t nnodes;

    SchedulerConfig cfg;

    /* Active execution context for node->compute(ctx, out) */
    uint64_t active_frame_index;
    uint32_t active_slot;
    uint32_t active_node_id;
} EngineContext;

enum {
    ENGINE_NODE_MOL_COM = 0,
    ENGINE_NODE_CELL_LIST = 1,
    ENGINE_NODE_DENSITY_PER_TYPE = 2,
    ENGINE_DEFAULT_NODE_COUNT = 3
};

int scheduler_init(
    EngineContext *ctx,
    const SystemTopology *system,
    const SchedulerConfig *cfg
);

void scheduler_destroy(EngineContext *ctx);

int scheduler_ingest_frame(
    EngineContext *ctx,
    uint64_t frame_index,
    const FrameInput *frame
);

int scheduler_execute_node(
    EngineContext *ctx,
    uint32_t node_id,
    uint64_t frame_index,
    PropertyOutput **out
);

int scheduler_ensure_mol_com(EngineContext *ctx, uint64_t frame_index);
int scheduler_ensure_cell_list(EngineContext *ctx, uint64_t frame_index);

int scheduler_run_density_per_type(
    EngineContext *ctx,
    uint64_t frame_index,
    const double **density_out,
    uint32_t *ntypes_out
);

int cell_list_query_neighbors(
    const FrameSlot *slot,
    uint32_t ref_mol,
    float cutoff,
    const DistancePolicy *policy,
    NeighborVisitorFn visitor,
    void *user_data
);

#ifdef __cplusplus
}
#endif

#endif /* ENGINE_TYPES_H */
