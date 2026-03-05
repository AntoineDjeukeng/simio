#include "engine_types.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

enum {
    ENGINE_OK = 0,
    ENGINE_ERR_ARG = -1,
    ENGINE_ERR_ALLOC = -2,
    ENGINE_ERR_CONFIG = -3,
    ENGINE_ERR_STATE = -4
};

static const uint32_t kDepsCellList[] = { ENGINE_NODE_MOL_COM };
static const uint32_t kDepsDensity[] = { ENGINE_NODE_MOL_COM, ENGINE_NODE_CELL_LIST };

static int compute_node_mol_com(EngineContext *ctx, PropertyOutput *out);
static int compute_node_cell_list(EngineContext *ctx, PropertyOutput *out);
static int compute_node_density_per_type(EngineContext *ctx, PropertyOutput *out);

static void *xcalloc(size_t n, size_t elem_size) {
    if (elem_size != 0 && n > SIZE_MAX / elem_size) {
        return NULL;
    }
    return calloc(n, elem_size);
}

static uint32_t frame_slot_index(const EngineContext *ctx, uint64_t frame_index) {
    return (uint32_t)(frame_index % ctx->frame_ring.k_slots);
}

static uint64_t window_start_for_frame(const EngineContext *ctx, uint64_t frame_index) {
    const uint64_t k = (uint64_t)ctx->frame_ring.k_slots;
    if (k == 0 || frame_index + 1 < k) {
        return 0;
    }
    return frame_index + 1 - k;
}

static int get_slot_for_frame(EngineContext *ctx, uint64_t frame_index, FrameSlot **slot_out) {
    FrameSlot *slot;
    uint32_t slot_idx;

    if (ctx == NULL || slot_out == NULL || ctx->frame_ring.k_slots == 0 || ctx->frame_ring.slots == NULL) {
        return ENGINE_ERR_ARG;
    }

    slot_idx = frame_slot_index(ctx, frame_index);
    slot = &ctx->frame_ring.slots[slot_idx];

    if (!slot->occupied || slot->frame_index != frame_index) {
        return ENGINE_ERR_STATE;
    }

    *slot_out = slot;
    return ENGINE_OK;
}

static void free_frame_ring(FrameRing *ring) {
    uint32_t i;

    if (ring == NULL || ring->slots == NULL) {
        return;
    }

    for (i = 0; i < ring->k_slots; ++i) {
        FrameSlot *slot = &ring->slots[i];
        free(slot->positions_xyz);
        free(slot->mol_com_xyz);
        free(slot->cell_list.cell_offsets);
        free(slot->cell_list.cell_ids);
        free(slot->cell_list.mol_cell);
    }

    free(ring->slots);
    ring->slots = NULL;
    ring->k_slots = 0;
    ring->frames_seen = 0;
}

static void free_workspace(Workspace *ws) {
    if (ws == NULL) {
        return;
    }

    free(ws->counts_tls);
    free(ws->fill_tls);
    free(ws->cell_totals);
    free(ws->bins_tls);
    free(ws->bins_reduced);
    memset(ws, 0, sizeof(*ws));
}

static void free_cache_pool(PropertyCachePool *cache) {
    if (cache == NULL) {
        return;
    }

    free(cache->system_entries);
    free(cache->frame_entries);
    free(cache->window_entries);
    free(cache->system_storage);
    free(cache->frame_storage);
    free(cache->window_storage);
    memset(cache, 0, sizeof(*cache));
}

static void bind_cache_storage(
    PropertyCacheEntry *entries,
    uint8_t *storage,
    size_t nentries,
    uint32_t stride
) {
    size_t i;

    for (i = 0; i < nentries; ++i) {
        entries[i].out.data = storage + (i * (size_t)stride);
        entries[i].out.bytes = 0;
        entries[i].out.capacity_bytes = (uint64_t)stride;
        entries[i].valid = 0;
        memset(&entries[i].key, 0, sizeof(entries[i].key));
    }
}

static int init_frame_ring(EngineContext *ctx) {
    uint32_t s;
    const SystemTopology *sys = ctx->system;
    const SchedulerConfig *cfg = &ctx->cfg;
    uint64_t ncells64;
    uint32_t ncells;

    if (cfg->ring_size == 0 || cfg->max_cell_items < sys->nmol) {
        return ENGINE_ERR_CONFIG;
    }

    if (cfg->cell_dims[0] == 0 || cfg->cell_dims[1] == 0 || cfg->cell_dims[2] == 0) {
        return ENGINE_ERR_CONFIG;
    }

    if (cfg->cell_size[0] <= 0.0f || cfg->cell_size[1] <= 0.0f || cfg->cell_size[2] <= 0.0f) {
        return ENGINE_ERR_CONFIG;
    }

    ncells64 = (uint64_t)cfg->cell_dims[0] * (uint64_t)cfg->cell_dims[1] * (uint64_t)cfg->cell_dims[2];
    if (ncells64 == 0 || ncells64 > cfg->max_cells) {
        return ENGINE_ERR_CONFIG;
    }
    ncells = (uint32_t)ncells64;

    ctx->frame_ring.k_slots = cfg->ring_size;
    ctx->frame_ring.frames_seen = 0;
    ctx->frame_ring.slots = (FrameSlot *)xcalloc(cfg->ring_size, sizeof(FrameSlot));
    if (ctx->frame_ring.slots == NULL) {
        return ENGINE_ERR_ALLOC;
    }

    for (s = 0; s < cfg->ring_size; ++s) {
        FrameSlot *slot = &ctx->frame_ring.slots[s];

        slot->positions_xyz = (float *)xcalloc((size_t)sys->natoms * 3U, sizeof(float));
        slot->mol_com_xyz = (float *)xcalloc((size_t)sys->nmol * 3U, sizeof(float));
        slot->cell_list.cell_offsets = (uint32_t *)xcalloc((size_t)ncells + 1U, sizeof(uint32_t));
        slot->cell_list.cell_ids = (uint32_t *)xcalloc((size_t)cfg->max_cell_items, sizeof(uint32_t));
        slot->cell_list.mol_cell = (uint32_t *)xcalloc((size_t)sys->nmol, sizeof(uint32_t));

        if (slot->positions_xyz == NULL || slot->mol_com_xyz == NULL ||
            slot->cell_list.cell_offsets == NULL || slot->cell_list.cell_ids == NULL ||
            slot->cell_list.mol_cell == NULL) {
            return ENGINE_ERR_ALLOC;
        }

        slot->occupied = 0;
        slot->mol_com_ready = 0;
        slot->mol_com_frame_index = 0;
        slot->cell_list.ready = 0;
        slot->cell_list.built_frame_index = 0;
        slot->cell_list.nitems = 0;

        slot->cell_list.nx = cfg->cell_dims[0];
        slot->cell_list.ny = cfg->cell_dims[1];
        slot->cell_list.nz = cfg->cell_dims[2];
        slot->cell_list.ncells = ncells;
        slot->cell_list.max_items = cfg->max_cell_items;
        slot->cell_list.cell_size[0] = cfg->cell_size[0];
        slot->cell_list.cell_size[1] = cfg->cell_size[1];
        slot->cell_list.cell_size[2] = cfg->cell_size[2];
        slot->cell_list.inv_cell_size[0] = 1.0f / cfg->cell_size[0];
        slot->cell_list.inv_cell_size[1] = 1.0f / cfg->cell_size[1];
        slot->cell_list.inv_cell_size[2] = 1.0f / cfg->cell_size[2];
        slot->cell_list.pbc[0] = cfg->pbc[0];
        slot->cell_list.pbc[1] = cfg->pbc[1];
        slot->cell_list.pbc[2] = cfg->pbc[2];
    }

    return ENGINE_OK;
}

static int init_workspace(EngineContext *ctx) {
    Workspace *ws = &ctx->workspace;
    uint32_t nthreads = ctx->cfg.nthreads;
    uint32_t nbins = ctx->cfg.max_bins;
    uint32_t ncells = ctx->frame_ring.slots[0].cell_list.ncells;

    if (nthreads == 0 || nbins == 0 || nbins < ctx->system->ntypes) {
        return ENGINE_ERR_CONFIG;
    }

    ws->nthreads = nthreads;
    ws->nbins = nbins;
    ws->ncells = ncells;

    ws->counts_tls = (uint32_t *)xcalloc((size_t)nthreads * ncells, sizeof(uint32_t));
    ws->fill_tls = (uint32_t *)xcalloc((size_t)nthreads * ncells, sizeof(uint32_t));
    ws->cell_totals = (uint32_t *)xcalloc(ncells, sizeof(uint32_t));
    ws->bins_tls = (double *)xcalloc((size_t)nthreads * nbins, sizeof(double));
    ws->bins_reduced = (double *)xcalloc(nbins, sizeof(double));

    if (ws->counts_tls == NULL || ws->fill_tls == NULL || ws->cell_totals == NULL ||
        ws->bins_tls == NULL || ws->bins_reduced == NULL) {
        return ENGINE_ERR_ALLOC;
    }

    return ENGINE_OK;
}

static int init_cache_pool(EngineContext *ctx) {
    PropertyCachePool *cache = &ctx->cache;
    uint32_t max_nodes = ctx->cfg.max_nodes;
    uint32_t ring = ctx->cfg.ring_size;
    uint32_t max_output_bytes = ctx->cfg.max_output_bytes;
    size_t n_frame_entries;
    size_t n_window_entries;
    size_t n_system_entries;

    if (max_nodes == 0 || ring == 0 || max_output_bytes == 0) {
        return ENGINE_ERR_CONFIG;
    }

    n_system_entries = (size_t)max_nodes;
    n_frame_entries = (size_t)max_nodes * ring;
    n_window_entries = (size_t)max_nodes * ring;

    cache->max_nodes = max_nodes;
    cache->ring_size = ring;
    cache->max_output_bytes = max_output_bytes;

    cache->system_entries = (PropertyCacheEntry *)xcalloc(n_system_entries, sizeof(PropertyCacheEntry));
    cache->frame_entries = (PropertyCacheEntry *)xcalloc(n_frame_entries, sizeof(PropertyCacheEntry));
    cache->window_entries = (PropertyCacheEntry *)xcalloc(n_window_entries, sizeof(PropertyCacheEntry));
    cache->system_storage = (uint8_t *)xcalloc(n_system_entries, max_output_bytes);
    cache->frame_storage = (uint8_t *)xcalloc(n_frame_entries, max_output_bytes);
    cache->window_storage = (uint8_t *)xcalloc(n_window_entries, max_output_bytes);

    if (cache->system_entries == NULL || cache->frame_entries == NULL || cache->window_entries == NULL ||
        cache->system_storage == NULL || cache->frame_storage == NULL || cache->window_storage == NULL) {
        return ENGINE_ERR_ALLOC;
    }

    bind_cache_storage(cache->system_entries, cache->system_storage, n_system_entries, max_output_bytes);
    bind_cache_storage(cache->frame_entries, cache->frame_storage, n_frame_entries, max_output_bytes);
    bind_cache_storage(cache->window_entries, cache->window_storage, n_window_entries, max_output_bytes);

    return ENGINE_OK;
}

static uint64_t hash_mix(uint64_t seed, uint64_t v) {
    seed ^= v + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
}

static int register_default_nodes(EngineContext *ctx) {
    PropertyNode *nodes;
    uint64_t h;

    if (ctx->cfg.max_nodes < ENGINE_DEFAULT_NODE_COUNT) {
        return ENGINE_ERR_CONFIG;
    }

    nodes = (PropertyNode *)xcalloc(ENGINE_DEFAULT_NODE_COUNT, sizeof(PropertyNode));
    if (nodes == NULL) {
        return ENGINE_ERR_ALLOC;
    }

    nodes[ENGINE_NODE_MOL_COM].node_id = ENGINE_NODE_MOL_COM;
    nodes[ENGINE_NODE_MOL_COM].name = "MolCOM";
    nodes[ENGINE_NODE_MOL_COM].scope = PROPERTY_SCOPE_FRAME;
    nodes[ENGINE_NODE_MOL_COM].params_hash = 0x1001ULL;
    nodes[ENGINE_NODE_MOL_COM].deps = NULL;
    nodes[ENGINE_NODE_MOL_COM].ndeps = 0;
    nodes[ENGINE_NODE_MOL_COM].output.label = "mol_com_xyz";
    nodes[ENGINE_NODE_MOL_COM].output.scope = PROPERTY_SCOPE_FRAME;
    nodes[ENGINE_NODE_MOL_COM].output.kind = OUTPUT_KIND_ARRAY_F32;
    nodes[ENGINE_NODE_MOL_COM].output.required_bytes = (uint64_t)ctx->system->nmol * 3ULL * sizeof(float);
    nodes[ENGINE_NODE_MOL_COM].compute = compute_node_mol_com;

    h = 0xC311ULL;
    h = hash_mix(h, ctx->cfg.cell_dims[0]);
    h = hash_mix(h, ctx->cfg.cell_dims[1]);
    h = hash_mix(h, ctx->cfg.cell_dims[2]);
    h = hash_mix(h, (uint64_t)ctx->cfg.pbc[0] | ((uint64_t)ctx->cfg.pbc[1] << 8) | ((uint64_t)ctx->cfg.pbc[2] << 16));

    nodes[ENGINE_NODE_CELL_LIST].node_id = ENGINE_NODE_CELL_LIST;
    nodes[ENGINE_NODE_CELL_LIST].name = "CellList";
    nodes[ENGINE_NODE_CELL_LIST].scope = PROPERTY_SCOPE_FRAME;
    nodes[ENGINE_NODE_CELL_LIST].params_hash = h;
    nodes[ENGINE_NODE_CELL_LIST].deps = kDepsCellList;
    nodes[ENGINE_NODE_CELL_LIST].ndeps = 1;
    nodes[ENGINE_NODE_CELL_LIST].output.label = "cell_list";
    nodes[ENGINE_NODE_CELL_LIST].output.scope = PROPERTY_SCOPE_FRAME;
    nodes[ENGINE_NODE_CELL_LIST].output.kind = OUTPUT_KIND_OPAQUE;
    nodes[ENGINE_NODE_CELL_LIST].output.required_bytes = sizeof(CellList);
    nodes[ENGINE_NODE_CELL_LIST].compute = compute_node_cell_list;

    h = 0xD351ULL;
    h = hash_mix(h, ctx->system->ntypes);
    h = hash_mix(h, ctx->workspace.nbins);
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].node_id = ENGINE_NODE_DENSITY_PER_TYPE;
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].name = "DensityPerType";
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].scope = PROPERTY_SCOPE_FRAME;
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].params_hash = h;
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].deps = kDepsDensity;
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].ndeps = 2;
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].output.label = "density_per_type";
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].output.scope = PROPERTY_SCOPE_FRAME;
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].output.kind = OUTPUT_KIND_ARRAY_F64;
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].output.required_bytes = (uint64_t)ctx->system->ntypes * sizeof(double);
    nodes[ENGINE_NODE_DENSITY_PER_TYPE].compute = compute_node_density_per_type;

    ctx->nodes = nodes;
    ctx->nnodes = ENGINE_DEFAULT_NODE_COUNT;
    return ENGINE_OK;
}

static void invalidate_overwritten_slot(EngineContext *ctx, uint32_t slot) {
    uint32_t node;

    for (node = 0; node < ctx->cache.max_nodes; ++node) {
        size_t idx = (size_t)node * ctx->cache.ring_size + slot;
        ctx->cache.frame_entries[idx].valid = 0;
        ctx->cache.window_entries[idx].valid = 0;
    }
}

static PropertyCacheEntry *cache_entry_for(
    EngineContext *ctx,
    const PropertyNode *node,
    uint64_t frame_index
) {
    uint32_t slot;

    if (node->node_id >= ctx->cache.max_nodes) {
        return NULL;
    }

    switch (node->scope) {
        case PROPERTY_SCOPE_SYSTEM:
            return &ctx->cache.system_entries[node->node_id];
        case PROPERTY_SCOPE_FRAME:
            slot = frame_slot_index(ctx, frame_index);
            return &ctx->cache.frame_entries[(size_t)node->node_id * ctx->cache.ring_size + slot];
        case PROPERTY_SCOPE_WINDOW:
            slot = (uint32_t)(window_start_for_frame(ctx, frame_index) % ctx->cache.ring_size);
            return &ctx->cache.window_entries[(size_t)node->node_id * ctx->cache.ring_size + slot];
        default:
            return NULL;
    }
}

static int cache_key_matches(
    const EngineContext *ctx,
    const PropertyNode *node,
    const PropertyCacheEntry *entry,
    uint64_t frame_index
) {
    uint64_t win_start;
    uint64_t win_end;

    if (!entry->valid) {
        return 0;
    }

    if (entry->key.node_id != node->node_id ||
        entry->key.scope != node->scope ||
        entry->key.params_hash != node->params_hash) {
        return 0;
    }

    if (node->scope == PROPERTY_SCOPE_SYSTEM) {
        return 1;
    }

    if (node->scope == PROPERTY_SCOPE_FRAME) {
        return entry->key.frame_index == frame_index;
    }

    win_start = window_start_for_frame(ctx, frame_index);
    win_end = win_start + ctx->frame_ring.k_slots - 1ULL;
    return entry->key.window_start == win_start && entry->key.window_end == win_end;
}

int scheduler_init(EngineContext *ctx, const SystemTopology *system, const SchedulerConfig *cfg) {
    int rc;

    if (ctx == NULL || system == NULL || cfg == NULL || system->mol_index == NULL ||
        system->mol_atom_offsets == NULL) {
        return ENGINE_ERR_ARG;
    }

    memset(ctx, 0, sizeof(*ctx));
    ctx->system = system;
    ctx->cfg = *cfg;

    rc = init_frame_ring(ctx);
    if (rc != ENGINE_OK) {
        scheduler_destroy(ctx);
        return rc;
    }

    rc = init_workspace(ctx);
    if (rc != ENGINE_OK) {
        scheduler_destroy(ctx);
        return rc;
    }

    rc = init_cache_pool(ctx);
    if (rc != ENGINE_OK) {
        scheduler_destroy(ctx);
        return rc;
    }

    rc = register_default_nodes(ctx);
    if (rc != ENGINE_OK) {
        scheduler_destroy(ctx);
        return rc;
    }

    return ENGINE_OK;
}

void scheduler_destroy(EngineContext *ctx) {
    if (ctx == NULL) {
        return;
    }

    free(ctx->nodes);
    ctx->nodes = NULL;
    ctx->nnodes = 0;

    free_cache_pool(&ctx->cache);
    free_workspace(&ctx->workspace);
    free_frame_ring(&ctx->frame_ring);

    ctx->system = NULL;
    memset(&ctx->cfg, 0, sizeof(ctx->cfg));
    ctx->active_frame_index = 0;
    ctx->active_slot = 0;
    ctx->active_node_id = 0;
}

int scheduler_ingest_frame(EngineContext *ctx, uint64_t frame_index, const FrameInput *frame) {
    FrameSlot *slot;
    uint32_t slot_idx;
    size_t nxyz;

    if (ctx == NULL || frame == NULL || frame->positions_xyz == NULL || ctx->system == NULL) {
        return ENGINE_ERR_ARG;
    }

    slot_idx = frame_slot_index(ctx, frame_index);
    slot = &ctx->frame_ring.slots[slot_idx];
    nxyz = (size_t)ctx->system->natoms * 3U;

    memcpy(slot->positions_xyz, frame->positions_xyz, nxyz * sizeof(float));
    slot->box[0] = frame->box[0];
    slot->box[1] = frame->box[1];
    slot->box[2] = frame->box[2];
    slot->time_ps = frame->time_ps;
    slot->step = frame->step;
    slot->frame_index = frame_index;
    slot->occupied = 1;

    slot->mol_com_ready = 0;
    slot->cell_list.ready = 0;
    slot->mol_com_frame_index = 0;
    slot->cell_list.built_frame_index = 0;
    slot->cell_list.nitems = 0;

    ctx->frame_ring.frames_seen = frame_index + 1;

    invalidate_overwritten_slot(ctx, slot_idx);
    return ENGINE_OK;
}

int scheduler_ensure_mol_com(EngineContext *ctx, uint64_t frame_index) {
    FrameSlot *slot;
    const SystemTopology *sys;
    uint32_t mol;
    int rc;

    if (ctx == NULL || ctx->system == NULL) {
        return ENGINE_ERR_ARG;
    }

    rc = get_slot_for_frame(ctx, frame_index, &slot);
    if (rc != ENGINE_OK) {
        return rc;
    }

    if (slot->mol_com_ready && slot->mol_com_frame_index == frame_index) {
        return ENGINE_OK;
    }

    sys = ctx->system;
    for (mol = 0; mol < sys->nmol; ++mol) {
        const uint32_t a0 = sys->mol_atom_offsets[mol];
        const uint32_t a1 = sys->mol_atom_offsets[mol + 1];
        float sx = 0.0f;
        float sy = 0.0f;
        float sz = 0.0f;
        float inv = 0.0f;
        uint32_t a;
        uint32_t nat;

        if (a1 <= a0 || a1 > sys->natoms) {
            return ENGINE_ERR_STATE;
        }

        nat = a1 - a0;
        inv = 1.0f / (float)nat;

        for (a = a0; a < a1; ++a) {
            const float *r = &slot->positions_xyz[(size_t)a * 3U];
            sx += r[0];
            sy += r[1];
            sz += r[2];
        }

        slot->mol_com_xyz[(size_t)mol * 3U + 0U] = sx * inv;
        slot->mol_com_xyz[(size_t)mol * 3U + 1U] = sy * inv;
        slot->mol_com_xyz[(size_t)mol * 3U + 2U] = sz * inv;
    }

    slot->mol_com_ready = 1;
    slot->mol_com_frame_index = frame_index;
    return ENGINE_OK;
}

static uint32_t wrap_cell_coord(int v, uint32_t n) {
    int m = (int)n;
    int r = v % m;
    if (r < 0) {
        r += m;
    }
    return (uint32_t)r;
}

static uint32_t flatten_cell(const CellList *cl, uint32_t ix, uint32_t iy, uint32_t iz) {
    return (iz * cl->ny + iy) * cl->nx + ix;
}

static uint32_t cell_index_from_xyz(const CellList *cl, const float box[3], const float xyz[3]) {
    uint32_t out[3];
    uint32_t n[3];
    uint32_t a;

    n[0] = cl->nx;
    n[1] = cl->ny;
    n[2] = cl->nz;

    for (a = 0; a < 3; ++a) {
        float c = xyz[a];
        float span = cl->cell_size[a] * (float)n[a];
        uint32_t ci;

        if (cl->pbc[a] && box[a] > 0.0f) {
            c = fmodf(c, box[a]);
            if (c < 0.0f) {
                c += box[a];
            }
        }

        if (c < 0.0f) {
            c = 0.0f;
        }
        if (c >= span) {
            c = span - 1e-6f;
        }

        ci = (uint32_t)floorf(c * cl->inv_cell_size[a]);
        if (ci >= n[a]) {
            ci = n[a] - 1U;
        }
        out[a] = ci;
    }

    return flatten_cell(cl, out[0], out[1], out[2]);
}

int scheduler_ensure_cell_list(EngineContext *ctx, uint64_t frame_index) {
    FrameSlot *slot;
    CellList *cl;
    Workspace *ws;
    const MoleculeIndex *idx;
    uint32_t nthreads;
    uint32_t ncells;
    uint32_t tid;
    uint32_t cell;
    int rc;

    if (ctx == NULL || ctx->system == NULL || ctx->system->mol_index == NULL) {
        return ENGINE_ERR_ARG;
    }

    rc = scheduler_ensure_mol_com(ctx, frame_index);
    if (rc != ENGINE_OK) {
        return rc;
    }

    rc = get_slot_for_frame(ctx, frame_index, &slot);
    if (rc != ENGINE_OK) {
        return rc;
    }

    cl = &slot->cell_list;
    ws = &ctx->workspace;
    idx = ctx->system->mol_index;
    nthreads = ws->nthreads;
    ncells = cl->ncells;

    if (cl->ready && cl->built_frame_index == frame_index) {
        return ENGINE_OK;
    }

    memset(ws->counts_tls, 0, (size_t)nthreads * ncells * sizeof(uint32_t));
    memset(ws->fill_tls, 0, (size_t)nthreads * ncells * sizeof(uint32_t));
    memset(ws->cell_totals, 0, (size_t)ncells * sizeof(uint32_t));
    memset(cl->cell_offsets, 0, ((size_t)ncells + 1U) * sizeof(uint32_t));

    /*
     * Phase A: per-thread, per-cell counts with canonical per-type strided loop.
     * In production this loop nest runs inside a real parallel region.
     */
    for (tid = 0; tid < nthreads; ++tid) {
        uint32_t t;
        for (t = 0; t < idx->ntypes; ++t) {
            const uint32_t n_type = idx->type_counts[t];
            const uint32_t *mol_ids = idx->type_mol_ids[t];
            uint32_t k;

            for (k = tid; k < n_type; k += nthreads) {
                const uint32_t mol = mol_ids[k];
                const float *r = &slot->mol_com_xyz[(size_t)mol * 3U];
                const uint32_t c = cell_index_from_xyz(cl, slot->box, r);
                if (mol >= ctx->system->nmol) {
                    return ENGINE_ERR_STATE;
                }
                ws->counts_tls[(size_t)tid * ncells + c] += 1U;
                cl->mol_cell[mol] = c;
            }
        }
    }

    /* Prefix sum over globally reduced per-cell totals. */
    cl->cell_offsets[0] = 0;
    for (cell = 0; cell < ncells; ++cell) {
        uint32_t total = 0;
        for (tid = 0; tid < nthreads; ++tid) {
            total += ws->counts_tls[(size_t)tid * ncells + cell];
        }
        ws->cell_totals[cell] = total;
        cl->cell_offsets[cell + 1U] = cl->cell_offsets[cell] + total;
    }

    cl->nitems = cl->cell_offsets[ncells];
    if (cl->nitems > cl->max_items) {
        return ENGINE_ERR_STATE;
    }

    /*
     * Precompute thread-private base offsets per cell.
     * fill_tls[tid][cell] becomes a thread-local write cursor.
     */
    for (cell = 0; cell < ncells; ++cell) {
        uint32_t base = cl->cell_offsets[cell];
        for (tid = 0; tid < nthreads; ++tid) {
            size_t at = (size_t)tid * ncells + cell;
            ws->fill_tls[at] = base;
            base += ws->counts_tls[at];
        }
    }

    /*
     * Phase B: same ownership pattern, fill flat cell_ids with no atomics.
     * In production this loop nest runs inside a real parallel region.
     */
    for (tid = 0; tid < nthreads; ++tid) {
        uint32_t t;
        for (t = 0; t < idx->ntypes; ++t) {
            const uint32_t n_type = idx->type_counts[t];
            const uint32_t *mol_ids = idx->type_mol_ids[t];
            uint32_t k;

            for (k = tid; k < n_type; k += nthreads) {
                const uint32_t mol = mol_ids[k];
                const uint32_t c = cl->mol_cell[mol];
                const size_t cursor_idx = (size_t)tid * ncells + c;
                const uint32_t write_idx = ws->fill_tls[cursor_idx]++;
                if (mol >= ctx->system->nmol) {
                    return ENGINE_ERR_STATE;
                }
                cl->cell_ids[write_idx] = mol;
            }
        }
    }

    cl->ready = 1;
    cl->built_frame_index = frame_index;
    return ENGINE_OK;
}

static float min_image_delta(float d, float box_len, uint8_t periodic) {
    if (!periodic || box_len <= 0.0f) {
        return d;
    }
    if (d > 0.5f * box_len) {
        d -= box_len;
    } else if (d < -0.5f * box_len) {
        d += box_len;
    }
    return d;
}

static void unpack_cell(const CellList *cl, uint32_t c, int *ix, int *iy, int *iz) {
    *ix = (int)(c % cl->nx);
    c /= cl->nx;
    *iy = (int)(c % cl->ny);
    c /= cl->ny;
    *iz = (int)c;
}

int cell_list_query_neighbors(
    const FrameSlot *slot,
    uint32_t ref_mol,
    float cutoff,
    const DistancePolicy *policy,
    NeighborVisitorFn visitor,
    void *user_data
) {
    const CellList *cl;
    const float *ri;
    float cutoff2;
    uint8_t periodic[3];
    int cx;
    int cy;
    int cz;
    int dz;
    int found = 0;

    if (slot == NULL || visitor == NULL || !slot->mol_com_ready || !slot->cell_list.ready) {
        return ENGINE_ERR_ARG;
    }

    cl = &slot->cell_list;
    if (ref_mol >= cl->nitems) {
        return ENGINE_ERR_ARG;
    }

    periodic[0] = policy ? policy->periodic[0] : cl->pbc[0];
    periodic[1] = policy ? policy->periodic[1] : cl->pbc[1];
    periodic[2] = policy ? policy->periodic[2] : cl->pbc[2];

    cutoff2 = cutoff * cutoff;
    ri = &slot->mol_com_xyz[(size_t)ref_mol * 3U];
    unpack_cell(cl, cl->mol_cell[ref_mol], &cx, &cy, &cz);

    for (dz = -1; dz <= 1; ++dz) {
        int dy;
        for (dy = -1; dy <= 1; ++dy) {
            int dx;
            for (dx = -1; dx <= 1; ++dx) {
                int nx = cx + dx;
                int ny = cy + dy;
                int nz = cz + dz;
                uint32_t ncell;
                uint32_t begin;
                uint32_t end;
                uint32_t p;

                if (nx < 0 || nx >= (int)cl->nx) {
                    if (!periodic[0]) {
                        continue;
                    }
                    nx = (int)wrap_cell_coord(nx, cl->nx);
                }
                if (ny < 0 || ny >= (int)cl->ny) {
                    if (!periodic[1]) {
                        continue;
                    }
                    ny = (int)wrap_cell_coord(ny, cl->ny);
                }
                if (nz < 0 || nz >= (int)cl->nz) {
                    if (!periodic[2]) {
                        continue;
                    }
                    nz = (int)wrap_cell_coord(nz, cl->nz);
                }

                ncell = flatten_cell(cl, (uint32_t)nx, (uint32_t)ny, (uint32_t)nz);
                begin = cl->cell_offsets[ncell];
                end = cl->cell_offsets[ncell + 1U];

                for (p = begin; p < end; ++p) {
                    const uint32_t nbr = cl->cell_ids[p];
                    const float *rj;
                    float dxr;
                    float dyr;
                    float dzr;
                    float r2;

                    if (nbr == ref_mol) {
                        continue;
                    }

                    rj = &slot->mol_com_xyz[(size_t)nbr * 3U];
                    dxr = min_image_delta(rj[0] - ri[0], slot->box[0], periodic[0]);
                    dyr = min_image_delta(rj[1] - ri[1], slot->box[1], periodic[1]);
                    dzr = min_image_delta(rj[2] - ri[2], slot->box[2], periodic[2]);
                    r2 = dxr * dxr + dyr * dyr + dzr * dzr;

                    if (r2 <= cutoff2) {
                        visitor(ref_mol, nbr, r2, user_data);
                        ++found;
                    }
                }
            }
        }
    }

    return found;
}

static int compute_node_mol_com(EngineContext *ctx, PropertyOutput *out) {
    FrameSlot *slot;
    int rc;

    rc = get_slot_for_frame(ctx, ctx->active_frame_index, &slot);
    if (rc != ENGINE_OK) {
        return rc;
    }

    rc = scheduler_ensure_mol_com(ctx, ctx->active_frame_index);
    if (rc != ENGINE_OK) {
        return rc;
    }

    out->data = slot->mol_com_xyz;
    out->bytes = (uint64_t)ctx->system->nmol * 3ULL * sizeof(float);
    out->capacity_bytes = out->bytes;
    return ENGINE_OK;
}

static int compute_node_cell_list(EngineContext *ctx, PropertyOutput *out) {
    FrameSlot *slot;
    int rc;

    rc = get_slot_for_frame(ctx, ctx->active_frame_index, &slot);
    if (rc != ENGINE_OK) {
        return rc;
    }

    rc = scheduler_ensure_cell_list(ctx, ctx->active_frame_index);
    if (rc != ENGINE_OK) {
        return rc;
    }

    out->data = &slot->cell_list;
    out->bytes = sizeof(CellList);
    out->capacity_bytes = sizeof(CellList);
    return ENGINE_OK;
}

static int compute_node_density_per_type(EngineContext *ctx, PropertyOutput *out) {
    FrameSlot *slot;
    const MoleculeIndex *idx;
    Workspace *ws;
    uint32_t nthreads;
    uint32_t tid;
    uint32_t t;
    uint64_t required;
    double *dst;
    int rc;

    rc = get_slot_for_frame(ctx, ctx->active_frame_index, &slot);
    if (rc != ENGINE_OK) {
        return rc;
    }

    idx = ctx->system->mol_index;
    ws = &ctx->workspace;
    nthreads = ws->nthreads;

    if (ws->nbins < idx->ntypes || out->data == NULL) {
        return ENGINE_ERR_STATE;
    }

    memset(ws->bins_tls, 0, (size_t)nthreads * ws->nbins * sizeof(double));
    memset(ws->bins_reduced, 0, (size_t)ws->nbins * sizeof(double));

    /*
     * Canonical parallel ownership pattern:
     * for each type t, thread tid handles indices k=tid, tid+nthreads, ...
     * This is written as explicit tid loops in the skeleton.
     */
    for (tid = 0; tid < nthreads; ++tid) {
        for (t = 0; t < idx->ntypes; ++t) {
            const uint32_t n_type = idx->type_counts[t];
            const uint32_t *mol_ids = idx->type_mol_ids[t];
            uint32_t k;

            for (k = tid; k < n_type; k += nthreads) {
                const uint32_t mol = mol_ids[k];
                if (mol >= ctx->system->nmol) {
                    return ENGINE_ERR_STATE;
                }
                (void)mol;
                ws->bins_tls[(size_t)tid * ws->nbins + t] += 1.0;
            }
        }
    }

    /* Reduction step after a conceptual barrier. */
    for (tid = 0; tid < nthreads; ++tid) {
        for (t = 0; t < idx->ntypes; ++t) {
            ws->bins_reduced[t] += ws->bins_tls[(size_t)tid * ws->nbins + t];
        }
    }

    required = (uint64_t)idx->ntypes * sizeof(double);
    if (out->capacity_bytes < required) {
        return ENGINE_ERR_CONFIG;
    }

    dst = (double *)out->data;
    memcpy(dst, ws->bins_reduced, (size_t)required);
    out->bytes = required;
    return ENGINE_OK;
}

int scheduler_execute_node(
    EngineContext *ctx,
    uint32_t node_id,
    uint64_t frame_index,
    PropertyOutput **out
) {
    PropertyNode *node;
    PropertyCacheEntry *entry;
    uint32_t d;
    int rc;

    if (ctx == NULL || node_id >= ctx->nnodes) {
        return ENGINE_ERR_ARG;
    }

    node = &ctx->nodes[node_id];
    entry = cache_entry_for(ctx, node, frame_index);
    if (entry == NULL) {
        return ENGINE_ERR_CONFIG;
    }

    if (cache_key_matches(ctx, node, entry, frame_index)) {
        if (out != NULL) {
            *out = &entry->out;
        }
        return ENGINE_OK;
    }

    for (d = 0; d < node->ndeps; ++d) {
        if (node->deps[d] >= ctx->nnodes) {
            return ENGINE_ERR_CONFIG;
        }
        rc = scheduler_execute_node(ctx, node->deps[d], frame_index, NULL);
        if (rc != ENGINE_OK) {
            return rc;
        }
    }

    ctx->active_frame_index = frame_index;
    ctx->active_slot = frame_slot_index(ctx, frame_index);
    ctx->active_node_id = node_id;

    entry->out.bytes = 0;
    rc = node->compute(ctx, &entry->out);
    if (rc != ENGINE_OK) {
        return rc;
    }

    entry->key.node_id = node->node_id;
    entry->key.scope = node->scope;
    entry->key.params_hash = node->params_hash;
    entry->key.frame_index = frame_index;
    entry->key.window_start = window_start_for_frame(ctx, frame_index);
    entry->key.window_end = entry->key.window_start + ctx->frame_ring.k_slots - 1ULL;
    entry->valid = 1;

    if (out != NULL) {
        *out = &entry->out;
    }
    return ENGINE_OK;
}

int scheduler_run_density_per_type(
    EngineContext *ctx,
    uint64_t frame_index,
    const double **density_out,
    uint32_t *ntypes_out
) {
    PropertyOutput *out = NULL;
    int rc = scheduler_execute_node(ctx, ENGINE_NODE_DENSITY_PER_TYPE, frame_index, &out);
    if (rc != ENGINE_OK) {
        return rc;
    }

    if (density_out != NULL) {
        *density_out = (const double *)out->data;
    }
    if (ntypes_out != NULL) {
        *ntypes_out = ctx->system->ntypes;
    }
    return ENGINE_OK;
}
