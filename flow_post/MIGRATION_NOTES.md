# Flow Postprocessing Migration Notes

## Reference methodology (unchanged physics)

The transport methodology remains the same:

1. Input trajectories of cumulative net crossings (`time_ps`, `SOL`, `NA`, `CL`)
2. Replica alignment on a common time grid (step-hold interpolation)
3. Long-time slope fit over `fit_last_ns` (default 45 ns)
4. Conversion to current, current density, and conductivity
5. Aggregation across fields/surfaces and Ohmic fits

Only the input/output layer changed.

## New reference workflow

```text
runs/rep_XX/gating_flux.csv
  -> extract_transport_input.py
  -> runs/rep_XX/transport_input.tsv
  -> reduce_counts.py
  -> fit_transport.py
```

`compile/` is no longer part of the intended default workflow.

## Center-plane mapping (fixed)

From `runs/rep_XX/gating_flux.csv`, use:

- `time_ps`
- `center_water_left`, `center_water_right`
- `center_na_left`, `center_na_right`
- `center_cl_left`, `center_cl_right`

Define:

- `center_water_net = center_water_right - center_water_left`
- `center_na_net    = center_na_right - center_na_left`
- `center_cl_net    = center_cl_right - center_cl_left`

Then cumulative:

- `SOL = cumsum(center_water_net)`
- `NA  = cumsum(center_na_net)`
- `CL  = cumsum(center_cl_net)`

and write `runs/rep_XX/transport_input.tsv` with columns exactly:

- `time_ps`
- `SOL`
- `NA`
- `CL`

## Commands

Inventory:

```bash
python3 flow_post/extract_transport_input.py /path/to/FIELD_00 --inventory-only
```

Extract normalized transport inputs:

```bash
python3 flow_post/extract_transport_input.py /path/to/FIELD_00
```

Reduce from runs:

```bash
python3 flow_post/reduce_counts.py /path/to/FIELD_00/runs --replica-glob "rep_*/transport_input.tsv"
```

Fit from runs:

```bash
python3 flow_post/fit_transport.py /path/to/FIELD_00/runs --replica-glob "rep_*/transport_input.tsv" --fit-last-ns 45
```

## Legacy compatibility (secondary)

Legacy `count_XX.dat` support remains optional:

- `extract_transport_input.py --legacy-count-outdir <dir>`
- `reduce_counts.py --legacy-indexed-counts ...`
- `fit_transport.py --legacy-indexed-counts ...`
