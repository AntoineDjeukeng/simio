#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# USER SETTINGS
# -----------------------------

SIMIO_BIN="/home/antoine/From_milan/compute/simio_reorg/bin/simio_xtc_all_properties"
TOPO_BIN="/home/antoine/From_milan/compute/simio_reorg/bin/simio_mini_gro_topology"
SIMIO_ROOT="$(cd "$(dirname "$SIMIO_BIN")/.." && pwd)"

RADICAL="react_neu_9_6"
GRO="${RADICAL}_20.gro"
NOTEBOOK="$SIMIO_ROOT/scripts/simio_single_run_report.ipynb"

# BASE_CONFIG="config_base.json"
BASE_CONFIG="$SIMIO_ROOT/scripts/config_base.json"
TOPO_JSON="topo_setup.json"

SOLUTION_NAMES="SOL,NA,CL"
CHANNEL_AXIS=2
CHANNEL_MIDDLE=1.35

REP_START=1
REP_END=20

# -----------------------------
# Generate topology JSON once
# -----------------------------
if [[ ! -f "$TOPO_JSON" ]]; then
    echo "[info] generating topology summary"
    "$TOPO_BIN" "$GRO" "$SOLUTION_NAMES" "$CHANNEL_AXIS" "$CHANNEL_MIDDLE" "$TOPO_JSON"
fi

mkdir -p runs configs

# -----------------------------
# Loop over replicas
# -----------------------------
for rep in $(seq -w "$REP_START" "$REP_END"); do

    xtc="${RADICAL}_${rep}.xtc"
    outdir="runs/rep_${rep}"
    cfg="configs/config_rep_${rep}.json"
    logfile="${outdir}/run.log"

    echo
    echo "-----------------------------"
    echo "Replica $rep"
    echo "xtc: $xtc"
    echo "out: $outdir"
    echo "-----------------------------"

    if [[ ! -f "$xtc" ]]; then
        echo "[warning] missing xtc: $xtc"
        continue
    fi

    mkdir -p "$outdir"

    # build config for this replica
    sed \
      -e "s|__XTC__|$xtc|g" \
      -e "s|__OUTDIR__|$outdir|g" \
      "$BASE_CONFIG" > "$cfg"

    # run analysis

    "$SIMIO_BIN" "$cfg" > "$logfile" 2>&1
    cp "$NOTEBOOK" "$outdir/simio_single_run_report.ipynb"

    echo "[done] rep $rep finished"

done