#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# USER SETTINGS
# -----------------------------

SIMIO_BIN="/home/antoine/From_milan/compute/simio_reorg/bin/simio_xtc_all_properties"
SIMIO_ROOT="$(cd "$(dirname "$SIMIO_BIN")/.." && pwd)"

RADICAL="${RADICAL:-react_neu_9_6}"
XTC_DIR="${XTC_DIR:-.}"
OUT_ROOT="${OUT_ROOT:-runs}"
CONFIG_DIR="${CONFIG_DIR:-configs}"
NOTEBOOK="$SIMIO_ROOT/scripts/simio_single_run_report.ipynb"

BASE_CONFIG="$SIMIO_ROOT/scripts/config_base.json"
TOPO_JSON="${TOPO_JSON:-topo_setup.json}"

REP_START="${REP_START:-1}"
REP_END="${REP_END:-20}"

# -----------------------------
# Required files
# -----------------------------
if [[ ! -x "$SIMIO_BIN" ]]; then
    echo "[error] missing executable: $SIMIO_BIN" >&2
    echo "Run make in $SIMIO_ROOT first." >&2
    exit 2
fi

if [[ ! -f "$BASE_CONFIG" ]]; then
    echo "[error] missing base config: $BASE_CONFIG" >&2
    exit 2
fi

if [[ ! -f "$TOPO_JSON" ]]; then
    echo "[error] missing topology JSON: $TOPO_JSON" >&2
    echo "Run this script from the directory containing topo_setup.json, or set TOPO_JSON=/absolute/path/topo_setup.json." >&2
    exit 2
fi

mkdir -p "$OUT_ROOT" "$CONFIG_DIR"

# -----------------------------
# Loop over replicas
# -----------------------------
for rep in $(seq -w "$REP_START" "$REP_END"); do

    xtc="${XTC_DIR}/${RADICAL}_${rep}.xtc"
    outdir="${OUT_ROOT}/rep_${rep}"
    cfg="${CONFIG_DIR}/config_rep_${rep}.json"
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
      -e "s|__TOPOLOGY_JSON__|$TOPO_JSON|g" \
      -e "s|__OUTDIR__|$outdir|g" \
      "$BASE_CONFIG" > "$cfg"

    # run analysis

    "$SIMIO_BIN" --config "$cfg" > "$logfile" 2>&1
    cp "$NOTEBOOK" "$outdir/simio_single_run_report.ipynb"

    echo "[done] rep $rep finished"

done
