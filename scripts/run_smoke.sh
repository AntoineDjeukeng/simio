#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <trajectory.xtc> <out_dir> [threads=1] [max_frames=200] [frame_begin=0] [frame_end=200]"
  exit 2
fi

XTC="$1"
OUT_DIR="$2"
THREADS="${3:-1}"
MAX_FRAMES="${4:-200}"
FRAME_BEGIN="${5:-0}"
FRAME_END="${6:-200}"

# Defaults match the canonical CLI usage printed by the binary.
GRID_CELL_NM=0.5
NSOL=5896
NNA=110
NCL=110
XMIN=7.11
XMAX=12.89
ZMIN=0.901
ZMAX=1.801
NX=100
R_CW=0.35
R_AW=0.38
R_OO=0.35
GATING_SEL=all
JUMP_KEEP_FRAMES=50

mkdir -p "$OUT_DIR"

./bin/simio_xtc_all_properties \
  "$XTC" \
  "$MAX_FRAMES" "$THREADS" \
  "$GRID_CELL_NM" "$NSOL" "$NNA" "$NCL" \
  "$XMIN" "$XMAX" \
  "$ZMIN" "$ZMAX" \
  "$NX" \
  "$R_CW" "$R_AW" "$R_OO" \
  "$GATING_SEL" \
  "$OUT_DIR" \
  "$JUMP_KEEP_FRAMES" \
  "$FRAME_BEGIN" "$FRAME_END" \
  > "$OUT_DIR/run.log" 2>&1

echo "OK: wrote outputs to $OUT_DIR"
