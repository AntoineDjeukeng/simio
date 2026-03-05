#!/usr/bin/env bash
set -euo pipefail

XTC="${1:-react_neu_9_6_20.xtc}"
FRAMES="${2:-200}"

OUTCSV="/tmp/coord_smoke.csv"
OUTLOG="/tmp/coord_smoke.log"

./build/simio_xtc_coord_x "$XTC" \
"$FRAMES" 4 0.5 5896 110 110 \
7.11 12.89 0.901 1.801 \
100 0.35 0.38 0.35 \
"$OUTCSV" > "$OUTLOG"

echo "=== INVARIANTS ==="
if grep -Eq "status=FAIL|oob=[1-9]" "$OUTLOG"; then
  echo "FAIL"
  exit 1
else
  echo "OK"
fi

echo "=== SUMMARY (outer20 vs mid60) ==="
python3 scripts/coord_summary.py "$OUTCSV"
