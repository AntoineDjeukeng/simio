#!/bin/bash
#PBS -N Flow_CDI_L1
#PBS -m e
#PBS -l nodes=1:ppn=24:mp4
#PBS -j oe

set -u

L_DIR="/data/antoine/Flow_CDI/single_H_7/L_1"

SIMIO_BIN="/home/antoine/simio_reorg1/bin/simio_xtc_all_properties"
BASE_CONFIG="/home/antoine/simio_reorg1/scripts/config_base.json"

CHARGES="negative neutral positive"
FIELDS="FIELD_00 FIELD_01 FIELD_02 FIELD_03"
REP_START=1
REP_END=20

cd "$L_DIR" || exit 1

for charge in $CHARGES; do
    charge_dir="$L_DIR/$charge"
    topo="$charge_dir/mini_topology.json"

    if [[ ! -f "$topo" ]]; then
        echo "[skip] missing topology: $topo"
        continue
    fi

    for field in $FIELDS; do
        field_dir="$charge_dir/$field"
        if [[ ! -d "$field_dir" ]]; then
            echo "[skip] missing field dir: $field_dir"
            continue
        fi

        if [[ -d "$field_dir/new_runs" ]]; then
            echo "[clean] removing existing: $field_dir/new_runs"
            rm -rf "$field_dir/new_runs"
        fi
        mkdir -p "$field_dir/new_runs"

        for rep_num in $(seq "$REP_START" "$REP_END"); do
            rep="$(printf "%02d" "$rep_num")"

            xtc="$(find "$field_dir" -maxdepth 1 -type f -name "*_${rep}.xtc" | head -n 1)"
            if [[ -z "$xtc" ]]; then
                echo "[skip] missing xtc for $charge/$field/rep_$rep"
                continue
            fi

            outdir="$field_dir/new_runs/rep_$rep"
            cfg="$outdir/config.json"
            logfile="$outdir/run.log"

            echo
            echo "charge: $charge  field: $field  rep: $rep"
            echo "xtc:  $xtc"
            echo "topo: $topo"
            echo "out:  $outdir"

            mkdir -p "$outdir"

            sed \
              -e "s|__XTC__|$xtc|g" \
              -e "s|__TOPOLOGY_JSON__|$topo|g" \
              -e "s|__OUTDIR__|$outdir|g" \
              "$BASE_CONFIG" > "$cfg"

            if ! "$SIMIO_BIN" --config "$cfg" > "$logfile" 2>&1; then
                echo "[error] failed: $charge/$field/rep_$rep"
                echo "[error] log tail: $logfile"
                tail -80 "$logfile" || true
                exit 1
            fi

            echo "[done] $charge/$field/rep_$rep"
        done
    done
done
