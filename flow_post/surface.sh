# from .../single_H_9/L_6/neutral
for fld in FIELD_00 FIELD_01 FIELD_02 FIELD_03; do
  # 1) runs/rep_XX/gating_flux.csv -> runs/rep_XX/transport_input.tsv
  python3 /data/antoine/Flow_CDI/flow_post/extract_transport_input.py "$fld"

  # 2) reduce from runs transport_input.tsv
  python3 /data/antoine/Flow_CDI/flow_post/reduce_counts.py "$fld/runs" --replica-glob "rep_*/transport_input.tsv"

  # 3) fit from runs transport_input.tsv
  python3 /data/antoine/Flow_CDI/flow_post/fit_transport.py "$fld/runs" --replica-glob "rep_*/transport_input.tsv" --fit-last-ns 45
done

python3 /data/antoine/Flow_CDI/flow_post/build_partial_master.py
python3 /data/antoine/Flow_CDI/flow_post/summarize_per_surface.py

# once you have at least one complete surface:
python3 /data/antoine/Flow_CDI/flow_post/ohmic_fit.py
