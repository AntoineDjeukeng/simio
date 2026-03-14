# from .../single_H_9/L_6/neutral
for fld in FIELD_00 FIELD_01 FIELD_02 FIELD_03; do
  python3 /data/antoine/Flow_CDI/flow_post/reduce_counts.py "$fld/compile"
  python3 /data/antoine/Flow_CDI/flow_post/fit_transport.py "$fld/compile" --fit-last-ns 45
done

python3 /data/antoine/Flow_CDI/flow_post/build_partial_master.py
python3 /data/antoine/Flow_CDI/flow_post/summarize_per_surface.py

# once you have at least one complete surface:
python3 /data/antoine/Flow_CDI/flow_post/ohmic_fit.py
