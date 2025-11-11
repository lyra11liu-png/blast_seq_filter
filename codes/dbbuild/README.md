python -m dbbuild.postprocess \
  --indir  /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --dbdir  /data1/liuyuxin/blast_seq_filter/original_data_lib/db \
  --wrap 80 \
  --skip-count 0 \
  --tmpdb /data1/liuyuxin/blast_seq_filter/original_data_lib/dedup_index.sqlite \
  --update-every 2000 \
  --keep-merged 0