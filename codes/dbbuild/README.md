python -m dbbuild.postprocess \
  --indir  /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --dbdir  /data1/liuyuxin/blast_seq_filter/original_data_lib/db \
  --categories bacteria \
  --wrap 80 \
  --tmpdb /data1/liuyuxin/blast_seq_filter/original_data_lib/dedup_index.sqlite \
  --skip-count 0 \
  --update-every 5000

python -m dbbuild.postprocess \
  --indir  /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --dbdir  /data1/liuyuxin/blast_seq_filter/original_data_lib/db \
  --wrap 80 \
  --tmpdb :memory: \
  --skip-count 0 \
  --update-every 5000 \
  --commit-every 200000

MMSEQS_SWITCH_TO_BLAST_BYTES=$((1024*1024*1024)) \
python dedup_build.py \
  --indir  /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --outdir /data1/liuyuxin/blast_seq_filter/original_data_lib/clean \
  --dbdir  /data1/liuyuxin/blast_seq_filter/original_data_lib/db \
  --tmp-root /data1/liuyuxin/_fasttmp/dedup_tmp \
  --jobs 14 \
  --threads-per-task 12 \
  --pid 0.995 --qcov 0.98 \
  --resume
