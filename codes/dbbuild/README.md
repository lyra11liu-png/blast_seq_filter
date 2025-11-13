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

export MMSEQS_SWITCH_TO_BLAST_BYTES=$((4*1024*1024*1024))
export BLAST_SUBJECT_DB=1
export BLAST_DB_MIN_SEQS=15000
export BLAST_DB_MIN_BYTES=$((400*1024*1024))
export BLAST_BATCH=512
export BLAST_DB_REFRESH=8000
export BLAST_WORD=48
export BLAST_UNGAPPED=1
export BLAST_SUBJECT_REFRESH=400
export BLAST_CLEAN_TMP=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1
python dedup_build.py \
  --indir  /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --outdir /data1/liuyuxin/blast_seq_filter/original_data_lib/clean \
  --dbdir  /data1/liuyuxin/blast_seq_filter/original_data_lib/db \
  --tmp-root /data1/liuyuxin/_fasttmp/dedup_tmp \
  --jobs 24 \
  --threads-per-task 6 \
  --pid 0.995 --qcov 0.98 \
  --resume
