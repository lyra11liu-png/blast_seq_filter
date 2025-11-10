export NCBI_API_KEY=1bd9501dcc57f92d38f6eb04e43c1bc21c09
export MAX_TOTAL=50000
python main.py \
  --csv /home/liuyuxin/projects/blast_seq_filter/data/name_resolution.csv \
  --taxid-col best_taxid --name-col best_name \
  --outdir /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --threads 4 --page-size 2000 --resume \
  --hybrid --virus-slen-min 200 --microbe-slen-min 1000 \
  --merge --dedup --build-db \
  --microbe-refseq-fallback 200




export NCBI_API_KEY=1bd9501dcc57f92d38f6eb04e43c1bc21c09
export MAX_TOTAL_VIRUSES=120000
export MAX_TOTAL_BACTERIA=60000
export MAX_TOTAL_FUNGI=60000
export MAX_TOTAL_ARCHAEA=40000

python main.py \
  --csv /home/liuyuxin/projects/blast_seq_filter/data/name_resolution.csv \
  --taxid-col best_taxid --name-col best_name \
  --outdir /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --dbdir  /data1/liuyuxin/blast_seq_filter/original_data_lib/db \
  --threads 12 --resume \
  --hybrid \
  --virus-slen-min 200 \
  --microbe-slen-min 1000 --microbe-refseq-fallback 200 \
  --tolerance 0.02 \
  --merge --dedup --build-db
