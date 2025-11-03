export NCBI_API_KEY=1bd9501dcc57f92d38f6eb04e43c1bc21c09
python main.py \
  --csv /home/liuyuxin/projects/blast_seq_filter/data/name_resolution.csv \
  --taxid-col best_taxid --name-col best_name \
  --outdir /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --threads 10 \
  --page-size 5000 \
  --resume \
  --merge --dedup --build-db \
  --hybrid \
  --virus-slen-min 200 \
  --microbe-slen-min 1000