export NCBI_API_KEY=1bd9501dcc57f92d38f6eb04e43c1bc21c09
python -m pathredowload.main \
  --csv /data1/liuyuxin/blast_seq_filter/original_data_lib/redow_pathogen_resolution.csv \
  --outdir /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --threads 6 \
  --qps 8 \
  --virus-slen-min 0 \
  --microbe-slen-min 0 \
  --viruses-include-genbank 1 \
  --microbes-include-genbank 1 \
  --refseq-min-records 1 \
  --only-complete 0 \
  --skip-existing 1 \
  --supplement-existing 1