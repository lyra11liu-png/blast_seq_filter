python run_pipeline.py \
  --indir /data1/liuyuxin/blast_seq_filter/w_filter_data/unfiltereddata/Rawdata \
  --outdir /data1/liuyuxin/blast_seq_filter/results \
  --human-index /data1/liuyuxin/blast_seq_filter/original_data_lib/db/human/human.mmi \
  --db-root /data1/liuyuxin/blast_seq_filter/original_data_lib/db \
  --db-groups viruses,bacteria,fungi,mycoplasma \
  --threads 16 \
  --preset map-hifi \
  --min-pident 90 \
  --min-align-len 150 \
  --max-hit-evalue 1e-30 \
  --parallel --jobs 10

python run_pipeline.py \
  --indir /data1/liuyuxin/blast_seq_filter/w_filter_data/unfiltereddata/Rawdata \
  --bam   /data1/liuyuxin/blast_seq_filter/w_filter_data/unfiltereddata/Rawdata/ATCC25922_LHG23935.hifi_reads.bam \
  --outdir /data1/liuyuxin/blast_seq_filter/host_removed_blast_results \
  --human-index /data1/liuyuxin/blast_seq_filter/original_data_lib/db/human/human.mmi \
  --db-root /data1/liuyuxin/blast_seq_filter/original_data_lib/db \
  --db-groups viruses,bacteria,fungi,mycoplasma \
  --threads 16 \
  --preset map-hifi \
  --min-pident 90 \
  --min-align-len 150 \
  --max-hit-evalue 1e-30