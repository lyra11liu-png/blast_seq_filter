
export NCBI_API_KEY=1bd9501dcc57f92d38f6eb04e43c1bc21c09
python codes/constructdatalib/preflight_taxonomy_check.py \
  --manifest data/pathogens_CN_WHO_merged.csv \
  --out data/name_resolution.csv

