# Host-depletion and multi-database BLAST classification pipeline

**Author:** lyra Liu  

A small Python pipeline for:

1. Removing human host reads from BAM files using **minimap2 + samtools**.
2. Classifying the remaining reads against multiple pathogen BLAST databases (viruses, bacteria, fungi, etc.).
3. Producing per-sample and overall summary tables that are easy to open in Excel.

The pipeline is designed for clinical sequencing data where host reads need to be removed before pathogen classification.

---

## Features

- **Host depletion from BAM**  
  Uses `samtools fastq` → `minimap2` mapping to a human reference → `samtools view -f4` to keep only unmapped (non-human) reads, then converts them to FASTA.

- **Multi-database pathogen BLAST**  
  Runs `blastn` against one or more databases (e.g. `viruses`, `bacteria`, `fungi`, …) and merges results.

- **Hit filtering by alignment quality**  
  Only keeps hits that pass user-defined thresholds:  
  - minimum identity (`--min-pident`)  
  - minimum alignment length (`--min-align-len`)  
  - maximum e-value (`--max-hit-evalue`)

- **Clean species names and grouped reports**  
  Uses optional `*.idmap.tsv` files to normalize species names and group hits by major category.

- **Per-sample logs and timing**  
  Each sample gets its own log file and timing breakdown; a global `runtime_summary.tsv` is also generated.

- **Optional parallel processing**  
  Process multiple BAM files in parallel with `--parallel --jobs N`.

---

## Repository structure

Key scripts:

- `utils.py`  
  - `run_cmd(...)`: small wrapper around `subprocess.run` with better error handling and optional logging.  
  - `Timer` context manager: measure wall-clock time for each step and report via logger.  
  - `default_logger_factory(...)`: create a simple logger that writes both to stdout and to a log file.

- `host_filter.py`  
  - `fastq_to_fasta(...)`: stream conversion from FASTQ to FASTA (low memory).  
  - `bam_to_nonhuman_fasta(...)`:  
    - `samtools fastq` converts BAM → FASTQ  
    - `minimap2` maps to human reference  
    - `samtools view -f4` keeps only unmapped reads  
    - `samtools fastq` exports non-human reads to FASTQ  
    - `fastq_to_fasta` converts non-human FASTQ → FASTA  

- `blast_classify.py`  
  - `run_multi_db_blast(...)`:
    - Iterates over multiple BLAST databases (groups like `viruses`, `bacteria`, `fungi`, etc.).  
    - Runs `blastn` for each group.  
    - Applies hit filters: `min_pident`, `min_align_len`, `max_hit_evalue`.  
    - Resolves taxonomy from optional `group.idmap.tsv` files.  
    - Writes:
      - `<prefix>.<group>.blast.tsv` – raw BLAST tabular output per group  
      - `<prefix>.all_hits.tsv` – all filtered hits merged  
      - `<prefix>.summary_by_species.tsv` – one row per `(group, kingdom, species)`  
      - `<prefix>.taxonomy_report.tsv` – Kraken-like text summary per sample

- `run_pipeline.py`  
  Main entry point. For one or many BAM files, it:
  1. Runs host filtering (`bam_to_nonhuman_fasta`) to generate `<sample>.nonhuman.fasta`.  
  2. Runs multi-DB BLAST (`run_multi_db_blast`) on the non-human FASTA.  
  3. Writes:
     - One subdirectory per sample with logs + BLAST outputs + summaries  
     - `runtime_summary.tsv` – per-sample timing table  
     - `taxonomy_summary.tsv` – merged taxonomy report for all samples

---

## Dependencies

### Software

- Python ≥ 3.9
- [minimap2](https://github.com/lh3/minimap2) (for host mapping)
- [samtools](http://www.htslib.org/) (for BAM and FASTQ processing)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (`blastn`, `makeblastdb`)

> All Python imports are from the standard library (`argparse`, `pathlib`, `subprocess`, `concurrent.futures`, etc.), so no extra Python packages are required.

### Reference data

- A minimap2 human index (`.mmi`) built from a human reference genome, e.g.:

  ```bash
  minimap2 -d human.mmi human_reference.fa

- One or more pathogen BLAST databases under a common directory, for example:

  ```text
  /path/to/db_root/
    ├── viruses.nsq / .nhr / .nin / ...
    ├── bacteria.nsq / .nhr / .nin / ...
    ├── fungi.nsq / .nhr / .nin / ...
    ├── mycoplasma.*
    ├── archaea.*
    ├── protozoa.*
    └── helminths.*


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