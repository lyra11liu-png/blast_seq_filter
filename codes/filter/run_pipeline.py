#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_pipeline.py

Main script: execute the following on all .bam files within a directory:
  - Remove human seqs using minimap2 + samtools
  - Perform blastn on remaining reads against multiple pathogen blast databases
  - Output all hits per sample & species statistics table

Support optional multi-process parallelism & basic computation time monitoring.

Author: lyra Liu
"""

from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from utils import default_logger_factory
from host_filter import bam_to_nonhuman_fasta
from blast_classify import run_multi_db_blast


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Host-depletion (human) + multi-DB pathogen BLAST classification pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--indir", required=True,
                   help="Enter the directory containing the .bam files.")
    p.add_argument("--bam", default=None,
                   help="Process only the specified single .bam file, providing the full path.")
    p.add_argument("--outdir", required=True,
                   help="Output directory structure: one subdirectory per sample.")

    p.add_argument("--human-index", required=True,
                   help="minimap2 human reference index .mmi path.")

    p.add_argument("--db-root", required=True,
                   help="Pathogen blast database directory.")
    p.add_argument("--db-groups", default="viruses,bacteria,fungi,mycoplasma,archaea,protozoa,helminths",
                   help=("The prefix names of the libraries to be used."))
    p.add_argument("--idmap-root", default=None,
                   help=("Directory containing the <group>.idmap.tsv file generated during categorized database creation;"
                         "If not specified, the default location will be the clean/subdirectory within the parent directory of db-root."))

    p.add_argument("--minimap2", default="minimap2",
                   help="minimap2 executable file name/path.")
    p.add_argument("--samtools", default="samtools",
                   help="samtools executable file name/path.")
    p.add_argument("--blastn", default="blastn",
                   help="blastn executable file name/path.")

    p.add_argument("--threads", type=int, default=8,
                   help="Numeber of threads used internally per sample.")
    p.add_argument("--preset", default="map-hifi",
                   help="minimap2 presets.")

    p.add_argument("--parallel", action="store_true",
                   help="Whether multiple .bam samples process in paralle.")
    p.add_argument("--jobs", type=int, default=2,
                   help="Number of processes in parallel.")

    p.add_argument("--min-pident", type=float, default=85.0,
                   help="Mimimum matching similarity.")
    p.add_argument("--min-align-len", type=int, default=80,
                   help="Minimum match length.")
    p.add_argument("--max-hit-evalue", type=float, default=1e-20,
                   help="Maximum hit e-value.")

    return p.parse_args()


def find_bam_files(indir: str) -> List[Path]:
    d = Path(indir)
    if not d.is_dir():
        raise SystemExit(f"Input directory not found: {indir}")
    bams = sorted(d.glob("*.bam"))
    if not bams:
        raise SystemExit(f"No BAM files found in: {indir}")
    return bams


def build_db_map(db_root: str, groups_str: str) -> Dict[str, str]:
    root = Path(db_root)
    if not root.is_dir():
        raise SystemExit(f"DB root directory not found: {db_root}")

    groups = [g.strip() for g in groups_str.split(",") if g.strip()]
    db_map: Dict[str, str] = {}

    for g in groups:
        prefix = root / g
        nhr = prefix.with_suffix(".nhr")
        nsq = prefix.with_suffix(".nsq")
        nal = prefix.with_suffix(".nal")
        nsq_parts = list(root.glob(f"{g}*.nsq"))

        if nhr.exists() or nsq.exists() or nal.exists() or nsq_parts:
            db_map[g] = str(prefix)
        else:
            print(f"[WARN] BLAST DB for group '{g}' not found under {db_root}, skip.")

    if not db_map:
        raise SystemExit("No valid BLAST DB found. Please check --db-root and --db-groups.")

    print("[INFO] Using BLAST DBs:")
    for g, p in db_map.items():
        print(f"  {g}: {p}")
    return db_map


def process_one_sample(
    bam_path: str,
    outdir: str,
    human_index: str,
    db_map: Dict[str, str],
    minimap2: str,
    samtools: str,
    blastn: str,
    threads: int,
    preset: str,
    idmap_root: Optional[str],
    min_pident: float,
    min_align_len: int,
    max_hit_evalue: float,
) -> Tuple[str, Dict[str, float]]:
    """
    Single sample processing function for the main process / subprocess.
    """
    bam = Path(bam_path)
    sample = bam.stem

    sample_outdir = Path(outdir) / sample
    sample_outdir.mkdir(parents=True, exist_ok=True)

    log_path = sample_outdir / f"{sample}.log"
    logger = default_logger_factory(str(log_path))

    logger(f"[SAMPLE] {sample}")
    logger(f"[INPUT BAM] {bam_path}")

    # Remove human reads & generate non-human .fasta files
    host_timings = bam_to_nonhuman_fasta(
        bam_path=str(bam),
        out_dir=str(sample_outdir),
        human_index=human_index,
        minimap2=minimap2,
        samtools=samtools,
        threads=threads,
        preset=preset,
        logger=logger,
    )
    nonhuman_fa = str(sample_outdir / f"{sample}.nonhuman.fasta")

    # Perform blastn searches on non-human fasta seqs across multiple databases.
    blast_timings = run_multi_db_blast(
        query_fa=nonhuman_fa,
        out_prefix=str(sample_outdir / sample),
        db_map=db_map,
        blastn=blastn,
        threads=threads,
        idmap_root=idmap_root,
        logger=logger,
        min_pident=min_pident,
        min_align_len=min_align_len,
        max_hit_evalue=max_hit_evalue,
    )

    # Summary of timings
    timings: Dict[str, float] = {}
    timings.update(host_timings)
    timings.update(blast_timings)
    return sample, timings


def write_runtime_summary(all_timings: Dict[str, Dict[str, float]], outdir: str) -> None:
    """
    Write the timing information 4 all samples into a single .tsv summary file 4 monitoring purposes.
    """
    out_path = Path(outdir) / "runtime_summary.tsv"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("sample\tstep\tseconds\n")
        for sample, tdict in sorted(all_timings.items()):
            for step, sec in sorted(tdict.items()):
                f.write(f"{sample}\t{step}\t{sec:.2f}\n")
    print(f"[RUNTIME] Summary written to: {out_path}")


def merge_taxonomy_reports(outdir: str) -> None:
    out_dir = Path(outdir)
    merged_path = out_dir / "taxonomy_summary.tsv"

    with open(merged_path, "w", encoding="utf-8") as fout:
        fout.write("Sample\tTaxonomy\n")
        for sub in sorted(out_dir.iterdir()):
            if not sub.is_dir():
                continue
            sample = sub.name
            tax_fp = sub / f"{sample}.taxonomy_report.tsv"
            if not tax_fp.is_file():
                continue
            with open(tax_fp, "r", encoding="utf-8") as fin:
                header = fin.readline()
                for line in fin:
                    fout.write(line)

    print(f"[TAXONOMY] Merged taxonomy table written to: {merged_path}")


def main() -> None:
    args = parse_args()
    
    if args.bam is not None:
        bam_path = Path(args.bam)
        if not bam_path.is_file():
            raise SystemExit(f"--bam file not found: {args.bam}")
        bams = [bam_path]
        print(f"[INFO] Only process single BAM file: {args.bam}")
    else:
        bams = find_bam_files(args.indir)
        print(f"Found {len(bams)} BAM files in {args.indir}")
        
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # Construct the group -> db_prefix mapping.
    db_map = build_db_map(args.db_root, args.db_groups)
    
    if args.idmap_root is not None:
        idmap_root = args.idmap_root
    else:
        idmap_root = str((Path(args.db_root).resolve().parent / "clean"))
    print(f"[INFO] Using idmap_root: {idmap_root}")

    all_timings: Dict[str, Dict[str, float]] = {}

    if args.parallel:
        jobs = max(1, int(args.jobs))
        print(f"Running in parallel with {jobs} processes")
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            future_to_bam = {
                ex.submit(
                    process_one_sample,
                    str(bam),
                    args.outdir,
                    args.human_index,
                    db_map,
                    args.minimap2,
                    args.samtools,
                    args.blastn,
                    args.threads,
                    args.preset,
                    idmap_root,
                    args.min_pident,
                    args.min_align_len,
                    args.max_hit_evalue,
                ): bam
                for bam in bams
            }
            for fut in as_completed(future_to_bam):
                bam = future_to_bam[fut]
                try:
                    sample, timings = fut.result()
                    all_timings[sample] = timings
                    print(f"[DONE] sample {sample}")
                except Exception as e:
                    print(f"[ERROR] sample {bam.name} failed: {e}")
    else:
        print("Running samples sequentially (no --parallel)")
        for bam in bams:
            try:
                sample, timings = process_one_sample(
                    str(bam),
                    args.outdir,
                    args.human_index,
                    db_map,
                    args.minimap2,
                    args.samtools,
                    args.blastn,
                    args.threads,
                    args.preset,
                    idmap_root,
                    args.min_pident,
                    args.min_align_len,
                    args.max_hit_evalue,
                )
                all_timings[sample] = timings
                print(f"[DONE] sample {sample}")
            except Exception as e:
                print(f"[ERROR] sample {bam.name} failed: {e}")

    # Create a comprehensive timeline & merge the tax table.
    write_runtime_summary(all_timings, args.outdir)
    merge_taxonomy_reports(args.outdir)


if __name__ == "__main__":
    main()
