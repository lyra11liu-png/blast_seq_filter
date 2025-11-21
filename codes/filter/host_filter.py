#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
host_filter.py
Use minimap2 + samtools to remove human host seqs from sequencing .bam files.
Export the remaining read as .fasta files 4 subsequent blast identification.

Author: lyra Liu
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, Callable

from utils import run_cmd, Timer


def fastq_to_fasta(fq_path: str, fa_path: str,
                   logger: Optional[Callable[[str], None]] = None) -> None:
    """
    Convert fastq to fasta (line-by-line streaming processing).
    """
    if logger:
        logger(f"[INFO] Convert FASTQ to FASTA: {fq_path} -> {fa_path}")
    n = 0
    with open(fq_path, "r", encoding="utf-8") as fin, \
            open(fa_path, "w", encoding="utf-8") as fout:
        while True:
            header = fin.readline()
            if not header:
                break  # EOF
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                raise ValueError(f"Unexpected EOF in FASTQ file: {fq_path}")
            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ header line: {header.rstrip()}")
            header = header[1:].strip()
            seq = seq.strip()
            fout.write(f">{header}\n{seq}\n")
            n += 1
    if logger:
        logger(f"[INFO] FASTQ records converted: {n}")


def bam_to_nonhuman_fasta(
    bam_path: str,
    out_dir: str,
    human_index: str,
    minimap2: str = "minimap2",
    samtools: str = "samtools",
    threads: int = 8,
    preset: str = "map-hifi",
    logger: Optional[Callable[[str], None]] = None,
) -> Dict[str, float]:
    """
    Perform the following steps on a single .bam file:
        .bam -> .fastaq -> filter non-human seqs -> .fasta
    """
    bam = Path(bam_path)
    out_dir_p = Path(out_dir)
    out_dir_p.mkdir(parents=True, exist_ok=True)

    sample = bam.stem
    prefix = out_dir_p / sample

    nonhuman_fq = f"{prefix}.nonhuman.fastq"
    nonhuman_fa = f"{prefix}.nonhuman.fasta"

    timings: Dict[str, float] = {}

    # samtools fastq -> minimap2 align -> samtools view filter unmapped -> samtools fastq output
    pipe_cmd = (
        f"{samtools} fastq -n {bam_path} | "
        f"{minimap2} -t {threads} -ax {preset} {human_index} - | "
        f"{samtools} view -b -f 4 - | "
        f"{samtools} fastq -n - > {nonhuman_fq}"
    )
    with Timer("bam_to_nonhuman_fastq", logger) as t:
        run_cmd(pipe_cmd, logger=logger, shell=True)
    timings["bam_to_nonhuman_fastq"] = t.elapsed

    with Timer("fastq_to_fasta", logger) as t:
        fastq_to_fasta(nonhuman_fq, nonhuman_fa, logger=logger)
    timings["fastq_to_fasta"] = t.elapsed

    if logger:
        logger(f"[DONE] Host filtering finished for sample {sample}")
        logger(f"[OUTPUT] non-human FASTA: {nonhuman_fa}")

    return timings
