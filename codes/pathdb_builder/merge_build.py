#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_build.py
"""

import os, shutil, subprocess, pathlib
from typing import Tuple
from util import log, fasta_records_iter, sha1seq

def merge_group_fastas(group_dir: pathlib.Path, merged_path: pathlib.Path) -> Tuple[int, int]:
    """
    Merge all species fasta files under the category & remove duplicates by accesion.
    """
    seen, n_in, n_out = set(), 0, 0
    with open(merged_path, "w", encoding="utf-8") as out:
        for fa in sorted(group_dir.glob("*.fasta")):
            if fa.name.endswith(".part"): continue
            for hdr, seq in fasta_records_iter(str(fa)):
                n_in += 1
                acc = hdr[1:].strip().split()[0]
                if acc in seen: continue
                seen.add(acc); n_out += 1
                out.write(hdr + "\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i+80] + "\n")
    return n_in, n_out

def merge_group_taxmaps(group_dir: pathlib.Path, merged_taxmap: pathlib.Path) -> int:
    n = 0
    with open(merged_taxmap, "w", encoding="utf-8") as out:
        for tm in sorted(group_dir.glob("*.taxmap.tsv")):
            if tm.name.endswith(".part"): continue
            with open(tm, "r", encoding="utf-8") as f:
                for ln in f:
                    out.write(ln); n += 1
    return n
                    
def has_cd_hit_est() -> str:
    path = shutil.which("cd-hit-est")
    return path or ""

def dedup_by_cd_hit(in_fa: pathlib.Path, out_fa: pathlib.Path, threads: int=8, c: float=0.99) -> bool:
    exe = has_cd_hit_est()
    if not exe: return False
    cmd = [exe, "-i", str(in_fa), "-o", str(out_fa), "-c", 
           str(c), "-n", "10", "-T", str(threads), "-M", "0"]
    log("Eliminate redundancy:" + " ".join(cmd))
    code = subprocess.call(cmd)
    return code == 0

def dedup_by_hash(in_fa: pathlib.Path, out_fa: pathlib.Path) -> int:
    """Simple content deduplication."""
    seen, n_in, n_out = set(), 0, 0
    with open(out_fa, "w", encoding="utf-8") as out:
        for hdr, seq in fasta_records_iter(str(in_fa)):
            n_in += 1
            key = sha1seq(seq)
            if key in seen: continue
            seen.add(key); n_out += 1
            out.write(hdr + "\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")
    log(f"Eliminate redundancy(hash):input={n_in} output={n_out}")
    return n_out

def ensure_makeblastdb(bin_hint: str=None) -> str:
    cand = bin_hint or os.getenv("MAKEBLASTDB") or shutil.which("makeblastdb")
    if not cand:
        raise FileNotFoundError("makeblastdb not found.")
    return cand

def build_blast_db(makeblastdb: str, fasta_in: pathlib.Path, taxmap: pathlib.Path,
                   out_prefix: pathlib.Path, title: str):
    cmd = [makeblastdb, "-dbtype", "nucl",
           "-in", str(fasta_in), "-parse_seqids",
           "-title", title, "-out", str(out_prefix)]
    if taxmap and taxmap.exists():
        cmd += ["-taxid_map", str(taxmap)]
    log("Database creation: " + " ".join(cmd))
    code = subprocess.call(cmd)
    if code != 0:
        raise RuntimeError(f"makeblastdb fail. (Failing code:{code})")