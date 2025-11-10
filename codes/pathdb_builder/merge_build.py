#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_build.py
"""

import os, shutil, subprocess, pathlib
from typing import Tuple
from util import log, fasta_records_iter, sha1seq

def merge_group_fastas(group_dir: pathlib.Path, merged_path: pathlib.Path, progress_every: int = 200000) -> Tuple[int, int]:
    """
    流式合并组内所有 FASTA 到 merged_path。
    不做“跨全组 accession 去重”（nuccore accession 全库唯一；大 set 会吃爆内存）。
    支持 .fa/.fasta/.fa.gz/.fasta.gz；按 80 列换行输出；定期打印进度。
    """
    import gzip
    def _iter_fa(p: pathlib.Path):
        gz = (p.suffix == ".gz") or str(p).endswith(".fa.gz") or str(p).endswith(".fasta.gz")
        opener = gzip.open if gz else open
        mode   = "rt" if gz else "r"
        with opener(p, mode, encoding="utf-8", errors="ignore") as f:
            hdr, seq = None, []
            for ln in f:
                if ln.startswith(">"):
                    if hdr is not None:
                        yield hdr, "".join(seq)
                    hdr, seq = ln.rstrip("\n"), []
                else:
                    seq.append(ln.strip())
            if hdr is not None:
                yield hdr, "".join(seq)

    files = sorted(list(group_dir.glob("*.fasta")) +
                   list(group_dir.glob("*.fa")) +
                   list(group_dir.glob("*.fasta.gz")) +
                   list(group_dir.glob("*.fa.gz")))
    files = [p for p in files if (not p.name.endswith(".part")) and (not p.exists() or p.stat().st_size > 0)]
    
    n_in, n_out = 0, 0
    merged_path.parent.mkdir(parents=True, exist_ok=True)
    with open(merged_path, "w", encoding="utf-8") as out:
        for fi, fa in enumerate(files, 1):
            if fa.name.endswith(".part"):
                continue
            if fa.exists() and fa.stat().st_size == 0:
                continue
            for hdr, seq in _iter_fa(fa):
                n_in += 1; n_out += 1
                out.write(hdr + "\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i+80] + "\n")
                if (n_out % progress_every) == 0:
                    log(f"[merge] {group_dir.name}: {n_out} seq written ... (file {fi}/{len(files)})")
    log(f"[merge] {group_dir.name} done: input={n_in}, output={n_out}, to={merged_path}")
    return n_in, n_out

def merge_group_taxmaps(group_dir: pathlib.Path, merged_taxmap: pathlib.Path, progress_every: int = 500000) -> int:
    """
    由组内 FASTA 直接生成 taxid_map（seqid \t taxid），避免“拼接 .taxmap.tsv”带来的重复与暴涨。
    依赖文件名形如：{taxid}_{Name}.fasta(.gz)
    """
    import gzip, re
    def _iter_hdrs(p: pathlib.Path):
        gz = (p.suffix == ".gz") or str(p).endswith(".fa.gz") or str(p).endswith(".fasta.gz")
        opener = gzip.open if gz else open
        mode   = "rt" if gz else "r"
        with opener(p, mode, encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if ln.startswith(">"):
                    yield ln[1:].strip().split()[0]  # accession/seqid

    files = sorted(list(group_dir.glob("*.fasta")) +
                   list(group_dir.glob("*.fa")) +
                   list(group_dir.glob("*.fasta.gz")) +
                   list(group_dir.glob("*.fa.gz")))

    merged_taxmap.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(merged_taxmap, "w", encoding="utf-8") as out:
        for fa in files:
            if fa.name.endswith(".part"):
                continue
            m = re.match(r"(\d+)_", fa.name)
            if not m:
                log(f"[taxmap] skip (no taxid prefix): {fa.name}")
                continue
            taxid = m.group(1)
            for seqid in _iter_hdrs(fa):
                out.write(f"{seqid}\t{taxid}\n")
                n += 1
                if (n % progress_every) == 0:
                    log(f"[taxmap] {group_dir.name}: {n} lines ...")
    log(f"[taxmap] {group_dir.name} done: {n} lines -> {merged_taxmap}")
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

def filter_taxmap_by_fasta(fa_path: pathlib.Path, taxmap_in: pathlib.Path, taxmap_out: pathlib.Path) -> int:
    """仅保留 taxmap 中与 fa_path 内 seqid 匹配的行。"""
    keep = set()
    for hdr, _ in fasta_records_iter(str(fa_path)):
        keep.add(hdr[1:].strip().split()[0])
    n = 0
    with open(taxmap_in, "r", encoding="utf-8") as fin, open(taxmap_out, "w", encoding="utf-8") as fout:
        for ln in fin:
            sid = ln.split("\t", 1)[0]
            if sid in keep:
                fout.write(ln); n += 1
    log(f"[taxmap] filtered by {fa_path.name}: kept {n} lines")
    return n
