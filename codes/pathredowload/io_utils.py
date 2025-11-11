# -*- coding: utf-8 -*-
import os, re

def ensure_dir(d: str):
    os.makedirs(d, exist_ok=True)

def sanitize_name(s: str) -> str:
    s = (s or "NA")
    s = re.sub(r"[^\w\-\.\+]+", "_", s.strip())
    s = re.sub(r"_+", "_", s)
    return s.strip("_")

def out_path(outdir: str, category: str, taxid: str, name: str) -> str:
    fname = f"{taxid}_{sanitize_name(name)}.fasta"
    d = os.path.join(outdir, category.lower())
    ensure_dir(d)
    return os.path.join(d, fname)

def read_fasta_ids(path: str) -> set:
    """读取已有 fasta 的 accession（取 '>' 后首个空格前 token）"""
    ids = set()
    if not os.path.isfile(path): return ids
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fi:
            for ln in fi:
                if ln.startswith(">"):
                    tok = ln[1:].split()[0]
                    if tok: ids.add(tok)
    except Exception:
        pass
    return ids
