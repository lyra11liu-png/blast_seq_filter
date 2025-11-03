#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
util.py
General utility function set:
- Unified log output
- Directory creation & filename sanitization
- Count fasta lines
- Safely count complete records in .fasta.part files
- Iterator
- Simple seqs content hashing
"""

import os, re, hashlib, pathlib
from datetime import datetime

def now(): 
    """Return the current time as strings."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
def log(msg): 
    """Output log."""
    print(f"[{now()}] {msg}", flush=True)

def ensure_dir(p: pathlib.Path):
    """Ensure path exists."""
    p.mkdir(parents=True, exist_ok=True); return p
    
def sanitize(s: str) -> str:
    """Clean strings for safe filenames."""
    s = (s or "").strip()
    s = re.sub(r"\s+", "_", s)
    return re.sub(r"[^A-Za-z0-9_.-]", "", s)[:80] or "unnamed"

def count_lines(path: str) -> int:
    """Count the number of lines in a text file."""
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            return sum(1 for _ in f)
    except FileNotFoundError:
        return 0
    
def count_fasta_records(path: str) -> int:
    """Count the number of > in .fasta."""
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            return sum(1 for ln in f if ln.startswith(">"))
    except FileNotFoundError:
        return 0
    
def fasta_records_iter(path: str):
    """Iterate througn fasta records."""
    hdr, seq = None, []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq)
                hdr, seq = ln.rstrip("\n"), []
            else:
                seq.append(ln.strip())
        if hdr is not None:
            yield hdr, "".join(seq)
            
def count_fasta_records_safe(part_path: str) -> int:
    """Count .fasta.part numbers."""
    if not os.path.exists(part_path): return 0
    ends_with_nl = True
    try:
        with open(part_path, "rb") as fb:
            fb.seek(-1, os.SEEK_END)
            ends_with_nl = (fb.read(1) == b"\n")
    except Exception:
        ends_with_nl = True
        
    count, last_hdr_off, off = 0, None, 0
    with open(part_path, "rb") as fb:
        for ln in fb:
            if ln.startswith(b">"):
                last_hdr_off = off
                count += 1
            off += len(ln)
            
    if (not ends_with_nl) and (last_hdr_off is not None):
        with open(part_path, "rb+") as fw:
            fw.truncate(last_hdr_off)
        count = max(0, count - 1)
    return count

def sha1seq(s: str) -> str:
    return hashlib.sha1(s.encode("ascii","ignore")).hexdigest()