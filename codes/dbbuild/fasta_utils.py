# -*- coding: utf-8 -*-
import gzip, hashlib

def open_auto(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")

def iter_fasta_records(path):
    """yield (header, seq)，header不含'>'; seq为拼接后的大写序列（去空白）"""
    with open_auto(path) as f:
        h, seq = None, []
        for ln in f:
            if ln.startswith(">"):
                if h is not None:
                    yield h, "".join(seq).upper()
                h = ln.strip()[1:]
                seq = []
            else:
                seq.append(ln.strip())
        if h is not None:
            yield h, "".join(seq).upper()

def count_fasta_records(path) -> int:
    """快速计数：统计'>'行数（兼容.gz）"""
    n = 0
    with open_auto(path) as f:
        for ln in f:
            if ln.startswith(">"):
                n += 1
    return n

def write_fasta_record(fo, header, seq, wrap=80):
    fo.write(">" + header + "\n")
    if wrap and wrap > 0:
        for i in range(0, len(seq), wrap):
            fo.write(seq[i:i+wrap] + "\n")
    else:
        fo.write(seq + "\n")

def seq_md5(seq: str) -> str:
    # 已在 iter 中转大写，这里再保底一次
    return hashlib.md5(seq.upper().encode("utf-8")).hexdigest()

def primary_id(header: str) -> str:
    return header.split()[0]

def unique_header(header: str, used_ids: set) -> str:
    """确保 -parse_seqids 唯一：若重复则追加|dupN"""
    pid = primary_id(header)
    if pid not in used_ids:
        used_ids.add(pid)
        return header
    n = 2
    base = header
    while True:
        cand = f"{base}|dup{n}"
        pid2 = primary_id(cand)
        if pid2 not in used_ids:
            used_ids.add(pid2)
            return cand
        n += 1
