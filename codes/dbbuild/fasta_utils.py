# -*- coding: utf-8 -*-
import gzip, hashlib

def open_auto(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")

def iter_fasta_records(path):
    """yield (header, seq)；header不含'>'; seq为拼接后的大写序列"""
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

def write_fasta_record(fo, header, seq, wrap=80):
    fo.write(">" + header + "\n")
    if wrap and wrap > 0:
        for i in range(0, len(seq), wrap):
            fo.write(seq[i:i+wrap] + "\n")
    else:
        fo.write(seq + "\n")

def md5_blob(seq: str) -> bytes:
    """返回 16 字节的 MD5 BLOB，用于 SQLite BLOB 主键"""
    return hashlib.md5(seq.encode("utf-8")).digest()

def primary_id(header: str) -> str:
    return header.split()[0]

def unique_header(header: str, used_ids: set) -> str:
    pid = primary_id(header)
    if pid not in used_ids:
        used_ids.add(pid)
        return header
    n = 2
    base = header
    while True:
        cand = f"{base}|dup{n}"
        if primary_id(cand) not in used_ids:
            used_ids.add(primary_id(cand))
            return cand
        n += 1
