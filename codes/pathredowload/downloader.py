# -*- coding: utf-8 -*-
import sys, os
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List
from tqdm import tqdm
from .queries import build_query
from .ncbi import esearch_history, efetch_stream
from .io_utils import out_path, read_fasta_ids

def _category_slen_min(cat: str, cfg: Dict) -> int:
    if cat.lower() == "viruses":
        return int(cfg["virus_slen_min"])
    return int(cfg["microbe_slen_min"])

def _refseq_only(cat: str, cfg: Dict) -> bool:
    microbe_cats = {"bacteria","fungi","archaea","mycoplasma","protozoa","helminths"}
    return (cat.lower() in microbe_cats) and (not bool(cfg.get("microbes_include_genbank", 0)) is True)

def _streams_to_lines(chunks):
    """把 bytes chunk 转成逐行文本（utf-8，保留换行）"""
    buf = ""
    for ch in chunks:
        buf += ch.decode("utf-8", errors="ignore")
        while True:
            i = buf.find("\n")
            if i == -1: break
            line = buf[:i+1]
            buf = buf[i+1:]
            yield line
    if buf:
        if not buf.endswith("\n"): buf += "\n"
        yield buf

def _append_unique_from_stream(chunks, fo, existing_ids: set) -> Dict[str,int]:
    """
    将 EFETCH 的 fasta 文本流解析为记录，若 accession 未出现过则追加写入。
    返回 {"appended":x, "total":y}
    """
    total = appended = 0
    keep = False
    cur = []
    cur_id = ""
    for line in _streams_to_lines(chunks):
        if line.startswith(">"):
            if cur:
                total += 1
                if keep:
                    fo.write("".join(cur).encode("utf-8"))
                    appended += 1
                cur.clear()
            cur_id = line[1:].split()[0]
            keep = (cur_id not in existing_ids)
            if keep:
                existing_ids.add(cur_id)
            cur = [line]
        else:
            cur.append(line)
    if cur:
        total += 1
        if keep:
            fo.write("".join(cur).encode("utf-8"))
            appended += 1
    return {"appended": appended, "total": total}

def _one_taxid(row: Dict, cfg: Dict) -> Dict:
    taxid = str(row["best_taxid"])
    name  = str(row.get("best_name","NA"))
    cat   = str(row["category"])
    out   = out_path(cfg["outdir"], cat, taxid, name)

    # 已有文件策略
    if os.path.isfile(out) and os.path.getsize(out) > 0:
        if cfg["skip_existing"] and not cfg["supplement_existing"]:
            return {"taxid": taxid, "category": cat, "name": name, "status": "skip(exists)"}

    term = build_query(
        taxid, cat,
        slen_min=_category_slen_min(cat, cfg),
        only_complete=cfg["only_complete"],
        refseq_only=_refseq_only(cat, cfg),
        viruses_include_genbank=cfg["viruses_include_genbank"]
    )

    try:
        count, webenv, qk = esearch_history(term, cfg["api_key"], cfg["timeout"], cfg["retries"])
        if count == 0:
            return {"taxid": taxid, "category": cat, "name": name, "status": "no_records", "term": term}

        page = int(cfg["page_size"])
        written = 0
        if os.path.isfile(out) and os.path.getsize(out) > 0 and cfg["supplement_existing"]:
            # 补充模式：读取已有 accession 集合，增量追加
            existing = read_fasta_ids(out)
            with open(out, "ab") as fo:
                for start in range(0, count, page):
                    stat = _append_unique_from_stream(
                        efetch_stream(webenv, qk, cfg["api_key"], start, min(page, count-start),
                                      cfg["timeout"], cfg["retries"]),
                        fo, existing
                    )
                    written += stat["appended"]
            return {"taxid": taxid, "category": cat, "name": name, "status": f"supplement:+{written} of {count}"}
        else:
            # 新建模式：直接整流写入（不覆盖现有大文件）
            if os.path.exists(out) and not cfg["supplement_existing"]:
                # 有文件但选择了跳过
                return {"taxid": taxid, "category": cat, "name": name, "status": "skip(exists)"}
            with open(out, "wb") as fo:
                for start in range(0, count, page):
                    for chunk in efetch_stream(webenv, qk, cfg["api_key"], start, min(page, count-start),
                                               cfg["timeout"], cfg["retries"]):
                        fo.write(chunk)
                        written += len(chunk)
            return {"taxid": taxid, "category": cat, "name": name, "status": f"ok:{count}"}
    except Exception as e:
        return {"taxid": taxid, "category": cat, "name": name, "status": f"error:{e}"}

def download_all(df, cfg: Dict):
    rows: List[Dict] = df[["best_taxid","best_name","category"]].drop_duplicates().to_dict(orient="records")
    pbar = tqdm(total=len(rows), desc="Downloading TaxIDs", unit="taxid")
    stats = {"ok":0,"skip":0,"no_records":0,"error":0,"supplement":0}
    misses = []  # 记录未下载到的项（不建空文件）

    def _wrap(row):
        res = _one_taxid(row, cfg)
        pbar.update(1)
        s = res["status"]
        if s.startswith("ok"): stats["ok"] += 1
        elif s.startswith("skip"): stats["skip"] += 1
        elif s.startswith("no_records"):
            stats["no_records"] += 1
            misses.append(res)
        elif s.startswith("supplement"):
            stats["supplement"] += 1
        else:
            stats["error"] += 1
        return res

    with ThreadPoolExecutor(max_workers=int(cfg["threads"])) as ex:
        futs = [ex.submit(_wrap, r) for r in rows]
        for f in as_completed(futs):
            res = f.result()
            print(f"[{res['category']}] {res['taxid']} {res['name']} -> {res['status']}", flush=True)

    pbar.close()

    # 写出 not_found.tsv（仅一份简要清单）
    if misses:
        nf = os.path.join(cfg["outdir"], "not_found.tsv")
        with open(nf, "w", encoding="utf-8") as fo:
            fo.write("category\ttaxid\tname\tterm\n")
            for m in misses:
                fo.write(f"{m['category']}\t{m['taxid']}\t{m['name']}\t{m.get('term','')}\n")
        print(f"[INFO] Wrote not-found list -> {nf}", file=sys.stderr)

    print("Summary:", stats, file=sys.stderr)
