#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main.py
Command-line entry:
  - Read .csv
  - Parallelize downloads 4 each taxid
  - Write statistics tables
"""

import csv, re, time, pathlib, argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from util import log, ensure_dir
from downloader import download_one
from merge_build import (merge_group_fastas, merge_group_taxmaps,
                         dedup_by_cd_hit, dedup_by_hash, ensure_makeblastdb, build_blast_db)

GROUPS = ["viruses", "bacteria", "fungi", "archaea", "other"]

def parse_taxid(raw: str):
    if not raw: return []
    return [p for p in re.split(r"[|,;\s]+", raw) if p.isdigit()]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="包含 TaxID 的 CSV（如 best_taxid 列可含多 ID）")
    ap.add_argument("--taxid-col", default="best_taxid", help="CSV 中 TaxID 列名")
    ap.add_argument("--name-col", default="best_name", help="CSV 中名称列名（用于命名）")
    ap.add_argument("--outdir", required=True, help="输出根目录")
    ap.add_argument("--threads", type=int, default=4, help="并行下载的线程数")

    ap.add_argument("--no-subtree", action="store_true", help="不使用 txid[Subtree]（默认使用）")
    ap.add_argument("--no-genomic-only", action="store_true", help="不限制 biomol_genomic（默认限制）")
    ap.add_argument("--no-refseq-only", action="store_true", help="不限制 RefSeq（默认限制）")
    ap.add_argument("--slen-min", type=int, default=1000, help="序列长度下限（默认 1000）")
    ap.add_argument("--slen-max", type=int, default=100000000, help="序列长度上限（默认 1e8）")
    ap.add_argument("--page-size", type=int, default=5000, help="efetch/esummary 每页数量")
    ap.add_argument("--hybrid", action="store_true",
                    help="按门类自动策略：Viruses=放开RefSeq+低阈值；微生物=RefSeq-only+高阈值")
    ap.add_argument("--virus-slen-min", type=int, default=200,
                    help="hybrid模式下：病毒序列长度下限（默认200）")
    ap.add_argument("--microbe-slen-min", type=int, default=1000,
                    help="hybrid模式下：细菌/真菌/古菌序列长度下限（默认1000）")
    ap.add_argument("--microbe-refseq-fallback", type=int, default=200,
                help="微生物物种在 RefSeq-only 下命中条目 <该阈值时，自动放宽到 GenBank+RefSeq（默认200）")

    ap.add_argument("--resume", action="store_true", help="断点续传并跳过已完整文件")
    ap.add_argument("--tolerance", type=float, default=0.02, help="判定完整的容忍度（默认 2%）")
    ap.add_argument("--no-verify-skip", action="store_true", help="不做完整性校验；已有文件也重下（谨慎）")

    # Merge & database creation.
    ap.add_argument("--merge", action="store_true", help="按门类合并（并按 accession 去重）")
    ap.add_argument("--dedup", action="store_true", help="合并后再做序列去冗余（优先 cd-hit-est 0.99）")
    ap.add_argument("--build-db", action="store_true", help="用 makeblastdb 生成核酸库")
    ap.add_argument("--makeblastdb", default=None, help="makeblastdb 路径（不填将从环境中查找）")
    ap.add_argument("--dbdir", default=None,
                help="可选：把 makeblastdb 的输出前缀放到这个目录（默认放各 group 子目录）")
    
    args = ap.parse_args()
    out_base = ensure_dir(pathlib.Path(args.outdir).resolve())
    db_dir = pathlib.Path(args.dbdir).resolve() if args.dbdir else None
    if db_dir: ensure_dir(db_dir)
    
    # ===== Collect tasks =====
    tasks = []
    with open(args.csv, newline="", encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for row in rd:
            label = (row.get(args.name_col) or row.get("latin_input") or "").strip()
            raw = (row.get(args.taxid_col) or "").strip()
            for tx in parse_taxid(raw):
                tasks.append((tx, label))
    log(f"Totally {len(tasks)} tasks.")
    
    # ===== Parallelize downloading =====
    results = []
    if tasks:
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = [ex.submit(
                download_one, tx, label, out_base,
                subtree=not args.no_subtree,
                genomic_only=not args.no_genomic_only,
                refseq_only=not args.no_refseq_only,
                slen_min=args.slen_min, slen_max=args.slen_max,
                page=args.page_size,
                resume=args.resume,
                verify_skip=not args.no_verify_skip,
                tolerance=args.tolerance,
                hybrid=args.hybrid,
                virus_slen_min=args.virus_slen_min,
                microbe_slen_min=args.microbe_slen_min,
                microbe_refseq_fallback=args.microbe_refseq_fallback
            ) for tx,label in tasks]
            for fu in as_completed(futs):
                try:
                    results.append(fu.result())
                except Exception as e:
                    log(f"[Tasks Fail:] {e}")

        log(f"Download phase completed, time taken: {time.time()-t0:.1f}s.")
        
    # ===== Stastics tables =====
    stat = out_base / "download_stats.tsv"
    with open(stat, "w", encoding="utf-8") as f:
        f.write("taxid\tgroup\tsci_name\trank\tfile\ttaxmap\tcount\texpected\tskipped\tseconds\tincomplete\n")
        for r in sorted(results, key=lambda x: (x["group"], int(x["taxid"]))):
            f.write("\t".join([
                str(r["taxid"]), r["group"], r.get("sci_name",""), r.get("rank",""),
                r["file"], r["taxmap"], str(r["count"]), str(r.get("expected",0)),
                str(r.get("skipped",False)), f'{r["seconds"]:.1f}', str(r.get("incomplete",False))
            ]) + "\n")
    log(f"Wrote stastics:{stat}")
    
    miss = out_base / "incomplete_or_missing.tsv"
    with open(miss, "w", encoding="utf-8") as f:
        f.write("taxid\tgroup\tsci_name\trank\tfile\texpected\tcount\n")
        for r in sorted(results, key=lambda x: (x["group"], int(x["taxid"]))):
            expd = int(r.get("expected",0)); cnt = int(r.get("count",0) or 0)
            inc  = bool(r.get("incomplete", False))
            if inc or (expd>0 and cnt==0 and not r.get("skipped",False)):
                f.write("\t".join([
                    str(r["taxid"]), r["group"], r.get("sci_name",""), r.get("rank",""),
                    r.get("file",""), str(expd), str(cnt)
                ]) + "\n")
    log(f"Loss list:{miss}")
    
    # ===== Merge & datebase creation =====
    if args.merge or args.build_db or args.dedup:
        mk = ensure_makeblastdb(args.makeblastdb) if args.build_db else None
        for g in GROUPS:
            gdir = out_base / g
            if not gdir.exists():
                continue

            merged = gdir / f"{g}_all.fasta"
            merged_tax = gdir / f"{g}_all.taxmap.tsv"
            log(f"[Merge] {g} -> {merged.name} / {merged_tax.name}")

            # 1) 合并
            n_in, n_out = merge_group_fastas(gdir, merged)
            n_map = merge_group_taxmaps(gdir, merged_tax)
            log(f"  Merge finished: input={n_in}, output={n_out}; taxmap={n_map}")

            # 2) 去重（cd-hit-est 优先，缺失则回退哈希去重）
            dedup_out = merged
            tax_for_db = merged_tax
            if args.dedup:
                dedup_out = gdir / f"{g}_all.dedup.fasta"
                if not dedup_by_cd_hit(merged, dedup_out, threads=max(1, args.threads), c=0.99):
                    log("cd-hit-est not found, fallback to hash-dedup")
                    _ = dedup_by_hash(merged, dedup_out)

                # taxmap 只保留 dedup 后仍存在的 seqid
                from merge_build import filter_taxmap_by_fasta
                filtered_tax = gdir / f"{g}_all.dedup.taxmap.tsv"
                filter_taxmap_by_fasta(dedup_out, merged_tax, filtered_tax)
                tax_for_db = filtered_tax

            # 3) 建库
            if args.build_db:
                if db_dir:
                    prefix = db_dir / f"{g}_blastdb"
                else:
                    prefix = gdir / f"{g}_blastdb"
                title  = f"{g}_nucl"
                build_blast_db(mk, dedup_out, tax_for_db, prefix, title)

    log("ALL DONE!")

if __name__ == "__main__":
    main()