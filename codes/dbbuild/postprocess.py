#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
合并每类 FASTA → 精确去重（内容完全相同才视为重复）→ 分类别建 BLAST 库
提供“进度条监控”：先计数总记录，再按记录推进去重进度条；makeblastdb 显示阶段进度（开始/完成）

示例：
  python -m path_post.postprocess \
    --indir  /data1/liuyuxin/blast_seq_filter/original_data_lib \
    --dbdir  /data1/liuyuxin/blast_seq_filter/original_data_lib/db \
    --wrap 80 \
    --skip-count 0 \
    --tmpdb /data1/liuyuxin/blast_seq_filter/original_data_lib/dedup_index.sqlite
"""
import os, argparse, sqlite3
from typing import List, Tuple
from tqdm import tqdm
from .categories import CATEGORIES
from .fasta_utils import (
    iter_fasta_records, count_fasta_records, write_fasta_record,
    seq_md5, unique_header
)
from .external import need, run

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def list_fasta_files(cat_dir: str) -> List[str]:
    files = []
    for root, _, fs in os.walk(cat_dir):
        for fn in fs:
            if fn.lower().endswith((".fasta",".fa",".fna",".fasta.gz",".fa.gz",".fna.gz")):
                files.append(os.path.join(root, fn))
    return files

def count_records_all(files: List[str]) -> int:
    """第一阶段：统计所有文件总记录数，用于进度条上限"""
    total = 0
    for fp in tqdm(files, desc="[count] files", unit="file"):
        total += count_fasta_records(fp)
    return total

def dedup_and_merge(cat: str, files: List[str], out_fa: str, wrap: int, tmpdb: str,
                    total_records: int, update_every: int = 1000) -> Tuple[int,int]:
    """
    第二阶段：按记录推进去重进度条。
    返回 (kept, total_processed)
    """
    if not files:
        return 0, 0

    # sqlite 做“已见”索引（md5为主键），大规模更稳
    conn = sqlite3.connect(tmpdb)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("CREATE TABLE IF NOT EXISTS seen (md5 TEXT PRIMARY KEY)")
    conn.commit()

    used_ids = set()
    kept = 0; proc = 0
    ensure_dir(os.path.dirname(out_fa))
    tmp_out = out_fa + ".part"

    with open(tmp_out, "w", encoding="utf-8") as fo, \
         tqdm(total=total_records if total_records>0 else None,
              desc=f"[{cat}] dedup", unit="rec", mininterval=0.2) as bar:
        for fp in files:
            for h, s in iter_fasta_records(fp):
                proc += 1
                m = seq_md5(s)
                try:
                    conn.execute("INSERT INTO seen(md5) VALUES (?)", (m,))
                    uh = unique_header(h, used_ids)
                    write_fasta_record(fo, uh, s, wrap=wrap)
                    kept += 1
                except sqlite3.IntegrityError:
                    pass
                # 进度条推进与摘要（避免每条都 set_postfix 造成开销）
                bar.update(1)
                if (proc % update_every) == 0:
                    bar.set_postfix(kept=kept, dup=proc-kept, refresh=True)
        bar.set_postfix(kept=kept, dup=proc-kept, refresh=True)

    conn.commit(); conn.close()
    os.replace(tmp_out, out_fa)  # 原子替换
    return kept, proc

def make_db(in_fa: str, db_prefix: str, title: str, cat: str):
    need("makeblastdb")
    print(f"[{cat}] makeblastdb -> {db_prefix}.*")
    run([
        "makeblastdb",
        "-in", in_fa,
        "-dbtype", "nucl",
        "-parse_seqids",
        "-blastdb_version", "5",
        "-title", title,
        "-out", db_prefix
    ])
    print(f"[{cat}] makeblastdb done -> {db_prefix}.*")

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, help="下载好的根目录（包含各 category 子目录）")
    ap.add_argument("--dbdir", required=True, help="makeblastdb 输出目录")
    ap.add_argument("--categories", nargs="*", default=CATEGORIES)
    ap.add_argument("--wrap", type=int, default=80, help="FASTA 输出换行宽度，0=不换行")
    ap.add_argument("--tmpdb", default=":memory:", help="sqlite 去重索引（:memory: 用内存；给路径则落盘）")
    ap.add_argument("--keep-merged", type=int, default=1, help="建库后是否保留合并去重后的 fasta（1=保留，0=删除）")
    ap.add_argument("--skip-count", type=int, default=0,
                    help="跳过预计数：总量未知时也会显示计数进度（但无百分比/ETA）")
    ap.add_argument("--update-every", type=int, default=1000,
                    help="进度条摘要每处理多少条刷新一次")
    return ap.parse_args()

def main():
    args = parse_args()
    ensure_dir(args.dbdir)

    for cat in args.categories:
        cat_dir = os.path.join(args.indir, cat)
        if not os.path.isdir(cat_dir):
            print(f"[SKIP] no such category dir: {cat_dir}")
            continue

        files = list_fasta_files(cat_dir)
        if not files:
            print(f"[{cat}] no fasta files under {cat_dir}")
            continue

        merged_fa = os.path.join(args.indir, f"{cat}.fa")
        db_prefix = os.path.join(args.dbdir, cat)
        title = f"{cat} dedup exact"

        # 计数（可跳过）
        total = 0
        if int(args.skip_count) == 0:
            total = count_records_all(files)
            print(f"[{cat}] total records = {total:,}")

        # 去重合并（按记录推进进度条）
        kept, proc = dedup_and_merge(
            cat, files, merged_fa, args.wrap, args.tmpdb, total_records=total,
            update_every=int(args.update_every)
        )
        print(f"[{cat}] unique={kept:,} / processed={proc:,} / dup={proc-kept:,} -> {merged_fa}")

        # 建库
        if kept > 0:
            make_db(merged_fa, db_prefix, title, cat)
        else:
            print(f"[{cat}] nothing to build.")

        if not int(args.keep_merged):
            try: os.remove(merged_fa)
            except Exception: pass

if __name__ == "__main__":
    main()