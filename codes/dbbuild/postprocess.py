#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
合并每类 FASTA → 精确去重（序列内容完全相同才算重复）→ 分类别建 BLAST 库
带“总进度 + 文件级进度条”，并兼容老版 SQLite。
"""

import os, argparse, sqlite3
from typing import List, Tuple
from tqdm import tqdm

from .categories import CATEGORIES
from .fasta_utils import (
    iter_fasta_records, write_fasta_record, md5_blob, unique_header
)
from .external import need, run

SQLITE_VER = tuple(int(x) for x in sqlite3.sqlite_version.split("."))
USE_RETURNING = SQLITE_VER >= (3, 35, 0)

def ensure_dir(p): os.makedirs(p, exist_ok=True)

def list_fasta_files(cat_dir: str) -> List[str]:
    files = []
    for root, _, fs in os.walk(cat_dir):
        for fn in fs:
            if fn.lower().endswith((".fasta",".fa",".fna",".fasta.gz",".fa.gz",".fna.gz")):
                files.append(os.path.join(root, fn))
    files.sort()
    return files

def count_records(path: str) -> int:
    n = 0
    if path.endswith(".gz"):
        import gzip
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if ln.startswith(">"): n += 1
    else:
        with open(path, "rt", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if ln.startswith(">"): n += 1
    return n

def count_records_all(files: List[str]) -> int:
    total = 0
    for fp in tqdm(files, desc="[count] files", unit="file"):
        total += count_records(fp)
    return total

def _open_sqlite(tmpdb: str) -> sqlite3.Connection:
    # autocommit 模式，手动 BEGIN/COMMIT 组成小事务，避免巨大事务卡顿
    conn = sqlite3.connect(tmpdb, isolation_level=None)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    try:
        conn.execute("PRAGMA mmap_size=536870912")  # 512MB 映射（可选）
    except Exception:
        pass
    conn.execute("CREATE TABLE IF NOT EXISTS seen (md5 BLOB PRIMARY KEY) WITHOUT ROWID")
    return conn

def relpath(p, base):
    try:
        return os.path.relpath(p, base)
    except Exception:
        return p

def dedup_and_merge(cat: str, files: List[str], out_fa: str, wrap: int, tmpdb: str,
                    total_records: int, update_every: int = 5000,
                    commit_every: int = 100_000, show_file_bar: bool = True,
                    base_dir_for_rel: str = "") -> Tuple[int,int]:
    """核心去重流程：总进度 + 文件级进度条"""
    if not files:
        return 0, 0

    conn = _open_sqlite(tmpdb)
    cur = conn.cursor()

    used_ids = set()
    kept = 0; proc = 0; since_commit = 0

    ensure_dir(os.path.dirname(out_fa))
    tmp_out = out_fa + ".part"

    # 总进度条
    main_bar = tqdm(total=total_records if total_records>0 else None,
                    desc=f"[{cat}] dedup", unit="rec", mininterval=0.2, position=0, leave=True)

    with open(tmp_out, "w", encoding="utf-8") as fo:
        cur.execute("BEGIN")
        for fp in files:
            # 文件级进度条
            f_total = count_records(fp) if total_records > 0 else 0
            f_desc = f"[{cat}] {relpath(fp, base_dir_for_rel) if base_dir_for_rel else os.path.basename(fp)}"
            f_bar = None
            if show_file_bar:
                f_bar = tqdm(total=f_total if f_total>0 else None,
                             desc=f_desc[:100], unit="rec", mininterval=0.2, position=1, leave=False)

            for h, s in iter_fasta_records(fp):
                proc += 1
                dig = md5_blob(s)

                if USE_RETURNING:
                    cur.execute(
                        "INSERT INTO seen(md5) VALUES (?) "
                        "ON CONFLICT(md5) DO NOTHING RETURNING 1",
                        (sqlite3.Binary(dig),)
                    )
                    is_new = (cur.fetchone() is not None)
                else:
                    cur.execute("INSERT OR IGNORE INTO seen(md5) VALUES (?)",
                                (sqlite3.Binary(dig),))
                    cur.execute("SELECT changes()")
                    is_new = (cur.fetchone()[0] == 1)

                if is_new:
                    write_fasta_record(fo, unique_header(h, used_ids), s, wrap=wrap)
                    kept += 1

                since_commit += 1
                if since_commit >= commit_every:
                    cur.execute("COMMIT")
                    cur.execute("BEGIN")
                    since_commit = 0

                # 推进进度
                main_bar.update(1)
                if f_bar: f_bar.update(1)
                if (proc % update_every) == 0:
                    main_bar.set_postfix(kept=kept, dup=proc-kept, sqlite=f"{SQLITE_VER}", refresh=True)

            if f_bar:
                f_bar.set_postfix(file_done=True, refresh=True)
                f_bar.close()

        cur.execute("COMMIT")
        main_bar.set_postfix(kept=kept, dup=proc-kept, sqlite=f"{SQLITE_VER}", refresh=True)

    main_bar.close()
    conn.close()
    os.replace(tmp_out, out_fa)
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
    ap.add_argument("--indir", required=True)
    ap.add_argument("--dbdir", required=True)
    ap.add_argument("--categories", nargs="*", default=CATEGORIES)
    ap.add_argument("--wrap", type=int, default=80)
    ap.add_argument("--tmpdb", default=":memory:", help=":memory:（最快）或 /dev/shm/xxx.sqlite 或 某磁盘路径")
    ap.add_argument("--keep-merged", type=int, default=1)
    ap.add_argument("--skip-count", type=int, default=0)
    ap.add_argument("--update-every", type=int, default=5000)
    ap.add_argument("--commit-every", type=int, default=100000)
    ap.add_argument("--no-file-progress", action="store_true", help="不显示每个文件的进度条")
    return ap.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.dbdir, exist_ok=True)

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

        total = 0
        if int(args.skip_count) == 0:
            total = count_records_all(files)
            print(f"[{cat}] total records = {total:,}")

        kept, proc = dedup_and_merge(
            cat, files, merged_fa, args.wrap, args.tmpdb, total_records=total,
            update_every=int(args.update_every),
            commit_every=int(args.commit_every),
            show_file_bar=(not args.no_file_progress),
            base_dir_for_rel=args.indir
        )
        print(f"[{cat}] unique={kept:,} / processed={proc:,} / dup={proc-kept:,} -> {merged_fa}")

        if kept > 0:
            make_db(merged_fa, db_prefix, title, cat)
        else:
            print(f"[{cat}] nothing to build.")

        if not int(args.keep_merged):
            try: os.remove(merged_fa)
            except Exception: pass

if __name__ == "__main__":
    main()
