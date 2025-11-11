#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, os, sys, pandas as pd
from .downloader import download_all
from .classify import attach_or_infer_category
from . import ncbi as ncbi_mod

def parse_args():
    ap = argparse.ArgumentParser(description="Download pathogen sequences by TaxID with category folders.")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)

    # 长度阈值（0=不限制）
    ap.add_argument("--virus-slen-min", type=int, default=0)
    ap.add_argument("--microbe-slen-min", type=int, default=0)

    # 数据源策略
    ap.add_argument("--viruses-include-genbank", type=int, default=1)
    ap.add_argument("--microbes-include-genbank", type=int, default=0)  # 微生物默认 RefSeq-only
    ap.add_argument("--microbes-genbank-fallback", type=int, default=1,  # RefSeq 太少自动补 GenBank
                    help="微生物 RefSeq 记录数 < refseq-min-records 时，回退追加 GenBank")
    ap.add_argument("--refseq-min-records", type=int, default=1,
                    help="回退阈值：RefSeq 少于此数则补 GenBank")
    ap.add_argument("--only-complete", type=int, default=0)

    # 既有文件策略
    ap.add_argument("--skip-existing", type=int, default=0)
    ap.add_argument("--supplement-existing", type=int, default=1)

    # 网络与分页
    ap.add_argument("--page-size", type=int, default=10000)
    ap.add_argument("--timeout", type=int, default=60)
    ap.add_argument("--retries", type=int, default=6)
    ap.add_argument("--qps", type=float, default=8.0,  # 有 API key 时建议 6~10
                    help="全局限速（请求/秒），防止被限流")
    return ap.parse_args()

def main():
    args = parse_args()
    api_key = os.environ.get("NCBI_API_KEY", "").strip()
    if not api_key:
        print("WARN: 未设置 NCBI_API_KEY：更容易触发限流。", file=sys.stderr)

    # 初始化全局限速器（对 esearch/efetch 生效）
    ncbi_mod.init_rate_limiter(qps=float(args.qps))

    df = pd.read_csv(args.csv)
    df = attach_or_infer_category(df)
    keep = {"Bacteria","Fungi","Archaea","Mycoplasma","Viruses","Protozoa","Helminths"}
    df = df[df["category"].isin(keep)].copy()

    if "best_taxid" not in df.columns and "taxid" in df.columns:
        df = df.rename(columns={"taxid":"best_taxid"})
    if "best_name" not in df.columns and "name" in df.columns:
        df = df.rename(columns={"name":"best_name"})
    if "best_taxid" not in df.columns: raise SystemExit("CSV 里缺 best_taxid/taxid 列")
    if "best_name" not in df.columns: df["best_name"] = df.get("latin_input", "NA")

    df = df[~df["best_taxid"].isna()].copy()
    df["best_taxid"] = df["best_taxid"].astype(str).str.extract(r"(\d+)")[0]
    df = df.dropna(subset=["best_taxid"])

    config = {
        "virus_slen_min": int(args.virus_slen_min),
        "microbe_slen_min": int(args.microbe_slen_min),
        "viruses_include_genbank": bool(int(args.viruses_include_genbank)),
        "microbes_include_genbank": bool(int(args.microbes_include_genbank)),
        "microbes_genbank_fallback": bool(int(args.microbes_genbank_fallback)),
        "refseq_min_records": int(args.refseq_min_records),
        "only_complete": bool(int(args.only_complete)),
        "skip_existing": bool(int(args.skip_existing)),
        "supplement_existing": bool(int(args.supplement_existing)),
        "threads": int(args.threads),
        "page_size": int(args.page_size),
        "timeout": int(args.timeout),
        "retries": int(args.retries),
        "api_key": api_key,
        "outdir": args.outdir
    }
    download_all(df, config)

if __name__ == "__main__":
    main()
