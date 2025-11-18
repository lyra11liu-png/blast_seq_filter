#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_pipeline.py
主脚本：对一个目录下所有 BAM 文件执行：
    1) minimap2 + samtools 去除人类序列
    2) 对剩余 reads 在多个病原体 BLAST 库上做 BLASTn
    3) 输出每个样本的所有命中(all_hits)和物种统计表(summary_by_species)
       以及 taxonomy_report.tsv + 合并后的 taxonomy_summary.tsv

支持可选的多进程并行，以及简单的计算时间监控。
"""

from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple, Optional

from utils import default_logger_factory
from host_filter import bam_to_nonhuman_fasta
from blast_classify import run_multi_db_blast


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Host-depletion (human) + multi-DB pathogen BLAST classification pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--indir", required=True,
                   help="输入 BAM 文件所在目录（会自动寻找 *.bam）")
    p.add_argument("--bam", default=None,
                   help="只处理指定的单个 BAM 文件（可选），给出完整路径")
    p.add_argument("--outdir", required=True,
                   help="输出结果目录，每个样本一个子目录")

    # 人类索引
    p.add_argument("--human-index", required=True,
                   help="minimap2 用的人类参考索引 .mmi 路径")

    # 病原体 BLAST 数据库根目录 + 使用的前缀列表
    p.add_argument("--db-root", required=True,
                   help="病原体 BLAST 库所在目录（里面有 viruses.nhr/bacteria.nhr 等）")
    p.add_argument("--db-groups", default="viruses,bacteria,fungi,mycoplasma,archaea,protozoa,helminths",
                   help=("要使用的库前缀名，逗号分隔；"
                         "例如 'viruses,bacteria,fungi,mycoplasma'"))
    p.add_argument("--idmap-root", default=None,
                   help=("可选：分类建库生成的 <group>.idmap.tsv 所在目录；"
                         "若不指定，将默认使用 db-root 的上一级目录下的 clean/ 子目录"))

    # 外部工具
    p.add_argument("--minimap2", default="minimap2",
                   help="minimap2 可执行文件名称或路径")
    p.add_argument("--samtools", default="samtools",
                   help="samtools 可执行文件名称或路径")
    p.add_argument("--blastn", default="blastn",
                   help="blastn 可执行文件名称或路径")

    p.add_argument("--threads", type=int, default=8,
                   help="每个样本内部使用的线程数（传给 minimap2 / blastn）")
    p.add_argument("--preset", default="map-hifi",
                   help="minimap2 预设（PacBio HiFi 用 map-hifi，ONT 可用 map-ont 等）")

    p.add_argument("--parallel", action="store_true",
                   help="是否并行处理多个 BAM 样本")
    p.add_argument("--jobs", type=int, default=2,
                   help="并行时的进程数（样本数并行）")

    # ===== 新增：命中筛选阈值 =====
    p.add_argument("--min-pident", type=float, default=85.0,
                   help="最小比对相似度（百分数，例如 85.0）")
    p.add_argument("--min-align-len", type=int, default=80,
                   help="最小比对长度（bp），短于该长度的命中将被丢弃")
    p.add_argument("--max-hit-evalue", type=float, default=1e-20,
                   help="最大命中 e-value（越小越严格），大于该值的命中将被丢弃")

    return p.parse_args()


def find_bam_files(indir: str) -> List[Path]:
    d = Path(indir)
    if not d.is_dir():
        raise SystemExit(f"Input directory not found: {indir}")
    bams = sorted(d.glob("*.bam"))
    if not bams:
        raise SystemExit(f"No BAM files found in: {indir}")
    return bams


def build_db_map(db_root: str, groups_str: str) -> Dict[str, str]:
    """
    根据 db_root + 组名列表构造 db_map: group -> db_prefix
    """
    root = Path(db_root)
    if not root.is_dir():
        raise SystemExit(f"DB root directory not found: {db_root}")

    groups = [g.strip() for g in groups_str.split(",") if g.strip()]
    db_map: Dict[str, str] = {}

    for g in groups:
        prefix = root / g
        nhr = prefix.with_suffix(".nhr")
        nsq = prefix.with_suffix(".nsq")
        nal = prefix.with_suffix(".nal")
        nsq_parts = list(root.glob(f"{g}*.nsq"))

        if nhr.exists() or nsq.exists() or nal.exists() or nsq_parts:
            db_map[g] = str(prefix)
        else:
            print(f"[WARN] BLAST DB for group '{g}' not found under {db_root}, skip.")

    if not db_map:
        raise SystemExit("No valid BLAST DB found. Please check --db-root and --db-groups.")

    print("[INFO] Using BLAST DBs:")
    for g, p in db_map.items():
        print(f"  {g}: {p}")
    return db_map


def process_one_sample(
    bam_path: str,
    outdir: str,
    human_index: str,
    db_map: Dict[str, str],
    minimap2: str,
    samtools: str,
    blastn: str,
    threads: int,
    preset: str,
    idmap_root: Optional[str],
    min_pident: float,
    min_align_len: int,
    max_hit_evalue: float,
) -> Tuple[str, Dict[str, float]]:
    """
    单样本处理函数，用于主进程或子进程。
    """
    bam = Path(bam_path)
    sample = bam.stem

    sample_outdir = Path(outdir) / sample
    sample_outdir.mkdir(parents=True, exist_ok=True)

    log_path = sample_outdir / f"{sample}.log"
    logger = default_logger_factory(str(log_path))

    logger(f"[SAMPLE] {sample}")
    logger(f"[INPUT BAM] {bam_path}")

    # 1. 去除人类 reads，生成 non-human FASTA
    host_timings = bam_to_nonhuman_fasta(
        bam_path=str(bam),
        out_dir=str(sample_outdir),
        human_index=human_index,
        minimap2=minimap2,
        samtools=samtools,
        threads=threads,
        preset=preset,
        logger=logger,
    )
    nonhuman_fa = str(sample_outdir / f"{sample}.nonhuman.fasta")

    # 2. 对 non-human FASTA 在多个库上做 BLASTn（输出所有“通过筛选”的命中）
    blast_timings = run_multi_db_blast(
        query_fa=nonhuman_fa,
        out_prefix=str(sample_outdir / sample),
        db_map=db_map,
        blastn=blastn,
        threads=threads,
        idmap_root=idmap_root,
        logger=logger,
        min_pident=min_pident,
        min_align_len=min_align_len,
        max_hit_evalue=max_hit_evalue,
    )

    # 汇总 timings
    timings: Dict[str, float] = {}
    timings.update(host_timings)
    timings.update(blast_timings)
    return sample, timings


def write_runtime_summary(all_timings: Dict[str, Dict[str, float]], outdir: str) -> None:
    """
    把所有样本的耗时信息写到一个 TSV 汇总文件中，便于监控。
    """
    out_path = Path(outdir) / "runtime_summary.tsv"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("sample\tstep\tseconds\n")
        for sample, tdict in sorted(all_timings.items()):
            for step, sec in sorted(tdict.items()):
                f.write(f"{sample}\t{step}\t{sec:.2f}\n")
    print(f"[RUNTIME] Summary written to: {out_path}")


def merge_taxonomy_reports(outdir: str) -> None:
    """
    把每个样本子目录下的 *.taxonomy_report.tsv 合并成一个总表：
        outdir/taxonomy_summary.tsv
    方便直接在 Excel 里打开。
    """
    out_dir = Path(outdir)
    merged_path = out_dir / "taxonomy_summary.tsv"

    with open(merged_path, "w", encoding="utf-8") as fout:
        fout.write("Sample\tTaxonomy\n")
        for sub in sorted(out_dir.iterdir()):
            if not sub.is_dir():
                continue
            sample = sub.name
            tax_fp = sub / f"{sample}.taxonomy_report.tsv"
            if not tax_fp.is_file():
                continue
            with open(tax_fp, "r", encoding="utf-8") as fin:
                header = fin.readline()  # 跳过每个文件自己的表头
                for line in fin:
                    fout.write(line)

    print(f"[TAXONOMY] Merged taxonomy table written to: {merged_path}")


def main() -> None:
    args = parse_args()
    
    if args.bam is not None:
        bam_path = Path(args.bam)
        if not bam_path.is_file():
            raise SystemExit(f"--bam file not found: {args.bam}")
        bams = [bam_path]
        print(f"[INFO] Only process single BAM file: {args.bam}")
    else:
        bams = find_bam_files(args.indir)
        print(f"Found {len(bams)} BAM files in {args.indir}")
        
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # 构造 group -> db_prefix 映射（一次搞好，传给每个进程）
    db_map = build_db_map(args.db_root, args.db_groups)
    
    # 推断 idmap_root：若用户没显式指定，默认假设为 "<db-root>/../clean"
    if args.idmap_root is not None:
        idmap_root = args.idmap_root
    else:
        idmap_root = str((Path(args.db_root).resolve().parent / "clean"))
    print(f"[INFO] Using idmap_root: {idmap_root}")

    all_timings: Dict[str, Dict[str, float]] = {}

    if args.parallel:
        # 多进程并行按样本跑
        jobs = max(1, int(args.jobs))
        print(f"Running in parallel with {jobs} processes")
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            future_to_bam = {
                ex.submit(
                    process_one_sample,
                    str(bam),
                    args.outdir,
                    args.human_index,
                    db_map,
                    args.minimap2,
                    args.samtools,
                    args.blastn,
                    args.threads,
                    args.preset,
                    idmap_root,
                    args.min_pident,
                    args.min_align_len,
                    args.max_hit_evalue,
                ): bam
                for bam in bams
            }
            for fut in as_completed(future_to_bam):
                bam = future_to_bam[fut]
                try:
                    sample, timings = fut.result()
                    all_timings[sample] = timings
                    print(f"[DONE] sample {sample}")
                except Exception as e:
                    print(f"[ERROR] sample {bam.name} failed: {e}")
    else:
        # 顺序逐个样本跑，便于调试
        print("Running samples sequentially (no --parallel)")
        for bam in bams:
            try:
                sample, timings = process_one_sample(
                    str(bam),
                    args.outdir,
                    args.human_index,
                    db_map,
                    args.minimap2,
                    args.samtools,
                    args.blastn,
                    args.threads,
                    args.preset,
                    idmap_root,
                    args.min_pident,
                    args.min_align_len,
                    args.max_hit_evalue,
                )
                all_timings[sample] = timings
                print(f"[DONE] sample {sample}")
            except Exception as e:
                print(f"[ERROR] sample {bam.name} failed: {e}")

    # 写出整体耗时表 + 合并 tax 表
    write_runtime_summary(all_timings, args.outdir)
    merge_taxonomy_reports(args.outdir)


if __name__ == "__main__":
    main()
