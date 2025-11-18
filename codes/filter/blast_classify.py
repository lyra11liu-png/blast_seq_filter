#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
blast_classify.py  —— 多库 BLAST 分类（输出所有命中）

对非人 reads FASTA 依次在多个 BLAST 库上跑 blastn：
    - 保留各库原始 BLAST 表：{out_prefix}.<group>.blast.tsv
    - 合并所有库命中到一个大表：{out_prefix}.all_hits.tsv
    - 基于所有命中按 group + kingdom + species 统计：
      {out_prefix}.summary_by_species.tsv 和 {out_prefix}.taxonomy_report.tsv
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, Callable, Tuple, Set, List
import re   # <<< 新增：用于清洗 147573_Piedraia_hortae 这类名字

from utils import run_cmd, Timer


# ======================================================================
# 新增：把 147573_Piedraia_hortae -> Piedraia hortae 之类的“干净物种名”
# ======================================================================
def _normalize_species_name(raw_unit: str, orig_header: str, fallback: str) -> str:
    """
    根据 idmap 中的 unit / orig_header / sseqid 生成“干净”的物种名。

    优先级：
      1) unit（例如 147573_Piedraia_hortae）
      2) orig_header 的第一个 token
      3) fallback（通常是 sseqid）

    规则：
      - 若形如 "数字_后缀"，去掉前面的数字和下划线；
      - 把剩余的 "_" 换成空格。
    """
    name = raw_unit or (orig_header.split()[0] if orig_header else fallback)
    name = name.strip()
    if not name:
        return "Unknown"

    # 去掉数字前缀：147573_Piedraia_hortae -> Piedraia_hortae
    m = re.match(r"^(\d+)_([A-Za-z].*)$", name)
    if m:
        name = m.group(2)

    # 把剩余的 "_" 改成空格：Piedraia_hortae -> Piedraia hortae
    return name.replace("_", " ")


def run_multi_db_blast(
    query_fa: str,
    out_prefix: str,
    db_map: Dict[str, str],
    blastn: str = "blastn",
    threads: int = 8,
    task: str = "megablast",
    evalue: float = 1e-10,        # 传给 blastn 的 -evalue 上限
    max_target_seqs: int = 1,
    idmap_root: Optional[str] = None,
    logger: Optional[Callable[[str], None]] = None,
    top_n: int = 20,
    # ====== 新增：命中筛选阈值 ======
    min_pident: float = 85.0,     # 最小相似度（%）
    min_align_len: int = 80,      # 最短对齐长度（bp）
    max_hit_evalue: float = 1e-20 # 命中 e-value 上限
) -> Dict[str, float]:
    """
    在多个 BLAST 库上进行比对，并输出所有“符合条件”的命中物种。
    只有同时满足：
        pident >= min_pident
        align_len >= min_align_len
        evalue <= max_hit_evalue
    的命中才会进入 all_hits / summary / taxonomy。
    """
    query = Path(query_fa)
    out_prefix_p = Path(out_prefix)
    out_dir = out_prefix_p.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    timings: Dict[str, float] = {}
    
    # ==== 预加载 idmap: group -> {short_id -> (category, unit, orig_header)} ====
    idmaps: Dict[str, Dict[str, Tuple[str, str, str]]] = {}
    if idmap_root is not None:
        id_root = Path(idmap_root)
        if not id_root.is_dir():
            if logger:
                logger(f"[WARN] idmap_root directory not found: {id_root}, ignore idmap.")
        else:
            for group in db_map.keys():
                idmap_path = id_root / f"{group}.idmap.tsv"
                if not idmap_path.is_file():
                    if logger:
                        logger(f"[WARN] idmap for group '{group}' not found: {idmap_path}")
                    continue
                mapping: Dict[str, Tuple[str, str, str]] = {}
                with open(idmap_path, "r", encoding="utf-8") as f:
                    header = f.readline()  # 跳过表头: short_id category unit file orig_header
                    for line in f:
                        line = line.rstrip("\n")
                        if not line:
                            continue
                        parts = line.split("\t")
                        if len(parts) < 5:
                            continue
                        short_id, category, unit, src_file, orig_header = parts[:5]
                        mapping[short_id] = (category, unit, orig_header)
                idmaps[group] = mapping
            if logger and idmaps:
                logger("[INFO] loaded idmap for groups: " + ", ".join(sorted(idmaps.keys())))

    # 合并所有“通过筛选”的命中
    all_hits_tsv = f"{out_prefix}.all_hits.tsv"
    # 统计 (group, kingdom, species) -> set(read_id)
    species_reads: Dict[Tuple[str, str, str], Set[str]] = {}

    # 这里的 outfmt 不再要 sskingdoms / sscinames，而是自己用 idmap 填物种
    fmt = "6 qseqid sseqid pident length evalue bitscore"

    with open(all_hits_tsv, "w", encoding="utf-8") as fout_all:
        fout_all.write(
            "read_id\tgroup\tkingdom\tspecies\tpident\talign_len\tevalue\tbitscore\n"
        )

        # ========= 1. 依次在多个库上跑 BLAST，边跑边合并 =========
        for group, db_prefix in db_map.items():
            blast_tsv = f"{out_prefix}.{group}.blast.tsv"
            cmd = [
                blastn,
                "-task", task,
                "-db", db_prefix,
                "-query", str(query),
                "-outfmt", fmt,
                "-max_target_seqs", str(max_target_seqs),
                "-num_threads", str(threads),
                "-evalue", str(evalue),
                "-out", blast_tsv,
            ]
            step_name = f"blast_{group}"
            with Timer(step_name, logger) as t:
                run_cmd(cmd, logger=logger, shell=False)
            timings[step_name] = t.elapsed

            n_lines = 0
            with open(blast_tsv, "r", encoding="utf-8") as fin:
                for line in fin:
                    line = line.strip()
                    if not line:
                        continue
                    n_lines += 1
                    fields = line.split("\t")
                    if len(fields) < 6:
                        raise ValueError(
                            f"Unexpected BLAST outfmt columns ({len(fields)}) in {blast_tsv}: {line}"
                        )
                    qseqid, sseqid, pident_s, length_s, evalue_s, bitscore_s = fields[:6]

                    # 数值化
                    try:
                        pident = float(pident_s)
                    except ValueError:
                        pident = 0.0
                    try:
                        alen = int(length_s)
                    except ValueError:
                        alen = 0
                    try:
                        evalue_f = float(evalue_s)
                    except ValueError:
                        evalue_f = 1.0
                    try:
                        bitscore = float(bitscore_s)
                    except ValueError:
                        bitscore = 0.0

                    # ====== 新增：命中筛选 ======
                    # 例如：默认 min_pident=85，min_align_len=80，max_hit_evalue=1e-20
                    # 会把你之前 61 bp 的短命中 / e-value 相对较大的命中过滤掉。
                    if (pident < min_pident) or (alen < min_align_len) or (evalue_f > max_hit_evalue):
                        continue
                        
                    # === 用 idmap / group 补充 kingdom + 物种名（已经清洗） ===
                    if idmaps:
                        m = idmaps.get(group)
                        if m and sseqid in m:
                            cat2, unit2, orig2 = m[sseqid]
                            kingdom = (cat2 or group or "Unknown")   # 大类：bacteria/viruses...
                            species = _normalize_species_name(unit2, orig2, sseqid)
                        else:
                            kingdom = group
                            species = _normalize_species_name("", "", sseqid)
                    else:
                        # 没有 idmap，就退化成 group + sseqid 清洗
                        kingdom = group
                        species = _normalize_species_name("", "", sseqid)

                    # 写入 all_hits 表（注意：这里已经是“干净物种名”了）
                    fout_all.write(
                        f"{qseqid}\t{group}\t{kingdom}\t{species}\t"
                        f"{pident:.2f}\t{alen}\t{evalue_f:.2e}\t{bitscore:.1f}\n"
                    )

                    # 用 set 记录每个物种命中的 read_id，用于后面统计 unique_reads
                    key = (group, kingdom, species)
                    species_reads.setdefault(key, set()).add(qseqid)

            if logger:
                logger(f"[INFO] lines in {blast_tsv}: {n_lines}")

    # ========= 2. 统计物种：summary_by_species + taxonomy_report =========
    summary_tsv   = f"{out_prefix}.summary_by_species.tsv"
    taxonomy_tsv  = f"{out_prefix}.taxonomy_report.tsv"
    sample_name   = out_prefix_p.name
    
    from collections import defaultdict
    group_counts: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
    
    with Timer("summarize", logger) as t:
        # 根据 unique_reads 排序（命中的 read 数，从高到低）
        sorted_items = sorted(
            species_reads.items(),
            key=lambda kv: (-len(kv[1]), kv[0][0], kv[0][1], kv[0][2])
        )

        # 1) summary_by_species.tsv
        with open(summary_tsv, "w", encoding="utf-8") as fout:
            fout.write("group\tkingdom\tspecies\tunique_reads\n")
            for (group, kingdom, species), read_ids in sorted_items:
                cnt = len(read_ids)
                fout.write(f"{group}\t{kingdom}\t{species}\t{cnt}\n")
                group_counts[group].append((species, cnt))
        
        # 2) taxonomy_report.tsv（每个样本一行，类似 Kraken 风格汇总）
        GROUP_LABEL = {
            "bacteria":   "Bacteria",
            "viruses":    "Viruses",
            "fungi":      "Fungi",
            "mycoplasma": "Mycoplasma",
            "archaea":    "Archaea",
            "protozoa":   "Protozoa",
            "helminths":  "Helminths",
        }
        
        with open(taxonomy_tsv, "w", encoding="utf-8") as ftax:
            ftax.write("Sample\tTaxonomy\n")
            # 按 db_map 里库的顺序来输出，这样 Bacteria / Fungi / Viruses 顺序稳定
            for g in db_map.keys():
                label = GROUP_LABEL.get(g, g.title())
                species_list = group_counts.get(g)

                if not species_list:
                    line = f"{label}: No results."
                else:
                    parts = []
                    for sp, cnt in species_list:
                        # 这里的 sp 已经是“干净物种名”，例如 "Listeria ivanovii"
                        parts.append(f"{sp}({cnt})")
                    line = f"{label}: " + " | ".join(parts)

                # 写到 taxonomy_report.tsv
                ftax.write(f"{sample_name}\t{line}\n")
                # 同步打印到日志 / 终端
                if logger:
                    logger(f"[TAXONOMY] {sample_name}\t{line}")

        # 额外打印 top N 物种
        if logger:
            total_hits = sum(len(v) for v in species_reads.values())
            logger(f"[INFO] species with >=1 hit: {len(species_reads)}")
            logger(f"[INFO] total (group,species) hit sets: {total_hits}")
            logger("[TOP SPECIES] group\tkingdom\tspecies\tunique_reads")
            for (group, kingdom, species), read_ids in sorted_items[:top_n]:
                logger(f"[TOP] {group}\t{kingdom}\t{species}\t{len(read_ids)}")

    timings["summarize"] = t.elapsed
    return timings
