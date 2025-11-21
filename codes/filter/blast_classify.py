#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
blast_classify.py

For non-human reads, run blastn sequentially across multiple blast database:
    - Retain original blast tables 4 each database
    - Merge all hits from all databases into one large table
    - Generate statistics based on all hits grouped by group & kingdom & species

Author: lyra Liu
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, Callable, Tuple, Set, List
import re
from utils import run_cmd, Timer


def _normalize_species_name(raw_unit: str, orig_header: str, fallback: str) -> str:
    """
    Generate clean species names based on the unit, orig_header, sseqid fields in the idmap.

    Priority:
      1) unit (e.g., 147573_Piedraia_hortae)
      2) first token of orig_header
      3) fallback (sseqid)

    Rules:
      - If the pattern is "numer_fuffix", remove the preceding number & underscore
      - Replace the remaining "_" with spaces
    """
    name = raw_unit or (orig_header.split()[0] if orig_header else fallback)
    name = name.strip()
    if not name:
        return "Unknown"

    # Remove the numerical prefix: 147573_Piedraia_hortae -> Piedraia_hortae
    m = re.match(r"^(\d+)_([A-Za-z].*)$", name)
    if m:
        name = m.group(2)

    # Replace the remaining "_" with spaces: Piedraia_hortae -> Piedraia hortae
    return name.replace("_", " ")


def run_multi_db_blast(
    query_fa: str,
    out_prefix: str,
    db_map: Dict[str, str],
    blastn: str = "blastn",
    threads: int = 8,
    task: str = "megablast",
    evalue: float = 1e-10,
    max_target_seqs: int = 1,
    idmap_root: Optional[str] = None,
    logger: Optional[Callable[[str], None]] = None,
    top_n: int = 20,
    # ===== Hit filter threshold ====
    min_pident: float = 85.0,     # minimum similarity (%)
    min_align_len: int = 80,      # minimum alignment length (bp)
    max_hit_evalue: float = 1e-20 # hit e-value upper limit
) -> Dict[str, float]:
    """
    Perform alignments against multiple blast databases output all "qualifying" hits species.
    Only hits satisfying all three conditions simultaneously will be included in the final report.
    """
    query = Path(query_fa)
    out_prefix_p = Path(out_prefix)
    out_dir = out_prefix_p.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    timings: Dict[str, float] = {}
    
    # ==== Preload idmap: group -> {short_id -> (category, unit, orig_header)} ====
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
                    header = f.readline()
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

    # Merge all hits that passed the screening
    all_hits_tsv = f"{out_prefix}.all_hits.tsv"

    species_reads: Dict[Tuple[str, str, str], Set[str]] = {}

    fmt = "6 qseqid sseqid pident length evalue bitscore"

    with open(all_hits_tsv, "w", encoding="utf-8") as fout_all:
        fout_all.write(
            "read_id\tgroup\tkingdom\tspecies\tpident\talign_len\tevalue\tbitscore\n"
        )

        # ===== 1. Run blast on multiple databases sequentially, merging results as they r processed. =====
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

                    # Numerical representation
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

                    if (pident < min_pident) or (alen < min_align_len) or (evalue_f > max_hit_evalue):
                        continue
                        
                    if idmaps:
                        m = idmaps.get(group)
                        if m and sseqid in m:
                            cat2, unit2, orig2 = m[sseqid]
                            kingdom = (cat2 or group or "Unknown")
                            species = _normalize_species_name(unit2, orig2, sseqid)
                        else:
                            kingdom = group
                            species = _normalize_species_name("", "", sseqid)
                    else:
                        kingdom = group
                        species = _normalize_species_name("", "", sseqid)

                    # Write to the all_hits table.
                    fout_all.write(
                        f"{qseqid}\t{group}\t{kingdom}\t{species}\t"
                        f"{pident:.2f}\t{alen}\t{evalue_f:.2e}\t{bitscore:.1f}\n"
                    )

                    # Use a set to record the read_id 4 each species hit, 4 later counting unique_reads.
                    key = (group, kingdom, species)
                    species_reads.setdefault(key, set()).add(qseqid)

            if logger:
                logger(f"[INFO] lines in {blast_tsv}: {n_lines}")

    # ==== 2. Statistical species: summary_by_species + taxonomy_report ====
    summary_tsv   = f"{out_prefix}.summary_by_species.tsv"
    taxonomy_tsv  = f"{out_prefix}.taxonomy_report.tsv"
    sample_name   = out_prefix_p.name
    
    from collections import defaultdict
    group_counts: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
    
    with Timer("summarize", logger) as t:
        # Sort by unique_reads (number of reads hit, descending order).
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
        
        # 2) taxonomy_report.tsv
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
            for g in db_map.keys():
                label = GROUP_LABEL.get(g, g.title())
                species_list = group_counts.get(g)

                if not species_list:
                    line = f"{label}: No results."
                else:
                    parts = []
                    for sp, cnt in species_list:
                        parts.append(f"{sp}({cnt})")
                    line = f"{label}: " + " | ".join(parts)

                # Write to taxonomy_report.tsv.
                ftax.write(f"{sample_name}\t{line}\n")
                # Synchronous printing to log.
                if logger:
                    logger(f"[TAXONOMY] {sample_name}\t{line}")

        if logger:
            total_hits = sum(len(v) for v in species_reads.values())
            logger(f"[INFO] species with >=1 hit: {len(species_reads)}")
            logger(f"[INFO] total (group,species) hit sets: {total_hits}")
            logger("[TOP SPECIES] group\tkingdom\tspecies\tunique_reads")
            for (group, kingdom, species), read_ids in sorted_items[:top_n]:
                logger(f"[TOP] {group}\t{kingdom}\t{species}\t{len(read_ids)}")

    timings["summarize"] = t.elapsed
    return timings
