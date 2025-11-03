#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
downloader.py
Pipeline:
  1. taxonomy.esummary determines output directory grouping
  2. construct term
  3. esearch_history retrieves count
  4. avoid redundant downloads
  5. resume interrupted downloads
  6. return a single result dictionary to the statics module
"""
import os, pathlib, time
from util import log, ensure_dir, sanitize, count_lines, count_fasta_records, count_fasta_records_safe
from taxonomy import map_group
from entrez import tax_esummary, esearch_history, efetch_fasta_history, esummary_taxmap

def build_term(taxid: str, subtree=True, genomic_only=True, refseq_only=False,
               slen_min=1000, slen_max=100000000) -> str:
    term = f"txid{taxid}[Subtree]" if subtree else f"txid{taxid}[Organism:exp]"
    if genomic_only: term += " AND biomol_genomic[PROP]"
    if refseq_only:  term += " AND srcdb_refseq[PROP]"
    if slen_min or slen_max:
        term += f" AND {int(slen_min)}:{int(slen_max)}[SLEN]"
    return term

def download_one(taxid: str, label: str, out_base: pathlib.Path,
                 subtree=True, genomic_only=True, refseq_only=False,
                 slen_min=1000, slen_max=100000000, page=5000,
                 resume=True, verify_skip=True, tolerance=0.02,
                 hybrid=False, virus_slen_min=200, microbe_slen_min=1000) -> dict:
    t0 = time.time()
    
    # 1. Metadata & grouping
    meta  = tax_esummary(taxid)
    sci   = meta.get("scientificname") or label or f"taxid{taxid}"
    group = map_group(meta, label)
    gdir  = ensure_dir(out_base / group)
    
    # 2. Filename & path
    safe_label = sanitize(sci)
    fa_path = gdir / f"{taxid}_{safe_label}.fasta"
    tm_path = gdir / f"{taxid}_{safe_label}.taxmap.tsv"
    fa_part = gdir / f"{taxid}_{safe_label}.fasta.part"
    tm_part = gdir / f"{taxid}_{safe_label}.taxmap.tsv.part"
    
    eff_refseq_only = refseq_only
    eff_slen_min    = slen_min
    if hybrid:
        if group == "viruses":
            eff_refseq_only = False
            eff_slen_min    = min(eff_slen_min, virus_slen_min)
        else:
            eff_refseq_only = False
            eff_slen_min    = max(eff_slen_min, microbe_slen_min)
        
    # 3. Construct term & search history
    term = build_term(taxid, subtree, genomic_only, eff_refseq_only, eff_slen_min, slen_max)
    expected, webenv, qk = esearch_history(term)
    
    # 4. If already complete, skip
    if fa_path.exists() and verify_skip:
        got = count_fasta_records(str(fa_path))
        if expected == 0 or got >= int(expected * (1 - tolerance)):
            log(f"[SKIP] txid {taxid} {sci} already have {got}, should have {expected}")
            return dict(taxid=taxid, group=group, file=str(fa_path), taxmap=str(tm_path),
                        count=got, expected=expected, skipped=True,
                        sci_name=sci, rank=meta.get("rank",""),
                        seconds=0.0, incomplete=False)
            
    log(f"[{taxid}] {sci} starts(expect {expected} lines)")
    
    # 5a. Update taxmap
    tm_start = 0
    if resume:
        if tm_path.exists() and count_lines(str(tm_path)) < expected and not tm_part.exists():
            os.replace(str(tm_path), str(tm_part))
        tm_start = count_lines(str(tm_part)) if tm_part.exists() else 0
    esummary_taxmap(webenv, qk, expected, out_part_path=str(tm_part), page=page, start=tm_start)
    if tm_part.exists():
        os.replace(str(tm_part), str(tm_path))
        
    # 5b. Add fasta
    fa_start = count_fasta_records_safe(str(fa_part)) if (resume and fa_part.exists()) else 0
    got = efetch_fasta_history(webenv, qk, expected, out_tmp_path=str(fa_part), page=page, start=fa_start)
    os.replace(str(fa_part), str(fa_path))
    
    # 6. Summary
    sec = time.time()-t0
    incomplete = (expected > 0 and got < int(expected * (1 - tolerance)))
    log(f"[Finish] txid {taxid} {sci} -> {fa_path.name} | lines={got}/{expected} | time {sec:.1f}s | group={group}")
    
    return dict(taxid=taxid, group=group, file=str(fa_path), taxmap=str(tm_path),
                count=got, expected=expected, skipped=False,
                sci_name=sci, rank=meta.get("rank",""),
                seconds=sec, incomplete=incomplete)