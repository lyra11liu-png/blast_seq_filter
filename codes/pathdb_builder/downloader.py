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
from util import log, ensure_dir, sanitize, count_lines, count_fasta_records, count_fasta_records_safe, reveal_hidden_part
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
                 hybrid=False, virus_slen_min=200, microbe_slen_min=1000,
                 microbe_refseq_fallback=200) -> dict:
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
    
    fa_part = reveal_hidden_part(fa_part)
    tm_part = reveal_hidden_part(tm_part)
    
    eff_refseq_only = refseq_only
    eff_slen_min    = slen_min
    if hybrid:
        if group == "viruses":
            eff_refseq_only = False
            eff_slen_min    = min(eff_slen_min, virus_slen_min)
        else:
            eff_refseq_only = True
            eff_slen_min    = max(eff_slen_min, microbe_slen_min)
        
    # 3. Construct term & search history
    term = build_term(taxid, subtree, genomic_only, eff_refseq_only, eff_slen_min, slen_max)
    expected, webenv, qk = esearch_history(term)
    log(f"[{taxid}] group={group} term={term}")
        
    if expected == 0:
        log(f"[{taxid}] No records match the term; skip creating FASTA/TaxMap.")
        # 清理历史遗留的空壳
        for p in (fa_part, tm_part, fa_path, tm_path):
            try:
                if os.path.exists(p) and os.path.getsize(p) == 0:
                    os.remove(p)
            except Exception:
                pass
        return dict(
            taxid=taxid, group=group, file="", taxmap="",
            count=0, expected=0, skipped=False,
            sci_name=sci, rank=meta.get("rank", ""),
            seconds=0.0, incomplete=False
        )

    def _cap_for(g: str) -> int:
        # 先看组专属，再看全局；环境变量留空或0表示不限
        envs = [f"MAX_TOTAL_{g.upper()}", "MAX_TOTAL"]
        for name in envs:
            v = (os.getenv(name) or "").strip()
            if v.isdigit() and int(v) > 0:
                return int(v)
        return 0

    CAP = _cap_for(group)
    if CAP and expected > CAP:
        log(f"[cap] limit total from {expected} to {CAP} ({group})")
        expected = CAP
    
    if hybrid and group != "viruses" and eff_refseq_only and expected < microbe_refseq_fallback:
        log(f"[{taxid}] RefSeq-only 命中 {expected} < {microbe_refseq_fallback}，回退到 GenBank+RefSeq")
        eff_refseq_only = False
        term = build_term(taxid, subtree, genomic_only, eff_refseq_only, eff_slen_min, slen_max) + \
           " AND NOT metagenome[Title] AND NOT environmental[Title] AND NOT uncultured[Title]"
        log(f"[{taxid}] fallback term={term}")
        expected, webenv, qk = esearch_history(term)
        if expected == 0:
            log(f"[{taxid}] fallback term hits 0; skip FASTA/TaxMap.")
            for p in (fa_part, tm_part, fa_path, tm_path):
                try:
                    if os.path.exists(p) and os.path.getsize(p) == 0:
                        os.remove(p)
                except Exception:
                    pass
            return dict(
                taxid=taxid, group=group, file="", taxmap="",
                count=0, expected=0, skipped=False,
                sci_name=sci, rank=meta.get("rank", ""),
                seconds=0.0, incomplete=False
            )
        CAP = _cap_for(group)
        if CAP and expected > CAP:
            log(f"[cap] limit total from {expected} to {CAP} ({group})")
            expected = CAP

        
    if group == "viruses" and not eff_refseq_only:
        term = term + " AND NOT wgs[Filter] AND NOT tsa[Filter] AND NOT mRNA[Filter]"
        expected, webenv, qk = esearch_history(term)
        log(f"[{taxid}] viruses term refined: {term} (count={expected})")
        if expected == 0:
            log(f"[{taxid}] fallback term hits 0; skip FASTA/TaxMap.")
            for p in (fa_part, tm_part, fa_path, tm_path):
                try:
                    if os.path.exists(p) and os.path.getsize(p) == 0:
                        os.remove(p)
                except Exception:
                    pass
            return dict(
                taxid=taxid, group=group, file="", taxmap="",
                count=0, expected=0, skipped=False,
                sci_name=sci, rank=meta.get("rank", ""),
                seconds=0.0, incomplete=False
            )
        CAP = _cap_for(group)
        if CAP and expected > CAP:
            log(f"[cap] limit total from {expected} to {CAP} ({group})")
            expected = CAP

    
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
    
    # ===== 5a. Update taxmap =====
    SKIP_TAXMAP = os.getenv("SKIP_TAXMAP", "0") == "1"
    if SKIP_TAXMAP:
        log(f"[{taxid}] skip esummary/taxmap stage (we will rebuild taxmap from FASTA later)")
    else:
        tm_start = 0
        if resume:
            if tm_path.exists() and count_lines(str(tm_path)) < expected and not tm_part.exists():
                os.replace(str(tm_path), str(tm_part))
            tm_start = count_lines(str(tm_part)) if tm_part.exists() else 0

        attempts, cur_page = 0, page
        while not SKIP_TAXMAP:
            try:
                esummary_taxmap(webenv, qk, expected, out_part_path=str(tm_part),
                                page=cur_page, start=tm_start, tag=f"[{taxid} {safe_label}]")
                break
            except Exception as e:
                attempts += 1
                if attempts > 5:
                    raise
                time.sleep(1.5 * attempts)
                expected, webenv, qk = esearch_history(term)
                tm_start = count_lines(str(tm_part)) if tm_part.exists() else 0
                cur_page = max(500, cur_page // 2)
                log(f"[{taxid}] esummary retry#{attempts} with page={cur_page}, start={tm_start} due to: {e}")

        if (not SKIP_TAXMAP) and tm_part.exists():
            os.replace(str(tm_part), str(tm_path))
        
    # ===== 5b. Add fasta =====
    fa_start = count_fasta_records_safe(str(fa_part)) if (resume and fa_part.exists()) else 0
    attempts, cur_page = 0, page
    while True:
        try: 
            got = efetch_fasta_history(
                webenv, qk, expected,
                out_tmp_path=str(fa_part), 
                page=cur_page, start=fa_start,
                tag=f"[{taxid} {safe_label}]"
            )
            break
        except Exception as e:
            attempts += 1
            if attempts > 5:
                raise
            time.sleep(1.5 * attempts)
            expected, webenv, qk = esearch_history(term)
            fa_start = count_fasta_records_safe(str(fa_part)) if (resume and fa_part.exists()) else 0
            cur_page = max(500, cur_page // 2)
            log(f"[{taxid}] efetch retry#{attempts} page={cur_page}, start={fa_start} due to: {e}")

    # 只有 .part 有内容才“晋升”为最终 .fasta；否则清理空壳，避免 0B
    if os.path.exists(fa_part):
        if os.path.getsize(fa_part) > 0:
            os.replace(str(fa_part), str(fa_path))
        else:
            try:
                os.remove(fa_part)
            except Exception:
                pass
            try:
                if os.path.exists(fa_path) and os.path.getsize(fa_path) == 0:
                    os.remove(fa_path)
            except Exception:
                pass
            
    # 6. Summary
    sec = time.time()-t0
    incomplete = (expected > 0 and got < int(expected * (1 - tolerance)))
    log(f"[Finish] txid {taxid} {sci} -> {fa_path.name} | lines={got}/{expected} | time {sec:.1f}s | group={group}")
    
    return dict(taxid=taxid, group=group, file=str(fa_path), taxmap=str(tm_path),
                count=got, expected=expected, skipped=False,
                sci_name=sci, rank=meta.get("rank",""),
                seconds=sec, incomplete=incomplete)