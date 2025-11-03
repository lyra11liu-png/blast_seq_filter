#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
entrez.py
Thin packaging 4 E-utilities
"""

from typing import Tuple
from util import log
from http_client import http_get, http_post

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def tax_esummary(taxid: str) -> dict:
    """Read taxonomy esummary."""
    r = http_get(f"{EUTILS}/esummary.fcgi", {
        "db":"taxonomy", "id":taxid, "retmode":"json"
    })
    d = r.json().get("result", {}).get(str(taxid), {})
    return {
        "taxid": str(d.get("taxid", "")),
        "scientificname": d.get("scientificname", ""),
        "rank": d.get("rank", ""),
        "division": (d.get("division") or ""),
        "lineage": (d.get("lineage") or "")
    }
    
def esearch_history(term: str) -> Tuple[int, str, str]:
    r = http_get(f"{EUTILS}/esearch.fcgi", {
        "db":"nuccore", "term":term, "retmode":"json",
        "usehistory":"y", "retmax":0
    })
    er = r.json().get("esearchresult", {})
    return int(er.get("count","0")), er.get("webenv",""), er.get("querykey","")

def efetch_fasta_history(webenv: str, query_key: str, count: int,
                         out_tmp_path: str, page: int=5000, start: int=0) -> int:
    """
    Pull fasta files by history page & write to out_tmp_path.
    """
    done, retstart = start, start
    mode = "ab" if start > 0 else "wb"
    with open(out_tmp_path, mode) as fout:
        while retstart < count:
            retmax = min(page, count - retstart)
            rr = http_post(f"{EUTILS}/efetch.fcgi", {
                "db":"nuccore","retmode":"text","rettype":"fasta",
                "webenv":webenv,"query_key":query_key,
                "retstart":retstart,"retmax":retmax
            }, stream=True)
            for chunk in rr.iter_content(chunk_size=1<<15):
                if chunk: fout.write(chunk)
            retstart += retmax
            done += retmax
            log(f"  efetch: {done}/{count} (this page {retmax})")
    return done

def esummary_taxmap(webenv: str, query_key: str, count: int,
                    out_part_path: str, page: int=5000, start: int=0) -> int:
    done, retstart = start, start
    mode = "a" if start > 0 else "w"
    with open(out_part_path, mode, encoding="utf-8") as out:
        while retstart < count:
            retmax = min(page, count - retstart)
            r = http_post(f"{EUTILS}/esummary.fcgi", {
                "db":"nuccore","retmode":"json",
                "webenv":webenv,"query_key":query_key,
                "retstart":retstart,"retmax":retmax
            })
            rs = r.json().get("result", {})
            for k, v in rs.items():
                if k == "uids": continue
                acc = v.get("caption")
                tx = v.get("taxid")
                if acc and tx:
                    out.write(f"{acc}\t{tx}\n")
                    done += 1
            retstart += retmax
    return done