#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
entrez.py
Thin packaging 4 E-utilities
"""

import time, requests
from urllib3.exceptions import ProtocolError, InvalidChunkLength
from typing import Tuple
from util import log
from http_client import http_get, http_post

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

class EutilsHistoryStale(RuntimeError):
    """E-utilities history(WebEnv) stale/broken or backend-side error."""
    pass

def _json_or_stale(response, what="esummary"):
    """Parse JSON; if非JSON/含error字段，则抛出可识别的异常，交给上层重建WebEnv。"""
    try:
        j = response.json()
    except Exception:
        txt = response.text[:200].replace("\n", " ")
        raise EutilsHistoryStale(f"{what} returned non-JSON or broken JSON: {txt}")
    if isinstance(j, dict) and ("error" in j):
        raise EutilsHistoryStale(f"{what} error: {j.get('error')}")
    return j

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
                         out_tmp_path: str, page: int=5000, start: int=0,
                         tag: str = "") -> int:
    """
    Pull fasta files by history page & write to out_tmp_path.
    tag: 便于多线程日志定位（例如 "[1824 Nocardia_asteroides] "）
    """
    if count <= start:
        if tag:
            log(f"{tag} efetch: {start}/{count} (nothing to do)")
        return start
    
    done, retstart = start, start
    mode = "ab" if start > 0 else "wb"
    with open(out_tmp_path, mode) as fout:
        while retstart < count:
            retmax = min(page, count - retstart)
            for attempt in range(5):
                try:
                    rr = http_post(f"{EUTILS}/efetch.fcgi", {
                        "db":"nuccore","retmode":"text","rettype":"fasta",
                        "WebEnv":webenv,"query_key":query_key,
                        "retstart":retstart,"retmax":retmax
                    }, stream=True)
                    for chunk in rr.iter_content(chunk_size=1<<15):
                        if chunk: fout.write(chunk)
                    rr.close()
                    break
                except (requests.exceptions.ChunkedEncodingError,
                        requests.exceptions.ConnectionError,
                        ProtocolError, InvalidChunkLength) as e:
                    time.sleep(1.5 * (attempt + 1))
                    if attempt == 3:
                        try:
                            r2 = http_post(f"{EUTILS}/efetch.fcgi",{
                                "db":"nuccore","retmode":"text","rettype":"fasta",
                                "WebEnv":webenv,"query_key":query_key,
                                "retstart":retstart,"retmax":retmax
                            }, stream=False)
                            fout.write(r2.content)
                        except Exception:
                            pass
                    if attempt == 4:
                        raise EutilsHistoryStale(f"efetch stream broken at start={retstart}, size={retmax}: {e}")
            retstart += retmax
            done += retmax
            if tag:
                log(f"{tag} efetch: {done}/{count} (this page {retmax})")
            else:
                log(f"  efetch: {done}/{count} (this page {retmax})")
    return done

def esummary_taxmap(webenv: str, query_key: str, count: int,
                    out_part_path: str, page: int=5000, start: int=0,
                    tag: str = "") -> int:
    done, retstart = start, start
    mode = "a" if start > 0 else "w"
    with open(out_part_path, mode, encoding="utf-8") as out:
        while retstart < count:
            retmax = min(page, count - retstart, 500)  # esummary 单页上限 500
            for attempt in range(5):
                try:
                    r = http_post(f"{EUTILS}/esummary.fcgi", {
                        "db":"nuccore","retmode":"json",
                        "WebEnv":webenv,"query_key":query_key,
                        "retstart":retstart,"retmax":retmax
                    })
                    rs = _json_or_stale(r, what="esummary").get("result", {})
                    break
                except Exception as e:
                    time.sleep(1.2 * (attempt + 1))
                    if attempt == 4:
                        raise EutilsHistoryStale(f"esummary failed at start={retstart}, size={retmax}: {e}")
            for k, v in rs.items():
                if k == "uids": continue
                acc = v.get("caption"); tx = v.get("taxid")
                if acc and tx:
                    out.write(f"{acc}\t{tx}\n"); done += 1
            retstart += retmax
            if tag and ((done - start) % 5000 == 0 or done >= count):
                log(f"{tag} esummary: {done}/{count} (this page {retmax})")
    return done
