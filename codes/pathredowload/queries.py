# -*- coding: utf-8 -*-
def build_query(taxid: str, category: str, *, slen_min: int,
                only_complete: bool, refseq_only: bool, viruses_include_genbank: bool) -> str:
    term = f"txid{taxid}[Subtree]"
    if slen_min and slen_min > 0:
        term += f" AND {slen_min}:100000000[SLEN]"
    term += " AND biomol_genomic[PROP]"

    if only_complete:
        term += " AND (complete[Title] OR complete[All Fields])"

    if category.lower() == "viruses":
        if viruses_include_genbank:
            term += " AND (srcdb_refseq[PROP] OR srcdb_genbank[PROP])"
        else:
            term += " AND srcdb_refseq[PROP]"
    else:
        if refseq_only:
            term += " AND srcdb_refseq[PROP]"
        else:
            term += " AND (srcdb_refseq[PROP] OR srcdb_genbank[PROP])"

    return term
