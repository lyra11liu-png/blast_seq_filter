# preflight taxonomy check.py

"""
Preflight check: resolve latin names to TaxID + rank + official name,
detect ambiguities/mismatches before downloading genomes.
    - Input CSV columns:
        group, chinese_names, latin_names, synonyms, source, level
    - Output:
        name_resolution.csv with columns:
            latin_input, synonyms_input, best_taxid, best_rank, best_name,
            status(OK/HIGHER_TAXON/AMBIGUOUS/NOT_FOUND), candidates_json, notes
"""

import os, time, json, re, argparse, pathlib, csv
import pandas as pd
import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.getenv("NCBI_API_KEY", "").strip()
HEADERS = {"User-Agent": "pathodb-preflight/1.0 (bioinformatics)"}

RENAMES_HINTS = {
    "Mpox virus": ["Monkeypox virus", "Monkey pox virus", "MPV"],
    "Human betaherpesvirus 5": ["Human cytomegalovirus", "HCMV"],
    "Human alphaherpesvirus 3": ["Varicella-zoster virus", "VZV"],
    "Human gammaherpesvirus 4": ["Epstein-Barr virus", "EBV"]
}

def load_manifest(path: str) -> pd.DataFrame:
    """Automatically detect csv/excel files & import them."""
    p = pathlib.Path(path)
    if not p.exists():
        raise SystemError(f"Manifest not found: {path}")
    if p.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)
    col_map = {}
    for c in list(df.columns):
        lc = c.strip().lower()
        if lc in ("latin_name", "latin names", "latin_names"):
            col_map[c] = "latin_name"
        elif lc in ("chinese_name", "chinese names", "chinese_names"):
            col_map[c] = "chinese_name"
        elif lc in ("synonym", "synonyms"):
            col_map[c] = "synonyms"
    if col_map:
        df = df.rename(columns=col_map)
    if "latin_name" not in df.columns:
        raise SystemError("Input must have a 'latin_name' column.")
    if "synonyms" not in df.columns:
        df["synonyms"] = ""

    df["latin_name"] = df["latin_name"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
    df["synonyms"]   = df["synonyms"].astype(str).str.replace(r"\s+", " ", regex=True).str.strip()
    return df

def esearch_tax(term: str, retmax=5):
    """Search 4 a term in NCBI taxonomy to return a list of TaxIDs."""
    params = {"db": "taxonomy", "retmode": "json", "retmax": str(retmax), "term": term}
    if API_KEY:
        params["api_key"] = API_KEY
    r = requests.get(f"{EUTILS}/esearch.fcgi", params=params, headers=HEADERS, timeout=30)
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", [])

def esummary_tax(taxid: str):
    """Retrieve summary information based on TaxID."""
    params = {"db":"taxonomy", "id": taxid, "retmode":"json"}
    if API_KEY:
        params["api_key"] = API_KEY
    r= requests.get(f"{EUTILS}/esummary.fcgi", params=params, headers=HEADERS, timeout=30)
    r.raise_for_status()
    res = r.json()
    uid = list(res.get("result", {}).get("uids", []))
    if not uid:
        return None
    rec = res["result"][uid[0]]
    return {
        "taxid": taxid,
        "scientificname": rec.get("scientificname"),
        "rank": rec.get("rank"),
        "genus": rec.get("genus"),
        "family": rec.get("family"),
        "division": rec.get("division"),
        "scientificname_authorship": rec.get("scientificname_authorship")
    }

def candidate_taxids_for_name(name: str):
    """Attempt strict scientific names matching; if unsuccessful, use All Names."""
    ids = esearch_tax(f'{name}[Scientific Name]')
    if ids:
        return ids, "scientific"
    ids = esearch_tax(f'{name}[All Names]')
    if ids:
        return ids, "allnames"
    return [], "none"

def resolve_one_name(latin_name: str, synonyms: str, sleep_between=0.2):
    """Parsing individual latin names."""
    latin_name = (latin_name or "").strip()
    syns = [s.strip() for s in re.split(r"[;,/|]", str(synonyms or "")) if s.strip()]
    if latin_name in RENAMES_HINTS:
        syns = list(dict.fromkeys(syns + RENAMES_HINTS[latin_name]))
    for official, alias_list in RENAMES_HINTS.items():
        if latin_name in alias_list:
            syns = list(dict.fromkeys([official] + syns + [latin_name]))
    tried = []
    hits = []
    ids, how = candidate_taxids_for_name(latin_name)
    tried.append((latin_name, how, ids))
    hits += [(latin_name, tid) for tid in ids]
    for s in syns:
        ids, how = candidate_taxids_for_name(s)
        tried.append((s, how, ids))
        hits += [(s, tid) for tid in ids]
    uniq = {}
    for nm, tid in hits:
        if tid not in uniq:
            uniq[tid] = {"by": nm}
    cand = []
    for tid in uniq.keys():
        info = esummary_tax(tid)
        if info:
            cand.append(info)
        time.sleep(sleep_between)
    species_like = [c for c in cand if str(c.get("rank","")).lower() in ("species","subspecies")]
    if len(species_like) == 1:
        best = species_like[0]; status = "OK"
    elif len(species_like) > 1:
        return None, "AMBIGUOUS", cand, tried
    else:
        if cand:
            best = cand[0]; status = "HIGHER_TAXON"
        else:
            return None, "NOT_FOUND", [], tried
    return best, status, cand, tried

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="path to pathogens_CN_WHO_merged.csv")
    ap.add_argument("--out", default="name_resolution.csv")
    ap.add_argument("--sleep", type=float, default=0.2, help="sleep seconds between esummary calls")  # 可调限速
    args = ap.parse_args()

    df = load_manifest(args.manifest)
    total = len(df)
    start = time.time()

    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["latin_input","synonyms_input","best_taxid","best_rank","best_name",
                    "status","candidates_json","note"])
        for i, r in df.iterrows():
            latin = r.get("latin_name", "")
            syns  = r.get("synonyms", "")
            try:
                best, status, cand, tried = resolve_one_name(latin, syns, sleep_between=args.sleep)
                if best:
                    row = [latin, syns, best["taxid"], best["rank"], best["scientificname"],
                           status, json.dumps(cand, ensure_ascii=False), ""]
                else:
                    row = [latin, syns, "", "", "", status,
                           json.dumps(cand, ensure_ascii=False), f"Tried: {tried}"]
            except Exception as e:
                row = [latin, syns, "", "", "", "ERROR", "[]", str(e)]
            w.writerow(row)
            fh.flush() 

            done = i + 1
            elapsed = time.time() - start
            rate = done / elapsed if elapsed > 0 else 0.0
            print(f"[{done}/{total}] {latin} -> {row[5]}  | {rate:.2f} rec/s", flush=True)

    print(f"wrote {args.out}")

if __name__ == "__main__":
    main()
