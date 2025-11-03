#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download_build_blastdb_v2.py

从 CSV 读取 TaxID（支持多ID用 | , ; 空白分隔），用 E-utilities 下载 nuccore 序列，
按病原体大类分目录（viruses/bacteria/fungi/archaea/other），合并并 makeblastdb 建库。
改进点：
- efetch 用 POST + 自适应分批，避免 414，自动重试
- 仅写 .part 临时文件，成功后 rename，杜绝 0KB
- 更稳健的分组（lineage 优先；fallback：名字里含 "virus" 也判为 viruses）
- 支持 reclassify（重分组迁移已有 fasta）与 fix-empty（清理 0KB/0条）
- 新增 verify/tolerance/redo-incomplete：校验“现有条数 vs 期望条数”，不完整可自动重下
- 生成 incomplete_or_missing.tsv 清单，便于核对遗漏

export NCBI_API_KEY=1bd9501dcc57f92d38f6eb04e43c1bc21c09
python download_build_db.py \
  --csv /home/liuyuxin/projects/blast_seq_filter/data/name_resolution.csv \
  --taxid-col best_taxid --name-col best_name \
  --outdir /data1/liuyuxin/blast_seq_filter/original_data_lib \
  --threads 8 \
  --fix-empty --reclassify \
  --verify --tolerance 0.02 --redo-incomplete \
  --merge --build-db
"""
import os, re, sys, json, time, pathlib, argparse, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

EUTILS   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY  = (os.getenv("NCBI_API_KEY") or "").strip()

def now(): return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
def log(msg): print(f"[{now()}] {msg}", flush=True)

# ---------------- HTTP ----------------
def _http_get(url, params, retry=5, sleep=0.5):
    if API_KEY: params = {**params, "api_key": API_KEY}
    qs = urllib.parse.urlencode(params, safe="()[],:|")
    full = f"{url}?{qs}"
    for i in range(retry):
        try:
            with urllib.request.urlopen(full, timeout=120) as r:
                return r.read()
        except Exception as e:
            time.sleep(sleep*(2**i))
            log(f"HTTP 重试 {i+1}/{retry}：{e}")
    raise RuntimeError(f"GET 失败：{full}")

def _http_post(url, params, retry=5, sleep=0.5):
    if API_KEY: params = {**params, "api_key": API_KEY}
    data = urllib.parse.urlencode(params, safe="()[],:|").encode("utf-8")
    for i in range(retry):
        try:
            req = urllib.request.Request(url, data=data, method="POST")
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read()
        except Exception as e:
            time.sleep(sleep*(2**i))
            log(f"HTTP(POST) 重试 {i+1}/{retry}：{e}")
    raise RuntimeError(f"POST 失败：{url}")

def esearch_count(term):
    """只取 nuccore 计数（不拉ID列表）"""
    raw = _http_get(f"{EUTILS}/esearch.fcgi",
                    {"db":"nuccore","term":term,"retmode":"json","retmax":0})
    j = json.loads(raw.decode("utf-8"))
    return int(j.get("esearchresult", {}).get("count", "0"))

# -------------- Taxonomy 分组 ----------
def tax_esummary(taxid: str):
    raw = _http_get(f"{EUTILS}/esummary.fcgi",
                    {"db":"taxonomy","id":taxid,"retmode":"json"})
    j = json.loads(raw.decode("utf-8"))
    d = j.get("result", {}).get(str(taxid), {})
    return {
        "taxid": str(d.get("taxid","")),
        "scientificname": d.get("scientificname",""),
        "rank": d.get("rank",""),
        "division": (d.get("division") or ""),
        "lineage": (d.get("lineage") or "")
    }

def map_group(meta: dict, name_hint: str="") -> str:
    """更稳健的分组：优先看 lineage，次看 division，最后看名字"""
    lin = meta.get("lineage","").lower()
    div = meta.get("division","").lower()
    hint = (name_hint or meta.get("scientificname","")).lower()

    def has(x): return x in lin or x in div
    if has("viruses") or "virus" in hint:
        return "viruses"
    if has("bacteria"):
        return "bacteria"
    if has("fungi") or "fungus" in lin:
        return "fungi"
    if has("archaea"):
        return "archaea"
    return "other"

# -------------- Entrez: search/fetch ----
def esearch_nuccore(term, retmax=10000):
    ids, start = [], 0
    while True:
        raw = _http_get(f"{EUTILS}/esearch.fcgi",
                        {"db":"nuccore","term":term,"retmode":"json",
                         "retstart": start, "retmax": retmax})
        er = json.loads(raw.decode("utf-8")).get("esearchresult", {})
        chunk = er.get("idlist", [])
        ids.extend(chunk)
        count = int(er.get("count","0"))
        start += retmax
        if start >= count or not chunk: break
    return ids

def efetch_fasta_post(ids, out_tmp_path):
    """POST + 自适应分批写 .part，成功后由调用方 rename"""
    B, done, total = 400, 0, len(ids)
    with open(out_tmp_path, "wb") as fout:
        i = 0
        while i < total:
            j = min(i + B, total)
            batch = ",".join(ids[i:j])
            try:
                raw = _http_post(
                    f"{EUTILS}/efetch.fcgi",
                    {"db": "nuccore", "id": batch, "rettype": "fasta", "retmode": "text"}
                )
                n_batch = j - i                 # 先算本批数量
                fout.write(raw)
                done += n_batch
                i = j
                log(f"    efetch 批完成：{done}/{total}（本批{n_batch}）")
                if B < 600:
                    B = min(600, B + 50)
            except Exception as e:
                if B <= 25:
                    raise
                oldB = B
                B = max(25, B // 2)
                log(f"    efetch 批失败（{e}），batch {oldB}->{B}，重试 …")
    return done

# -------------- IO/CSV/工具 -------------
def read_csv_rows(path):
    import csv
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f): yield r

def parse_taxids(s: str):
    if not s: return []
    return [p for p in re.split(r"[|\s,;]+", s.strip()) if p.isdigit()]

def sanitize(s: str):
    s = (s or "").strip()
    s = re.sub(r"\s+","_", s)
    return re.sub(r"[^A-Za-z0-9_.-]","", s)[:80] or "unnamed"

def count_fasta(path):
    try:
        with open(path,"r",encoding="utf-8",errors="ignore") as f:
            return sum(1 for ln in f if ln.startswith(">"))
    except: return 0

# -------------- 核心流程：处理一个 taxid ----
def handle_one_taxid(taxid: str, name_hint: str, out_base: pathlib.Path,
                     force=False, verify=False, tolerance=0.02, redo_incomplete=False):
    t0 = time.time()
    meta = tax_esummary(taxid)
    group = map_group(meta, name_hint)
    gdir = out_base / group
    gdir.mkdir(parents=True, exist_ok=True)

    label = sanitize(name_hint or meta.get("scientificname") or f"txid{taxid}")
    fa_path  = gdir / f"{taxid}_{label}.fasta"
    tmp_path = gdir / f".{taxid}_{label}.fasta.part"

    term = f"txid{taxid}[Organism:exp]"

    # --- 已存在：可选校验并决定是否重下 ---
    if fa_path.exists() and not force:
        nrec = count_fasta(fa_path)
        expected = esearch_count(term) if verify else 0
        incomplete = False
        if verify and expected > 0 and nrec < int(expected * (1 - tolerance)):
            incomplete = True
            log(f"[校验] txid {taxid} {fa_path.name} 现有={nrec}，应≈{expected}（<{(1-tolerance)*100:.0f}%），标记不完整")
            if redo_incomplete:
                try: os.remove(fa_path)
                except: pass
            else:
                return {"taxid":taxid,"group":group,"file":str(fa_path),
                        "count":nrec,"seconds":0.0,"skipped":True,
                        "sci_name":meta.get("scientificname",""),
                        "rank":meta.get("rank",""),
                        "expected":expected, "incomplete":True}
        else:
            log(f"[跳过] txid {taxid} 已存在：{fa_path.name}（{nrec} 条）")
            return {"taxid":taxid,"group":group,"file":str(fa_path),
                    "count":nrec,"seconds":0.0,"skipped":True,
                    "sci_name":meta.get("scientificname",""),
                    "rank":meta.get("rank",""),
                    "expected":expected, "incomplete":False}

    # --- 新下载或需重下 ---
    log(f"开始搜索：txid {taxid} | term='{term}'")
    ids = esearch_nuccore(term)
    log(f"  找到 {len(ids)} 条 nuccore 记录，开始 efetch->FASTA …")
    got = 0
    if ids:
        if tmp_path.exists():
            try: tmp_path.unlink()
            except: pass
        got = efetch_fasta_post(ids, str(tmp_path))
        os.replace(str(tmp_path), str(fa_path))
    sec = time.time()-t0
    log(f"完成：txid {taxid} -> {fa_path.name} | 条数={got} | 用时={sec:.1f}s | 分组={group}")
    return {"taxid":taxid,"group":group,"file":str(fa_path),
            "count":got,"seconds":sec,"skipped":False,
            "sci_name":meta.get("scientificname",""),
            "rank":meta.get("rank",""),
            "expected":len(ids), "incomplete":False}

# -------------- 合并 & 建库 ---------------
def merge_group_fastas(group_dir: pathlib.Path, merged_path: pathlib.Path):
    seen, n_in, n_out = set(), 0, 0
    with open(merged_path,"w",encoding="utf-8") as fout:
        for fa in sorted(group_dir.glob("*.fasta")):
            if fa.name.endswith(".part"): continue
            with open(fa,"r",encoding="utf-8",errors="ignore") as f:
                write=False
                for ln in f:
                    if ln.startswith(">"):
                        n_in += 1
                        acc = ln[1:].strip().split()[0]
                        write = acc not in seen
                        if write: seen.add(acc); n_out += 1; fout.write(ln)
                    else:
                        if write: fout.write(ln)
    return n_in, n_out

def ensure_makeblastdb(bin_hint):
    import shutil
    cand = bin_hint or os.getenv("MAKEBLASTDB") or shutil.which("makeblastdb")
    if not cand:
        raise FileNotFoundError("找不到 makeblastdb（可用 --makeblastdb 指定或设环境变量）。")
    return cand

def build_blast_db(makeblastdb, fasta_in, out_prefix, title=None, dbtype="nucl"):
    title = title or pathlib.Path(out_prefix).name
    cmd = f'"{makeblastdb}" -dbtype {dbtype} -in "{fasta_in}" -parse_seqids -title "{title}" -out "{out_prefix}"'
    log(f"建库命令：{cmd}")
    code = os.system(cmd)
    if code != 0: raise RuntimeError(f"makeblastdb 失败（退出码 {code}）")

# -------------- 维护工具：重分组 / 修复 0KB ----
def reclassify_all(out_base: pathlib.Path):
    """扫描所有 *.fasta，按 taxid 重新判组并搬移到正确目录"""
    moved = 0
    for fa in out_base.rglob("*.fasta"):
        if fa.name.endswith(".part"): continue
        m = re.match(r"^(\d+)_", fa.name)
        if not m: continue
        taxid = m.group(1)
        meta = tax_esummary(taxid)
        new_group = map_group(meta, fa.name)
        if fa.parent.name != new_group:
            tgt_dir = out_base / new_group
            tgt_dir.mkdir(parents=True, exist_ok=True)
            new_path = tgt_dir / fa.name
            log(f"[重分组] {fa}  ->  {new_path}")
            os.replace(str(fa), str(new_path))
            moved += 1
    log(f"[重分组完成] 共移动 {moved} 个文件。")

def fix_empty_fastas(out_base: pathlib.Path):
    """删除 0KB 或 0 条记录的 fasta"""
    removed = 0
    for fa in out_base.rglob("*.fasta"):
        if fa.stat().st_size == 0 or count_fasta(fa) == 0:
            log(f"[清理空文件] {fa}")
            try: fa.unlink(); removed += 1
            except: pass
    log(f"[清理完成] 共移除 {removed} 个空文件。")

# -------------- 主程序 --------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="包含 TaxID 的 CSV")
    ap.add_argument("--taxid-col", default="best_taxid")
    ap.add_argument("--name-col", default="best_name")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--build-db", action="store_true")
    ap.add_argument("--makeblastdb", default=None)
    ap.add_argument("--reclassify", action="store_true", help="对现有文件重算分组并迁移")
    ap.add_argument("--fix-empty", action="store_true", help="删除 0KB/0条记录的空文件")
    # 新增：校验相关
    ap.add_argument("--verify", action="store_true", help="校验已存在 fasta 的条目数 vs ESearch 计数")
    ap.add_argument("--tolerance", type=float, default=0.02, help="校验误差容忍度（默认 0.02=2%）")
    ap.add_argument("--redo-incomplete", action="store_true", help="对校验为不完整的条目自动重下")
    args = ap.parse_args()

    out_base = pathlib.Path(args.outdir).resolve()
    out_base.mkdir(parents=True, exist_ok=True)

    # 维护动作
    if args.fix_empty:    fix_empty_fastas(out_base)
    if args.reclassify:   reclassify_all(out_base)

    # 收集任务
    tasks = []
    import csv
    with open(args.csv, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            taxids = parse_taxids(row.get(args.taxid_col, ""))
            if not taxids: continue
            label = row.get(args.name_col, "") or row.get("latin_input","") or ""
            for tx in taxids:
                tasks.append((tx, label))

    log(f"共 {len(tasks)} 个 TaxID 任务。API_KEY={'SET' if API_KEY else 'NOT SET'}")
    t0 = time.time()
    results = []
    if tasks:
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = [ex.submit(handle_one_taxid, t, nm, out_base,
                              args.force, args.verify, args.tolerance, args.redo_incomplete)
                    for t,nm in tasks]
            for fu in as_completed(futs):
                results.append(fu.result())

    # 统计
    stat = out_base / "download_stats.tsv"
    with open(stat, "w", encoding="utf-8") as f:
        f.write("taxid\tgroup\tsci_name\trank\tfile\tcount\tseconds\tskipped\texpected\tincomplete\n")
        for r in sorted(results, key=lambda x:(x["group"], x["taxid"])):
            f.write("\t".join([
                str(r["taxid"]), r["group"], r.get("sci_name",""), r.get("rank",""),
                r["file"], str(r["count"]), f'{r["seconds"]:.1f}', str(r["skipped"]),
                str(r.get("expected",0)), str(r.get("incomplete",False))
            ]) + "\n")
    log(f"已写出统计：{stat}")

    # 未完成/缺失清单
    miss_path = out_base / "incomplete_or_missing.tsv"
    with open(miss_path, "w", encoding="utf-8") as f:
        f.write("taxid\tgroup\tsci_name\trank\tfile\texpected\tcount\n")
        for r in sorted(results, key=lambda x:(x["group"], x["taxid"])):
            expd = int(r.get("expected",0))
            cnt  = int(r.get("count",0) or 0)
            inc  = bool(r.get("incomplete", False))
            if inc or (expd>0 and cnt==0 and not r.get("skipped",False)):
                f.write("\t".join([
                    str(r["taxid"]), r["group"], r.get("sci_name",""), r.get("rank",""),
                    r.get("file",""), str(expd), str(cnt)
                ]) + "\n")
    log(f"未完成/缺失清单：{miss_path}（为空表示都正常）")

    # 合并 & 建库
    if args.merge or args.build_db:
        mk = ensure_makeblastdb(args.makeblastdb) if args.build_db else None
        for g in ["viruses","bacteria","fungi","archaea","other"]:
            gdir = out_base / g
            if not gdir.exists(): continue
            merged = gdir / f"{g}_all.fasta"
            log(f"[合并] {g} -> {merged.name}")
            n_in, n_out = merge_group_fastas(gdir, merged)
            log(f"  合并完成：入={n_in}，去重后={n_out}")
            if args.build_db:
                prefix = gdir / f"{g}_blastdb"
                build_blast_db(mk, str(merged), str(prefix), title=f"{g}_nucl")
    log(f"全部完成，用时 {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main()