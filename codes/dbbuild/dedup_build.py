#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, sys, gzip, time, shutil, sqlite3, hashlib, subprocess
from pathlib import Path
from typing import Iterable, Tuple, List, Optional, Dict, Union
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import threading
from tqdm import tqdm
import signal
STOP = False

def _sig_handler(sig, frame):
    global STOP
    if not STOP:
        STOP = True
        sys.stderr.write("\n[INTERRUPT] 捕获到中断信号：正在安全停止（会终止子进程），已完成的物种不会重跑。\n")
        for pid in list(_CHILD_PIDS):
            try: os.killpg(pid, signal.SIGTERM)
            except Exception: pass
            
signal.signal(signal.SIGINT, _sig_handler)
signal.signal(signal.SIGTERM, _sig_handler)

def spawn_log_watcher(log_file: Path, desc: str):
    stop_evt = threading.Event()
    bar = tqdm(total=100, desc=desc, unit="%", leave=False, position=1)
    def _watch():
        last = 0
        while not stop_evt.is_set():
            try:
                with open(log_file, "r", errors="ignore") as f:
                    f.seek(last)
                    for line in f:
                        # mmseqs 日志里常见“... 23% ...”这样的进度
                        m = re.search(r"(\d{1,3})\s*%", line)
                        if m:
                            p = max(0, min(100, int(m.group(1))))
                            bar.n = p; bar.refresh()
                    last = f.tell()
            except FileNotFoundError:
                pass
            time.sleep(0.5)
        bar.close()
    t = threading.Thread(target=_watch, daemon=True); t.start()
    return stop_evt

_DB_LOCK = threading.Lock()

FA_EXT = (".fa", ".fna", ".fasta", ".fas")
NOW = lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# -------------------- 基础工具 --------------------
# 允许的碱基集合（建库最稳妥只保留 A/C/G/T/N）
_ALLOWED = set("ACGTN")

def sanitize_and_flag(seq: str):
    """
    返回: (clean_seq, bad_pos, bad_chars)
      - clean_seq: 清洗后的序列（U->T，其它->N）
      - bad_pos:   原序列里非法字符的位置(0-based)
      - bad_chars: 这些非法字符的去重集合（列表）
    """
    raw = seq.upper().replace("U", "T")
    bad_pos = [i for i, ch in enumerate(raw) if ch not in _ALLOWED]
    if bad_pos:
        clean = "".join((ch if ch in _ALLOWED else "N") for ch in raw)
        bad_chars = sorted(set(raw[i] for i in bad_pos))
    else:
        clean = raw
        bad_chars = []
    return clean, bad_pos, bad_chars

def _sanitize_dna(seq: str) -> str:
    s = seq.upper().replace("U", "T")
    import re
    return re.sub(r"[^ACGTN]", "N", s)

def write_fasta_clean(fp: str, records):
    with open(fp, "w") as w:
        for h, s in records:
            cs = _sanitize_dna(s)
            w.write(">"+h+"\n")
            for i in range(0, len(cs), 60):
                w.write(cs[i:i+60]+"\n")
                
# ==== 短ID与精简标题工具 ====
MAX_STITLE_LEN = int(os.getenv("MAX_STITLE_LEN", "240"))

def _short_id(cat: str, idx: int, seq: str) -> str:
    """≤50字符的稳定短ID：<cat>_<7位序号>_<sha10>"""
    h10 = hashlib.sha1(seq.encode("utf-8")).hexdigest()[:10]
    base = re.sub(r"[^A-Za-z0-9_]", "_", cat)[:12]
    return f"{base}_{idx:07d}_{h10}"  # e.g. bacteria_0001234_ab12cd34ef

_ACC_RE = re.compile(r'\b([A-Z]{1,3}\d{2,8}\.\d+|[A-Z]{2}_\d+\.\d+)\b')

def _parse_safe_id(safe_id: str) -> Tuple[str, str]:
    """
    你的 merge 阶段把 header 写成：
        safe_id = "<unit>|<file>|<safe(orig)>"
    这里取出 unit 与 file（如果没写也容错）。
    """
    parts = safe_id.split("|", 2)
    unit = parts[0] if parts else ""
    srcf = parts[1] if len(parts) > 1 else ""
    return unit, srcf

def _compact_stitle(unit: str, srcf: str, orig: str) -> str:
    """
    生成精简 stitle："<物种> | file=<file> | acc=<ACC> | <注释>"
    - 去掉重复的物种名前缀
    - 压缩空白并截断到 MAX_STITLE_LEN
    """
    species = unit.replace("_", " ").strip()
    acc = ""
    m = _ACC_RE.search(orig)
    if m: acc = m.group(1)

    # 注释 = 去掉开头 accession、去掉与物种重复的冗余前缀
    anno = re.sub(r'^[^\s]+\s*', '', orig).strip()
    if anno.lower().startswith(species.lower()):
        anno = anno[len(species):].lstrip(" :,-")

    parts = [species]
    if srcf: parts.append(f"file={Path(srcf).name}")
    if acc:  parts.append(f"acc={acc}")
    if anno: parts.append(re.sub(r'\s+', ' ', anno))
    stitle = " | ".join(parts).strip()
    if len(stitle) > MAX_STITLE_LEN:
        stitle = stitle[:MAX_STITLE_LEN-3] + "..."
    return stitle

def need(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

_CHILD_PIDS = set()

def run_popen(cmd: List[str], log_fp: Optional[str] = None, cwd: Optional[str] = None) -> int:
    """启动外部程序并将其放入新会话；若收到 STOP，则整组 SIGTERM 结束。"""
    with open(log_fp, "a") if log_fp else open(os.devnull, "w") as logf:
        proc = subprocess.Popen(
            cmd, cwd=cwd, stdout=logf, stderr=logf,
            start_new_session=True  # 关键：让该进程成为新进程组
        )
        _CHILD_PIDS.add(proc.pid)
        try:
            while True:
                try:
                    return proc.wait(timeout=1.0)
                except subprocess.TimeoutExpired:
                    if STOP:
                        try:
                            os.killpg(proc.pid, signal.SIGTERM)  # 结束整个进程组
                        except Exception:
                            pass
        finally:
            _CHILD_PIDS.discard(proc.pid)

def run(cmd: List[str]) -> int:
    sys.stderr.write("[CMD] " + " ".join(cmd) + "\n")
    return subprocess.run(cmd).returncode

def ensure_dir(p: Union[str, Path]):
    Path(p).mkdir(parents=True, exist_ok=True)

def safe(s: str) -> str:
    # 仅允许字母/数字/下划线/点/连字符；其他字符一律替换为下划线
    allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-")
    return "".join((c if c in allowed else "_") for c in s)

def is_complete_header(h: str) -> bool:
    hl = h.lower()
    good = ("complete", "chromosome", "genome", "segment", "scaffold", "contig")
    bad  = ("partial", "fragment")
    return any(k in hl for k in good) and not any(k in hl for k in bad)

def pick_rep(headers: List[str], seqs: List[str]) -> int:
    idxs = list(range(len(headers)))
    comp = [i for i in idxs if is_complete_header(headers[i])]
    if comp: idxs = comp
    idxs.sort(key=lambda i: len(seqs[i]), reverse=True)
    return idxs[0]

def sha1(s: str) -> str:
    import hashlib
    return hashlib.sha1(s.encode("utf-8")).hexdigest()

# -------------------- FASTA IO --------------------
def read_fasta_lazy(fp: str) -> Iterable[Tuple[str, str]]:
    op = gzip.open if fp.endswith(".gz") else open
    with op(fp, "rt", encoding="utf-8", errors="ignore") as f:
        h, seq = None, []
        for line in f:
            if line.startswith(">"):
                if h is not None:
                    yield (h, "".join(seq).replace(" ", "").replace("\r", "").replace("\n", ""))
                h = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if h is not None:
            yield (h, "".join(seq).replace(" ", "").replace("\r", "").replace("\n", ""))

def write_fasta(fp: str, records: Iterable[Tuple[str, str]]):
    with open(fp, "w") as w:
        for h, s in records:
            w.write(">"+h+"\n")
            for i in range(0, len(s), 60):
                w.write(s[i:i+60]+"\n")

# -------------------- 目录发现（智能“物种单元”） --------------------
def has_subdir_with_fasta(cat_dir: Path) -> bool:
    for p in cat_dir.iterdir():
        if p.is_dir():
            if any(q.is_file() and (q.suffix in FA_EXT or str(q).endswith(".gz")) for q in p.rglob("*")):
                return True
    return False

def discover_species_units(cat_dir: Path) -> List[Path]:
    """
    若存在子目录含fasta：以“子目录”为物种单元；
    否则：以“根目录下每个 fasta 文件”为物种单元（避免把整个类别当成一个巨物种）。
    """
    units: List[Path] = []
    if has_subdir_with_fasta(cat_dir):
        seen = set()
        for fa in cat_dir.rglob("*"):
            if fa.is_file() and (fa.suffix in FA_EXT or str(fa).endswith(".gz")):
                seen.add(fa.parent)
        units = sorted(seen)
    else:
        units = sorted([fa for fa in cat_dir.iterdir()
                        if fa.is_file() and (fa.suffix in FA_EXT or str(fa).endswith(".gz"))])
    return units

def collect_fastas(unit: Path) -> List[Path]:
    """unit 可以是目录（返回其下所有 fasta），也可以是单个 fasta 文件（返回 [unit]）。"""
    if unit.is_file():
        return [unit]
    return sorted([p for p in unit.iterdir()
                   if p.is_file() and (p.suffix in FA_EXT or str(p).endswith(".gz"))])

def unit_name(unit: Path) -> str:
    return unit.stem if unit.is_file() else unit.name

# -------------------- SQLite 追踪 --------------------
def db_init(dbp: Path):
    # 关键：加 timeout，自动提交，busy_timeout；WAL 失败回退到 DELETE
    con = sqlite3.connect(
        dbp.as_posix(),
        timeout=600.0,                 # 等锁最长 600s
        check_same_thread=False,
        isolation_level=None           # autocommit，减少持有锁时间
    )
    con.execute("PRAGMA busy_timeout=600000;")  # 600s
    try:
        con.execute("PRAGMA journal_mode=WAL;") # 写放 WAL
    except sqlite3.OperationalError:
        con.execute("PRAGMA journal_mode=DELETE;")  # 回退
    con.execute("PRAGMA synchronous=NORMAL;")   # 速度/安全折中
    cur = con.cursor()
    cur.execute("""CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT)""")
    cur.execute("""CREATE TABLE IF NOT EXISTS species(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        category TEXT, species_dir TEXT,
        status TEXT, method TEXT,
        n_total INTEGER, n_after_exact INTEGER, n_kept INTEGER, n_dropped INTEGER,
        started_at TEXT, finished_at TEXT,
        UNIQUE(category, species_dir)
    )""")
    cur.execute("""CREATE TABLE IF NOT EXISTS seq_log(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        species_id INTEGER,
        header TEXT, seqlen INTEGER, hash TEXT,
        action TEXT, reason TEXT, ref_header TEXT
    )""")
    return con

def species_begin(con, category: str, sp_dir: str, method: str) -> int:
    with _DB_LOCK:
        cur = con.cursor()
        cur.execute("""INSERT INTO species(category, species_dir, status, method, started_at)
                       VALUES(?,?,?,?,?)
                       ON CONFLICT(category, species_dir)
                       DO UPDATE SET status=excluded.status, method=excluded.method, started_at=excluded.started_at""",
                    (category, sp_dir, "working", method, NOW()))
        con.commit()
        cur = con.cursor()
        cur.execute("SELECT id FROM species WHERE category=? AND species_dir=?", (category, sp_dir))
        return cur.fetchone()[0]

def species_done(con, sid: int, n_total: int, n_after_exact: int, n_kept: int, n_dropped: int):
    with _DB_LOCK:
        cur = con.cursor()
        cur.execute("""UPDATE species SET status=?, n_total=?, n_after_exact=?, n_kept=?, n_dropped=?, finished_at=?
                       WHERE id=?""",
                    ("done", n_total, n_after_exact, n_kept, n_dropped, NOW(), sid))
        con.commit()

def species_status(con, category: str, sp_dir: str) -> Optional[str]:
    with _DB_LOCK:
        cur = con.cursor()
        cur.execute("SELECT status FROM species WHERE category=? AND species_dir=?", (category, sp_dir))
        row = cur.fetchone()
        return row[0] if row else None

def log_seq(con, sid: int, header: str, seqlen: int, hsh: str, action: str, reason: str, ref_header: Optional[str]):
    with _DB_LOCK:
        cur = con.cursor()
        cur.execute("""INSERT INTO seq_log(species_id, header, seqlen, hash, action, reason, ref_header)
                       VALUES(?,?,?,?,?,?,?)""",
                    (sid, header, seqlen, hsh, action, reason, ref_header))
        con.commit()

def export_tsv(con, outdir: Path):
    with _DB_LOCK:
        kept = outdir/"kept.tsv"; dropped = outdir/"dropped.tsv"
        for fp, act in [(kept, "kept"), (dropped, "dropped")]:
            with open(fp, "w") as w:
                w.write("category\tspecies_dir\theader\tseqlen\treason\tref_header\n")
                cur = con.cursor()
                cur.execute("""SELECT s.category, s.species_dir, q.header, q.seqlen, q.reason, q.ref_header
                               FROM seq_log q JOIN species s ON q.species_id=s.id
                               WHERE q.action=?""", (act,))
                for row in cur.fetchall():
                    w.write("\t".join([str(x) if x is not None else "" for x in row])+"\n")

# -------------------- Phase 1: 精确去重 --------------------
def phase1_exact(records: List[Tuple[str,str]], sid: int, con, pbar_desc: str) -> Tuple[List[Tuple[str,str]], int]:
    headers = [h for h,_ in records]
    seqs    = [s for _,s in records]
    # hash 分桶
    buckets: Dict[str, List[int]] = {}
    bar = tqdm(range(len(seqs)), desc=pbar_desc, unit="seq", leave=False)
    for i in bar:
        hsh = sha1(seqs[i])
        buckets.setdefault(hsh, []).append(i)

    kept_idx = []
    dropped = 0
    for idxs in buckets.values():
        if len(idxs)==1:
            i = idxs[0]; kept_idx.append(i)
            log_seq(con, sid, headers[i], len(seqs[i]), sha1(seqs[i]), "kept", "unique", None)
        else:
            rep = pick_rep([headers[i] for i in idxs], [seqs[i] for i in idxs])
            rep_idx = idxs[rep]
            kept_idx.append(rep_idx)
            log_seq(con, sid, headers[rep_idx], len(seqs[rep_idx]), sha1(seqs[rep_idx]), "kept", "exact_rep", None)
            for i in idxs:
                if i == rep_idx: continue
                dropped += 1
                log_seq(con, sid, headers[i], len(seqs[i]), sha1(seqs[i]), "dropped", "exact_dup", headers[rep_idx])
    kept_idx = sorted(set(kept_idx))
    kept = [(headers[i], seqs[i]) for i in kept_idx]
    return kept, dropped

# -------------------- Phase 2: 近重复（优先 MMseqs2） --------------------
def mmseqs_near(in_fa: Path, out_fa: Path, tmpdir: Path,
                pid: float, qcov: float, sid: int, con, log_file: Path, threads: int, species_tag: str):
    if STOP:
        return
    ensure_dir(tmpdir)
    mm_in = tmpdir / "phase1.mmseqs.fa"
    recs = list(read_fasta_lazy(in_fa.as_posix()))
    write_fasta_clean(mm_in.as_posix(), recs)
    cmd = ["mmseqs","easy-cluster", 
           mm_in.as_posix(), (tmpdir/"cluster").as_posix(), tmpdir.as_posix(),
           "--dbtype","2", "--min-seq-id", str(pid), "-c", str(qcov),
           "--cov-mode","1","--threads", str(max(1, threads))]
    watcher = spawn_log_watcher(log_file, f"[{species_tag}] mmseqs")
    try:
        rc = run_popen(cmd, log_fp=log_file.as_posix())
    finally:
        watcher.set()
    if rc != 0:
        raise RuntimeError("mmseqs easy-cluster 失败，详见日志："+log_file.as_posix())

    # 找到 *_cluster.tsv
    clu_tsv = None
    for f in os.listdir(tmpdir):
        if f.endswith("_cluster.tsv"):
            clu_tsv = tmpdir / f
            break

    # 读入序列，建两个索引：完整ID 与 空格前token
    H, S = [], []
    for h, s in read_fasta_lazy(in_fa.as_posix()):
        H.append(h); S.append(s)
    idx_full = {H[i]: i for i in range(len(H))}
    idx_tok  = {H[i].split()[0]: i for i in range(len(H))}

    # 没有簇文件/空文件：直接原样输出
    if (clu_tsv is None) or (os.path.getsize(clu_tsv) == 0):
        write_fasta(out_fa.as_posix(), [(h, s) for h, s in zip(H, S)])
        for i in range(len(H)):
            log_seq(con, sid, H[i], len(S[i]), sha1(S[i]), "kept", "near_keep", None)
        return

    from collections import defaultdict
    clusters = defaultdict(set)
    with open(clu_tsv) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep, mem = parts[0], parts[1]
            clusters[rep].add(rep)
            clusters[rep].add(mem)

    chosen = []
    seen = set()
    for rep, members in clusters.items():
        ids = []
        for m in members:
            if m in idx_full:
                ids.append(idx_full[m])
            else:
                t = m.split()[0]
                if t in idx_tok:
                    ids.append(idx_tok[t])
        if not ids:
            # 这个簇里的ID都对不上，跳过，避免越界
            continue

        pick = pick_rep([H[i] for i in ids], [S[i] for i in ids])
        rep_i = ids[pick]
        chosen.append(rep_i)
        seen.update(ids)
        log_seq(con, sid, H[rep_i], len(S[rep_i]), sha1(S[rep_i]), "kept", "near_rep", None)
        for i in ids:
            if i == rep_i:
                continue
            reason = "near_dup" if len(S[i]) <= len(S[rep_i]) else "near_dup_longer_rep"
            log_seq(con, sid, H[i], len(S[i]), sha1(S[i]), "dropped", reason, H[rep_i])
            
    for i in range(len(H)):
        if i not in seen and i not in chosen:
            chosen.append(i)
            log_seq(con, sid, H[i], len(S[i]), sha1(S[i]), "kept", "near_keep", None)

    chosen = sorted(set(chosen))
    write_fasta(out_fa.as_posix(), [(H[i], S[i]) for i in chosen])

# === 替换原来的 blast_near(...) 整个函数 ===
def blast_near(in_fa: Path, out_fa: Path, tmpdir: Path,
               pid: float, qcov: float, threads: int,
               sid: int, con, pbar_desc: str):
    """
    用 batched BLAST 做“近重复”。
    临时 subject（kept.fa/kept_db）全部写入 tmpdir（来自 --tmp-root），
    若 makeblastdb 失败自动回退到 -subject FASTA。
    """
    ensure_dir(tmpdir)

    # 读取并按“更完整/更长优先”排序
    recs = list(read_fasta_lazy(in_fa.as_posix()))
    recs.sort(key=lambda x: (not is_complete_header(x[0]), -len(x[1])))

    # 统一的物种本地临时前缀（在 tmpdir 里）
    tmp_base = (tmpdir / "blast").as_posix()
    kept_file = tmp_base + ".kept.fa"
    dbpref    = tmp_base + ".keptdb"

    # —— 可调参数（环境变量覆盖）——
    BATCH = int(os.getenv("BLAST_BATCH", "256"))
    WORD  = int(os.getenv("BLAST_WORD",  "48"))
    UNGAP = os.getenv("BLAST_UNGAPPED", "1").lower() in ("1","true","yes")
    USE_DB_GLOBAL = os.getenv("BLAST_SUBJECT_DB", "0").lower() in ("1","true","yes")
    DB_REFRESH = int(os.getenv("BLAST_DB_REFRESH", "5000"))
    SUBJECT_REFRESH = int(os.getenv("BLAST_SUBJECT_REFRESH", "400"))
    DB_MIN_SEQS = int(os.getenv("BLAST_DB_MIN_SEQS", "10000"))
    DB_MIN_BYTES = int(os.getenv("BLAST_DB_MIN_BYTES", str(400*1024*1024)))
    CLEAN_TMP = os.getenv("BLAST_CLEAN_TMP", "1").lower() in ("1","true","yes")

    def _write_kept():
        write_fasta(kept_file, kept)

    def _make_db() -> bool:
        # 临时 DB 不用 -parse_seqids，鲁棒性更好
        try:
            for ext in (".nhr",".nin",".nsq"):
                p = dbpref + ext
                if os.path.exists(p): os.remove(p)
            proc = subprocess.run(
                ["makeblastdb","-in", kept_file, "-dbtype","nucl", "-out", dbpref],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            if proc.returncode != 0:
                sys.stderr.write(f"[WARN] makeblastdb(kept_db) 失败，回退 FASTA subject。"
                                 f" stderr: {proc.stderr.splitlines()[-1][:300]}\n")
                return False
            return True
        except Exception as e:
            sys.stderr.write(f"[WARN] makeblastdb(kept_db) 异常：{e}，回退 FASTA subject。\n")
            return False

    def _use_db_now() -> bool:
        try:
            kb = os.path.getsize(kept_file)
        except FileNotFoundError:
            kb = 0
        return USE_DB_GLOBAL and (len(kept) >= DB_MIN_SEQS or kb >= DB_MIN_BYTES)

    def _blast_batch(batch):
        if not batch:
            return {}
        qfa = tmp_base + ".q.fa"
        write_fasta(qfa, batch)
        cmd = [
            "blastn","-task","megablast",
            "-query", qfa,
            "-evalue","1e-20",
            "-max_hsps","1","-max_target_seqs","1",
            "-num_threads", str(max(1, threads)),
            "-dust","no","-soft_masking","false",
            "-perc_identity", f"{pid*100:.3f}",
            "-word_size", str(WORD),
            "-outfmt","6 qseqid sseqid length pident qlen slen"
        ]
        if UNGAP: cmd.append("-ungapped")
        if use_db_local and all(os.path.exists(dbpref+e) for e in (".nhr",".nin",".nsq")):
            cmd += ["-db", dbpref]
        else:
            cmd += ["-subject", kept_file]

        proc = subprocess.run(cmd, capture_output=True, text=True)
        hits = {}
        if proc.returncode == 0 and proc.stdout:
            for ln in proc.stdout.splitlines():
                parts = ln.split("\t")
                if len(parts) < 6: continue
                qid, sid_h, L, PID, QL, SL = parts[:6]
                L=float(L); PID=float(PID); QL=float(QL); SL=float(SL)
                if (PID/100.0) >= pid and (L/QL) >= qcov and SL >= QL:
                    hits[qid.split()[0]] = (sid_h, L, PID, QL, SL)
        try: os.remove(qfa)
        except: pass
        return hits

    # 初始化 kept
    kept: List[Tuple[str,str]] = []
    if not recs:
        write_fasta(out_fa.as_posix(), [])
        return
    kept.append(recs[0]); _write_kept()
    log_seq(con, sid, recs[0][0], len(recs[0][1]), sha1(recs[0][1]), "kept", "near_seed", None)

    use_db_local = False
    accepted = 1
    bar = tqdm(total=len(recs)-1, desc=pbar_desc, unit="seq", leave=False)
    batch: List[Tuple[str,str]] = []

    for h, s in recs[1:]:
        batch.append((h, s))
        if len(batch) >= BATCH:
            if not use_db_local and _use_db_now():
                use_db_local = _make_db()
            hits = _blast_batch(batch)
            for h2, s2 in batch:
                tk = h2.split()[0]
                if tk in hits:
                    log_seq(con, sid, h2, len(s2), sha1(s2), "dropped", "near_dup", hits[tk][0])
                else:
                    kept.append((h2, s2))
                    log_seq(con, sid, h2, len(s2), sha1(s2), "kept", "near_keep", None)
                    accepted += 1
                    if use_db_local and (accepted % DB_REFRESH == 0):
                        _write_kept(); use_db_local = _make_db()
                    elif (not use_db_local) and (accepted % SUBJECT_REFRESH == 0):
                        _write_kept()
                bar.update(1)
            batch.clear()

    if batch:
        if not use_db_local and _use_db_now():
            use_db_local = _make_db()
        hits = _blast_batch(batch)
        for h2, s2 in batch:
            tk = h2.split()[0]
            if tk in hits:
                log_seq(con, sid, h2, len(s2), sha1(s2), "dropped", "near_dup", hits[tk][0])
            else:
                kept.append((h2, s2))
                log_seq(con, sid, h2, len(s2), sha1(s2), "kept", "near_keep", None)
        bar.update(len(batch))

    write_fasta(out_fa.as_posix(), kept)

    # 可选清理 tmp
    if CLEAN_TMP:
        for suf in (".kept.fa",".keptdb.nhr",".keptdb.nin",".keptdb.nsq",".q.fa"):
            p = tmp_base + suf
            try: os.remove(p)
            except: pass

# -------------------- 合并(带进度) --------------------
def merge_with_progress(unit: Path, fas: List[Path], merged_fp: Path, species_tag: str, map_fp: Path):
    total_bytes = sum((fa.stat().st_size for fa in fas if fa.exists()))
    written = 0
    bar = tqdm(total=total_bytes if total_bytes>0 else None,
               desc=f"[{species_tag}] 合并文件", unit="B", unit_scale=True, leave=False)
    with open(merged_fp, "w") as w, open(map_fp, "w") as mw:
        mw.write("safe_id\torig_header\tsource_file\tspecies_unit\n")
        for idx, fa in enumerate(fas, 1):
            bar.set_postfix_str(f"{idx}/{len(fas)} {fa.name[:30]}")
            if str(fa).endswith(".gz"):
                op = gzip.open; mode = "rt"
            else:
                op = open; mode = "r"
            with op(fa, mode, encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith(">"):
                        orig  = line.strip()[1:]
                        safe_id = f"{safe(unit_name(unit))}|{safe(fa.name)}|{safe(orig)}"
                        w.write(">" + safe_id + " " + orig + "\n")
                        mw.write(f"{safe_id}\t{orig}\t{fa.name}\t{unit.as_posix()}\n")
                    else:
                        w.write(line)
                    written += len(line.encode("utf-8"))
                    if total_bytes>0: bar.update(len(line.encode("utf-8")))
    bar.close()

def write_cat_clean_with_short_ids(cat: str,
                                   dedup_fas: List[Path],
                                   out_fa: Path,
                                   map_tsv: Path):
    """
    把各物种 *.dedup.fa 聚合为 <cat>.clean.fa，写短ID（≤50），并输出两份报告：
      - <cat>.idmap.tsv  : 短ID -> 原始信息映射
      - <cat>.bad.tsv    : 发现非法字符的序列明细
      - <cat>.bad.fa     : 问题序列的【原始序列】合集（清洗前）
    """
    bad_tsv = out_fa.with_suffix(".bad.tsv")
    bad_fa  = out_fa.with_suffix(".bad.fa")
    BAD_PRINT = int(os.getenv("BAD_PRINT", "10"))  # 终端最多打印多少条问题记录

    # 估算进度
    total_bytes = 0
    for fa in dedup_fas:
        try: total_bytes += fa.stat().st_size
        except: pass
    bar = tqdm(total=total_bytes if total_bytes > 0 else None,
               desc=f"[{cat}] 聚合(clean.fa)", unit="B", unit_scale=True, leave=False)

    # 小工具：从 safe_id 拆出 unit 与源文件名（如果你已有同名函数，可删掉此实现）
    def _parse_safe_id(safe_id: str):
        parts = safe_id.split("|", 2)
        unit = parts[0] if len(parts) > 0 else "unit"
        srcf = parts[1] if len(parts) > 1 else "file"
        return unit, srcf

    # 生成短ID（如果你文件里已有 _short_id/_compact_stitle，请继续沿用原实现）
    def _short_id(cat: str, idx: int, seq: str) -> str:
        h = hashlib.sha1(seq.encode("utf-8")).hexdigest()[:10]
        sid = f"{safe(cat)}_{idx:08d}_{h}"
        return sid[:50]

    def _compact_stitle(unit: str, srcf: str, orig: str, maxlen: int = 180) -> str:
        base = f"{unit} {srcf}"
        if orig:
            base = f"{base} | {orig}"
        return (base[:maxlen]).rstrip()

    bad_count = 0
    printed   = 0
    idx = 0

    with open(out_fa, "w") as w_clean, \
         open(map_tsv, "w") as w_map, \
         open(bad_tsv, "w") as w_bad, \
         open(bad_fa, "w") as w_badfa:

        w_map.write("short_id\tcategory\tunit\tfile\torig_header\n")
        w_bad.write("short_id\tlen\tbad_cnt\tbad_chars\tfirst_10_pos\tcategory\tunit\tfile\torig_header\n")

        for fa in dedup_fas:
            consumed_bytes = 0
            for h, s in read_fasta_lazy(fa.as_posix()):
                idx += 1
                # 拆 safe_id 与原始标题
                parts = h.split(" ", 1)
                safe_id = parts[0]
                orig    = parts[1] if len(parts) > 1 else ""
                unit, srcf = _parse_safe_id(safe_id)

                # 清洗 + 记录问题
                clean, bad_pos, bad_chars = sanitize_and_flag(s)
                sid    = _short_id(cat, idx, clean)
                stitle = _compact_stitle(unit, srcf, orig)

                # 写 clean.fa（短ID + 精简标题）
                w_clean.write(">"+sid+" "+stitle+"\n")
                for i in range(0, len(clean), 60):
                    w_clean.write(clean[i:i+60]+"\n")

                # 写映射
                w_map.write(f"{sid}\t{cat}\t{unit}\t{srcf}\t{h}\n")

                # 如果有非法字符，写 bad.tsv 与 bad.fa（原始序列）
                if bad_pos:
                    bad_count += 1
                    first10 = ",".join(str(p) for p in bad_pos[:10])
                    chars   = "".join(bad_chars)
                    w_bad.write(f"{sid}\t{len(s)}\t{len(bad_pos)}\t{chars}\t{first10}\t{cat}\t{unit}\t{srcf}\t{h}\n")
                    w_badfa.write(">"+sid+" "+stitle+" [RAW]\n")
                    for i in range(0, len(s), 60):
                        w_badfa.write(s[i:i+60]+"\n")

                    # 终端打印前 BAD_PRINT 条
                    if printed < BAD_PRINT:
                        sys.stderr.write(
                            f"[BAD] {cat}/{unit}: {sid} len={len(s)} bad={len(bad_pos)} chars={chars} pos={first10}\n"
                        )
                        printed += 1

                # 粗略推进进度（以原长度近似）
                if total_bytes:
                    consumed_bytes += len(h) + 1 + len(s) + (len(s)//60 + 1)

            if total_bytes:
                try:
                    bar.update(min(consumed_bytes, fa.stat().st_size))
                except:
                    bar.update(consumed_bytes)

    bar.close()
    if bad_count:
        sys.stderr.write(f"[AGG] {cat}: 发现 {bad_count} 条含非法字符的序列；已写 {bad_tsv.name} / {bad_fa.name}\n")
    else:
        sys.stderr.write(f"[AGG] {cat}: 未发现非法字符序列。\n")

# -------------------- 建库 --------------------
def build_db(cat: str, cat_clean_fa: Path, dbdir: Path):
    ensure_dir(dbdir)
    outbase = dbdir/safe(cat)
    logfp   = outbase.as_posix()+".makedb.log"

    total = os.path.getsize(cat_clean_fa.as_posix())
    bar = tqdm(total=total, desc=f"[{cat}] 建库(makeblastdb)", unit="B", unit_scale=True, leave=False)

    def _produced_bytes():
        s = 0
        for ext in (".nhr",".nin",".nsq"):
            p = outbase.as_posix()+ext
            if os.path.exists(p): s += os.path.getsize(p)
        return s

    # 先尝试 parse_seqids（依赖我们写入的“短ID”）
    cmd = ["makeblastdb","-in",cat_clean_fa.as_posix(),"-dbtype","nucl",
           "-parse_seqids","-hash_index","-out",outbase.as_posix()]
    proc = subprocess.Popen(cmd, stdout=open(logfp,"w"),
                            stderr=subprocess.STDOUT, start_new_session=True)
    try:
        while True:
            rc = proc.poll()
            bar.n = min(_produced_bytes(), total)
            bar.refresh()
            if rc is not None:
                break
            time.sleep(0.5)
    finally:
        bar.close()

    if proc.returncode != 0:
        sys.stderr.write(f"[WARN] parse_seqids 失败，回退不解析ID：{cat}\n")
        # 清掉半成品
        for ext in (".nhr",".nin",".nsq"):
            try: os.remove(outbase.as_posix()+ext)
            except: pass
        rc = run(["makeblastdb","-in",cat_clean_fa.as_posix(),
                  "-dbtype","nucl","-out",outbase.as_posix()])
        if rc!=0: raise RuntimeError(f"makeblastdb 失败: {cat}")

    ali = need("blastdb_aliastool")
    if ali:
        subprocess.run([ali,"-dbtype","nucl","-dblist",outbase.as_posix(),
                        "-out",outbase.as_posix()+".alias"])
        nal = outbase.as_posix()+".nal"; nlq = outbase.as_posix()+".nlq"
        try:
            if os.path.islink(nlq) or os.path.exists(nlq): os.remove(nlq)
            if os.path.exists(nal): os.symlink(os.path.basename(nal), nlq)
        except Exception: pass

# -------------------- 处理一个物种单元 --------------------
def process_one_unit(args, con, category: str, unit: Path, method: str, outdir: Path, dbdir: Path):
    """
    单物种处理（供并行调用）。返回 (unit_name, out_fa_path 或 None, 统计dict)
    """
    unit_tag = unit_name(unit)
    cat_work = outdir/safe(category)
    ensure_dir(cat_work)
    sp_key = unit.as_posix()
    # 断点
    st = species_status(con, category, sp_key)
    out_fa = cat_work/(safe(unit_tag)+".dedup.fa")
    if args.resume and out_fa.exists():
        return (unit_tag, out_fa, {"skipped": True})

    sid = species_begin(con, category, sp_key, method)
    tmp_root_base = Path(args.tmp_root) if getattr(args, "tmp_root", None) else (cat_work/"_tmp")
    ensure_dir(tmp_root_base)
    tmp_root = tmp_root_base / safe(sp_key)
    ensure_dir(tmp_root)
    log_file = tmp_root/"species.log"

    # 收集输入
    fas = collect_fastas(unit)
    if not fas:
        species_done(con, sid, 0, 0, 0, 0)
        return (unit_tag, None, {"empty": True})

    merged = tmp_root/"merged.fa"
    map_fp = cat_work/"_maps"/(safe(unit_tag)+".idmap.tsv")
    ensure_dir(map_fp.parent)
    # 1) 合并（可视化进度）
    if not merged.exists():
        merge_with_progress(unit, fas, merged, f"{category}:{unit_tag}", map_fp)
    else:
        sys.stderr.write(f"[RESUME] 复用已存在的 merged.fa -> {merged}\n")

    # (B) phase1：若已存在 phase1.fa 就直接复用
    phase1 = tmp_root/"phase1.fa"
    if not phase1.exists():
        recs = list(read_fasta_lazy(merged.as_posix()))
        n_total = len(recs)
        kept1, dropped_exact = phase1_exact(recs, sid, con, pbar_desc=f"[{category}:{unit_tag}] 精确去重")
        write_fasta(phase1.as_posix(), kept1)
        n_after_exact = len(kept1)
    else:
        sys.stderr.write(f"[RESUME] 复用已存在的 phase1.fa -> {phase1}\n")
        n_total = sum(1 for _ in read_fasta_lazy(merged.as_posix()))
        n_after_exact = sum(1 for _ in read_fasta_lazy(phase1.as_posix()))

    # --- 根据 phase1 体积决定是否改走 BLAST ---
    size = phase1.stat().st_size
    thresh = int(os.environ.get("MMSEQS_SWITCH_TO_BLAST_BYTES", 1_000_000_000))  # 默认 1GiB
    method_local = method
    if method == "mmseqs" and size >= thresh:
        sys.stderr.write(f"[INFO] {category}/{unit_tag}: phase1={size/1e9:.2f} GB >= {thresh/1e9:.2f} GB; 改用 BLAST\n")
        method_local = "blast"

    # --- Phase2：按决定的方式执行 + MMseqs -> BLAST 回退 ---
    try:
        if method_local == "mmseqs":
            try:
                mmseqs_near(phase1, out_fa, tmp_root/"mm", args.pid, args.qcov,
                            sid, con, log_file, args.threads_per_task,
                            species_tag=f"{category}:{unit_tag}")
            except Exception as e:
                sys.stderr.write(f"[WARN] {category}/{unit_tag}: MMseqs失败，回退BLAST：{e}\n")
                blast_near(phase1, out_fa, tmp_root/"blast", args.pid, args.qcov, 
                           args.threads_per_task, sid, con,
                           pbar_desc=f"[{category}:{unit_tag}] 近重复(BLAST)")
        else:
            blast_near(phase1, out_fa, tmp_root/"blast", args.pid, args.qcov, 
                       args.threads_per_task, sid, con,
                       pbar_desc=f"[{category}:{unit_tag}] 近重复(BLAST)")
    except Exception as e:
        sys.stderr.write(f"[ERROR] {category}/{unit_tag}: {e}\n")
        species_done(con, sid, n_total, n_after_exact, 0, 0)
        return (unit_tag, None, {"error": str(e)})

    # 4) 统计与收尾
    n_kept = sum(1 for _ in read_fasta_lazy(out_fa.as_posix()))
    n_dropped = n_total - n_kept
    species_done(con, sid, n_total, n_after_exact, n_kept, n_dropped)
    return (unit_tag, out_fa, {"n_total": n_total, "n_kept": n_kept, "n_drop": n_dropped})

# -------------------- 主程序 --------------------
def main():
    import argparse
    ap = argparse.ArgumentParser(description="两层去重 + 分类建库 + 并行 + 细粒度进度")
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--dbdir", required=True)
    ap.add_argument("--jobs", type=int, default=1, help="并行物种数")
    ap.add_argument("--threads-per-task", type=int, default=4, help="每个物种内部线程数（mmseqs/blast）")
    ap.add_argument("--pid", type=float, default=0.995)
    ap.add_argument("--qcov", type=float, default=0.98)
    ap.add_argument("--categories", nargs="*")
    ap.add_argument("--tmp-root", type=str, default=None,
                help="把物种的临时目录放到该路径（建议 NVMe/SSD ）")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()

    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve()
    dbdir  = Path(args.dbdir).resolve()
    ensure_dir(outdir); ensure_dir(dbdir)

    if not (need("makeblastdb") and need("blastn")):
        sys.stderr.write("需要 NCBI BLAST+ (blastn/makeblastdb)。\n"); sys.exit(2)

    method = "mmseqs" if need("mmseqs") else "blast"
    sys.stderr.write(f"[INFO] 去重方式：{method.upper()}  (pid={args.pid}, qcov={args.qcov})\n")

    db_path = outdir/"dedup_index.sqlite"
    con = db_init(db_path)
    cur = con.cursor()
    
    def _exec_retry(con, sql, params=(), tries=20, base=0.1):
        for i in range(tries):
            try:
                con.execute(sql, params)
                return
            except sqlite3.OperationalError as e:
                if "locked" in str(e).lower():
                    time.sleep(base * (i+1))  # 线性退避
                    continue
                raise
        raise

    # main() 里改成：
    with _DB_LOCK:
        _exec_retry(con, "INSERT OR REPLACE INTO meta(key,value) VALUES(?,?)", ("method", method))
        _exec_retry(con, "INSERT OR REPLACE INTO meta(key,value) VALUES(?,?)", ("pid", str(args.pid)))
        _exec_retry(con, "INSERT OR REPLACE INTO meta(key,value) VALUES(?,?)", ("qcov", str(args.qcov)))
        _exec_retry(con, "INSERT OR REPLACE INTO meta(key,value) VALUES(?,?)", ("threads_per_task", str(args.threads_per_task)))

    cats = args.categories if args.categories else \
           [p.name for p in sorted(indir.iterdir()) if p.is_dir() and p.name not in ("db","clean","__pycache__")]
    if not cats:
        sys.stderr.write("未发现类别目录。\n"); sys.exit(1)

    for ci, cat in enumerate(cats, 1):
        sys.stderr.write(f"\n=== [{ci}/{len(cats)}] 类别: {cat} ===\n")
        cat_dir = indir/cat
        units = discover_species_units(cat_dir)
        if not units:
            sys.stderr.write(f"[WARN] {cat} 无物种单元，跳过。\n"); continue

        # 并行执行
        results = []
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
            fut2name = {ex.submit(process_one_unit, args, con, cat, u, method, outdir, dbdir): unit_name(u) for u in units}
            pending = set(fut2name.keys())
            bar = tqdm(total=len(units), desc=f"[{cat}] 物种完成", unit="sp", leave=True)

            while pending and not STOP:
                done, pending = wait(pending, timeout=0.5, return_when=FIRST_COMPLETED)
                for fut in done:
                    nm = fut2name[fut]
                    try:
                        nm, outfa, info = fut.result()
                    except Exception as e:
                        sys.stderr.write(f"[ERROR] {cat}/{nm}: {e}\n"); outfa = None
                    results.append(outfa)
                    bar.update(1)

            if STOP:
                for fut in pending:
                    fut.cancel()
            bar.close()

        cleaned = [p for p in results if isinstance(p, Path) and p and p.exists()]
        if cleaned:
            cat_clean = outdir/f"{safe(cat)}.clean.fa"
            cat_map   = outdir/f"{safe(cat)}.idmap.tsv"
            write_cat_clean_with_short_ids(cat, cleaned, cat_clean, cat_map)
            build_db(cat, cat_clean, dbdir)
            sys.stderr.write(f"[OK] {cat} 已建库：{(dbdir/safe(cat)).as_posix()}.*\n")
        else:
            sys.stderr.write(f"[WARN] {cat} 无可聚合结果，跳过建库。\n")


    export_tsv(con, outdir)
    sys.stderr.write(f"\n[REPORT] kept.tsv / dropped.tsv 导出完毕；SQLite：{(outdir/'dedup_index.sqlite').as_posix()}\n")
    sys.stderr.write("[DONE]\n")

if __name__ == "__main__":
    main()