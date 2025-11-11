# -*- coding: utf-8 -*-
import time, threading, requests
from typing import Iterable, Tuple

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ---- 全局限速器（防限流） ----
class _RateLimiter:
    def __init__(self, qps: float):
        self.min_gap = 1.0 / max(0.1, qps)
        self._lock = threading.Lock()
        self._next = 0.0
    def acquire(self):
        with self._lock:
            now = time.perf_counter()
            if now < self._next:
                time.sleep(self._next - now)
                now = time.perf_counter()
            self._next = now + self.min_gap

_limiter = _RateLimiter(8.0)        # 默认 8 req/s
_session = requests.Session()       # 复用连接

def init_rate_limiter(qps: float):
    global _limiter
    _limiter = _RateLimiter(qps)

def _sleep_backoff(i: int):
    time.sleep(min(60, 2 ** i))

def esearch_history(term: str, api_key: str, timeout: int, retries: int) -> Tuple[int, str, str]:
    """
    返回 (count, WebEnv, QueryKey)
    """
    params = {"db":"nuccore","term":term,"retmode":"json","usehistory":"y","retmax":0}
    if api_key: params["api_key"] = api_key
    for i in range(retries+1):
        try:
            _limiter.acquire()
            r = _session.get(f"{BASE}/esearch.fcgi", params=params, timeout=timeout)
            r.raise_for_status()
            j = r.json()
            count = int(j["esearchresult"]["count"])
            return count, j["esearchresult"]["webenv"], j["esearchresult"]["querykey"]
        except Exception:
            if i >= retries: raise
            _sleep_backoff(i)

def efetch_stream(webenv: str, query_key: str, api_key: str,
                  retstart: int, retmax: int, timeout: int, retries: int) -> Iterable[bytes]:
    """
    以 FASTA 文本流分页获取
    """
    params = {
        "db":"nuccore","query_key":query_key,"WebEnv":webenv,
        "rettype":"fasta","retmode":"text","retstart":retstart,"retmax":retmax
    }
    if api_key: params["api_key"] = api_key
    for i in range(retries+1):
        try:
            _limiter.acquire()
            with _session.get(f"{BASE}/efetch.fcgi", params=params, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                for chunk in r.iter_content(chunk_size=1<<16):
                    if chunk:
                        yield chunk
            return
        except Exception:
            if i >= retries: raise
            _sleep_backoff(i)
