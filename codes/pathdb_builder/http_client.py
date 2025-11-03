#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
http_client.py
Session encapsulation 4 requests:
- Connection pool + automatic retries
- Simple rate limiting
- Unified post interface with NCBI_API_KEY
"""

import os, time, threading, requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

API_KEY = (os.getenv("NCBI_API_KEY") or "").strip()

SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "pathodb-builder/3.0 (research)",
    "Accept-Encoding": "gzip, deflate"
})

# Automatic retry policy.
_retry = Retry(
    total=5, backoff_factor=0.5,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET","POST"]
)
_adapter = HTTPAdapter(max_retries=_retry, pool_connections=32, pool_maxsize=32)
SESSION.mount("https://", _adapter)
SESSION.mount("http://", _adapter)

_rate_lock = threading.Lock()
RPS = 8 if API_KEY else 3
_last_ts = 0.0
def _throttle():
    global _last_ts
    with _rate_lock:
        now = time.time()
        min_dt = 1.0 / RPS
        if now - _last_ts < min_dt:
            time.sleep(min_dt - (now - _last_ts))
        _last_ts = time.time()
        
def http_get(url: str, params: dict):
    if API_KEY: params = {**params, "api_key": API_KEY}
    _throttle()
    r = SESSION.get(url, params=params, timeout=180)
    r.raise_for_status()
    return r

def http_post(url: str, data: dict, stream: bool=False):
    if API_KEY: data = {**data, "api_key": API_KEY}
    _throttle()
    r = SESSION.post(url, data=data, timeout=300, stream=stream)
    r.raise_for_status()
    return r