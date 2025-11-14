#!/usr/bin/env python
# -*. coding: utf-8 -*-

"""
host_filter.py
Author: lyra LIU
Batch Bam.: removehuman seqs, extract unmapped reads after minimap2 alignment against human reference.
Support parallel processing, multiple log, timing, & resource monitoring.
"""
from __future__ import annotations
import os, sys, re, time, json, shutil, argparse, subprocess, datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    import psutil
except Exception:
    psutil = None
    
def which(x: str) -> str:
    p = shutil.which(x)
    if not p:
        raise FileExistsError(f"[ERR] Can't find dependence: {x}")