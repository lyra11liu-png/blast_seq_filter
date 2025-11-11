# -*- coding: utf-8 -*-
import shutil, subprocess

def need(cmd: str):
    if shutil.which(cmd) is None:
        raise SystemExit(f"ERROR: required command not found: {cmd}")

def run(cmd: list):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise SystemExit(f"[cmd failed] {' '.join(cmd)}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")
    return p
