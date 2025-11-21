#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
utils.py
General utility functions:
- Run external commands: minimap2/blastn/samtools
- Simple timer 4 monitoring computation time

Author: lyra Liu
"""
from __future__ import annotations
import subprocess
import time
from datetime import datetime
from typing import List, Optional, Callable

class CommandError(RuntimeError):
    pass

def run_cmd(cmd: List[str] | str,
            logger: Optional[Callable[[str], None]] = None,
            shell: bool = False) -> None:
    """
    Execute external commands & directly throm exceptions when errors occur.
    :params cmd: command to be executed.
    :params logger: log function, used 4 recording commands.
    :params shell: whether to execute using the shell. e.g. 4 complex commands requiring pipes.
    """
    if logger:
        logger(f"[CMD] {cmd}")
    try:
        p = subprocess.run(
            cmd,
            shell=shell,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
    except Exception as e:
        raise CommandError(f"Failed to run command: {cmd}\nError: {e}") from e
    
    if p.returncode != 0:
        msg = (f"[CMD FAILED] {cmd}\nRETURN CODE: {p.returncode}\n"
               f"STDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")
        if logger:
            logger(msg)
        raise CommandError(msg)
    else:
        if logger:
            logger(f"[CMD DONE] {cmd}")
            
class Timer:
    def __init__(self, name: str, logger: Optional[Callable[[str], None]] = None):
        self.name = name
        self.logger = logger
        self.elapsed: float = 0.0
        
    def __enter__(self):
        self.start = time.time()
        if self.logger:
            self.logger(f"[START] {self.name} at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.elapsed = time.time() - self.start
        if self.logger:
            self.logger(f"[END] {self.name} elapsed {self.elapsed:.1f} s")
        return False
    
def default_logger_factory(log_path: str) -> Callable[[str], None]:
    def _log(msg: str) -> None:
        line = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {msg}"
        print(line)
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    return _log