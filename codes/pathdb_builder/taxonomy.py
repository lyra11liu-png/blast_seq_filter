#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
taxonomy.py
Map species metadata to five major catalog categories: viruses/bacteria/fungi/archaea/other
"""

def map_group(meta: dict, name_hint: str="") -> str:
    lin = (meta.get("lineage","") + " " + meta.get("division","")).lower()
    hint = (name_hint or meta.get("scientificname","")).lower()
    if ("viruses" in lin) or ("virus" in hint): return "viruses"
    if "bacteria" in lin: return "bacteria"
    if ("fungi" in lin) or ("fungus" in lin): return "fungi"
    if "archaea" in lin: return "archaea"
    return "other"