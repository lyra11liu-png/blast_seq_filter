# -*- coding: utf-8 -*-
import json, re
import pandas as pd

# 简洁分类：CSV 若已有 category 列则复用；否则按名称/候选信息推断
FUNGI = {'candida','aspergillus','cryptococcus','coccidioides','histoplasma',
         'blastomyces','pneumocystis','mucor','rhizopus','trichophyton',
         'microsporum','epidermophyton','paracoccidioides','sporothrix',
         'penicillium','fusarium','scedosporium','cladosporium','alternaria',
         'nakaseomyces','pichia','exserohilum','hormonema'}
PROTOZOA = {'plasmodium','giardia','entamoeba','trichomonas','trypanosoma',
            'leishmania','toxoplasma','babesia','acanthamoeba','naegleria',
            'cryptosporidium','balantidium','cyclospora','isospora','blastocystis'}
HELMINTH = {'ascaris','taenia','echinococcus','schistosoma','trichinella',
            'strongyloides','ancylostoma','necator','enterobius','clonorchis',
            'opisthorchis','wuchereria','brugia','onchocerca','hymenolepis',
            'diphyllobothrium','fasciola'}
MYCO = {'mycoplasma','ureaplasma','acholeplasma','spiroplasma','mesoplasma','phytoplasma'}
ARCH = {'methanobrevibacter','methanosphaera','methanobacterium','halobacterium','sulfolobus','archaeoglobus'}

def _first_division(cjson):
    try:
        cand = json.loads(cjson) if isinstance(cjson, str) else None
        if cand and isinstance(cand, list):
            return (cand[0].get('division') or '').lower()
    except Exception:
        pass
    return ''

def _infer_one(text: str, division: str) -> str:
    t = (text or "").lower()
    d = (division or "").lower()
    if "virus" in t or " phage" in t or t.endswith("viridae") or "viruses" in d:
        return "Viruses"
    if any(g in t.split() for g in MYCO) or "tenericutes" in d or "mollicute" in t:
        return "Mycoplasma"
    if "archaea" in d or any(g in t for g in ARCH):
        return "Archaea"
    if any(g in t for g in FUNGI) or "fungi" in d:
        return "Fungi"
    if any(g in t for g in PROTOZOA):
        return "Protozoa"
    if any(g in t for g in HELMINTH):
        return "Helminths"
    if any(k in d for k in ["bacteria","proteobacter","firmicute","actinobacter",
                             "bacteroid","chlamyd","spirochete","cyanobacter"]):
        return "Bacteria"
    # 常见细菌属关键词兜底
    if re.search(r"(staphylo|strepto|enterococcus|escherichia|salmonella|shigella|vibrio|legionella|neisseria|listeria|brucella|yersinia|pseudomonas|acinetobacter|klebsiella|mycobacter|campylobacter|clostrid|corynebacter|borrelia|treponema|leptospira)", t):
        return "Bacteria"
    return "Bacteria"  # 兜底

def attach_or_infer_category(df: pd.DataFrame) -> pd.DataFrame:
    if "category" in df.columns:
        return df
    name_col = "best_name" if "best_name" in df.columns else ("latin_input" if "latin_input" in df.columns else None)
    parts = df[name_col].astype(str) if name_col else pd.Series([""]*len(df))
    div = df["candidates_json"].apply(_first_division) if "candidates_json" in df.columns else pd.Series([""]*len(df))
    df = df.copy()
    df["category"] = [ _infer_one(t, d) for t,d in zip(parts, div) ]
    return df
