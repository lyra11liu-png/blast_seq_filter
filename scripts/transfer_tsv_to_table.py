#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tsv_to_tax_table.py

功能：
  - 输入：两列的 .tsv 文件
        Sample\tTaxonomy
    其中 Taxonomy 形如：
        "Bacteria: Escherichia coli(309536) | Shigella boydii(2685) | ..."
        "Fungi: No results."
  - 输出：类似图 2 的表格，每一行一个 group：
        Sample  Group      Taxon1                           Taxon2 ...
        样本1   Bacteria   Escherichia coli(309536)        Shigella boydii(2685)
        样本1   Fungi      No results.
        样本1   Viruses    HSV-1 2(710)                    ...

用法示例：
    python tsv_to_tax_table.py \
        --in-tsv  input.tsv \
        --out-xlsx  output.xlsx
"""

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_taxonomy(sample: str, tax_str: str) -> dict:
    """
    把一行 "Sample, Taxonomy" 解析成一个 dict：
        {'Sample': sample, 'Group': 'Bacteria', 'Taxon1': 'xxx', 'Taxon2': 'yyy', ...}
    """
    tax_str = str(tax_str).strip()

    # 拆出 Group 和 后面的内容
    m = re.match(r"([^:]+):\s*(.*)", tax_str)
    if m:
        group = m.group(1).strip()
        rest = m.group(2).strip()
    else:
        # 没有冒号就全部当内容
        group = ""
        rest = tax_str

    row = {
        "Sample": sample,
        "Group": group if group else "Unknown",
    }

    # 没有命中
    if rest == "" or rest.lower().startswith("no results"):
        row["Taxon1"] = "No results."
        return row

    # 用竖线分隔多个物种
    items = [x.strip() for x in rest.split("|") if x.strip()]

    for i, item in enumerate(items, start=1):
        row[f"Taxon{i}"] = item

    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-tsv", required=True, help="输入的 tsv 文件，包含 Sample 和 Taxonomy 两列")
    ap.add_argument("--out-xlsx", required=True, help="输出的 Excel 文件，比如 result_tax_table.xlsx")
    args = ap.parse_args()

    in_path = Path(args.in_tsv)
    out_path = Path(args.out_xlsx)

    df = pd.read_csv(in_path, sep="\t", dtype=str)

    rows = []
    for _, r in df.iterrows():
        sample = r.get("Sample", r.iloc[0])
        tax_str = r.get("Taxonomy", r.iloc[1])
        rows.append(parse_taxonomy(sample, tax_str))

    out_df = pd.DataFrame(rows)

    out_df = out_df.sort_values(["Sample", "Group"], ignore_index=True)

    # 导出 Excel
    out_df.to_excel(out_path, index=False)
    print(f"✅ 写出：{out_path}")


if __name__ == "__main__":
    main()
