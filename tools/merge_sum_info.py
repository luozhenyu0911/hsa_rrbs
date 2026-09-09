#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
合并多个两列(tab分隔)的统计文件，并转置。

输入：多个文件，每个文件格式为 "字段名<TAB>数值"
输出：一行一个字段（作为列名），一列一个样本（作为行名）

用法:
    python merge_transpose.py file1.txt file2.txt file3.txt -o merged.tsv
"""

import argparse
import sys
from collections import OrderedDict


def parse_file(path):
    """解析单个文件，返回 OrderedDict {字段名: 数值}"""
    data = OrderedDict()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            # 按第一个 tab 分割成两列
            parts = line.split("\t", 1)
            if len(parts) < 2:
                # 没有 tab 的行跳过或保留原始
                parts = [parts[0], ""]
            key = parts[0].strip()
            val = parts[1].strip()
            # 同名字段保留第一个（或可用 += 合并，这里保留首个）
            if key not in data:
                data[key] = val
    return data


def main():
    parser = argparse.ArgumentParser(
        description="合并多个两列(tab)文件并转置"
    )
    parser.add_argument("-i", "--inputf", required=True, help="the input files list: per line one file name")
    parser.add_argument("-o", "--output", default="merged.tsv", help="输出文件路径")
    args = parser.parse_args()

    # 用文件名（去掉后缀）作为样本名列名
    import os
    files_list = [f.strip() for f in open(args.inputf, "r", encoding="utf-8").readlines()]
    sample_names = [os.path.splitext(os.path.basename(p.rstrip('_sum_info.txt')))[0] for p in files_list]

    # 按顺序解析每个文件
    sample_data = OrderedDict()
    field_order = []
    for name, path in zip(sample_names, files_list):
        d = parse_file(path)
        sample_data[name] = d
        # 记录首次出现的字段顺序
        for k in d:
            if k not in field_order:
                field_order.append(k)

    # 写出：第一行是表头（字段名），第一列是样本名，转置布局
    with open(args.output, "w", encoding="utf-8") as out:
        # 表头：空一格 + 各样本名
        # out.write("\t" + "\t".join(sample_names) + "\n")
        for field in field_order:
            row = [field]
            for name in sample_names:
                row.append(sample_data[name].get(field, ""))
            out.write("\t".join(row) + "\n")

    print(f"已合并 {len(files_list)} 个文件，输出: {args.output}")
    print(f"字段数: {len(field_order)}, 样本数: {len(sample_names)}")

if __name__ == "__main__":
    main()
