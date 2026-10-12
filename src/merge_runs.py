#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""合并多轮对接结果并去重

  对每个区域：以新一批结果替换旧一批中同编号的分子，输出合并后的数据集
  期望：每点 109 分子 × 9 姿态 = 981 条记录

做的是**文本级拼接**（不做任何坐标或属性改写），保证姿态文件逐字节可靠。
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import gzip
import os
import re
import sys
from collections import Counter

BASE = ROOT
FIX_IDS = {'016', '022', '041', '076', '077', '079', '080', '081', '082', '083', '084'}
DROP_IDS = {'089'}
SITES = ('U1', 'U3', 'U10')
SEP = '$$$$\n'


def split_records(path):
    with gzip.open(path, 'rb') as f:
        txt = f.read().decode('utf-8', 'replace')
    recs = []
    for p in txt.split(SEP):
        if p.strip() == '':
            continue
        recs.append(p + SEP)
    return recs


def rid(rec):
    return rec.split('\n', 1)[0].strip()


total = 0
for s in SITES:
    old_p = os.path.join(BASE, 'results', 'screen_%s.sdf.gz' % s)
    fix_p = os.path.join(BASE, 'results', 'refined_%s.sdf.gz' % s)
    out_p = os.path.join(BASE, 'results', 'merged_%s.sdf.gz' % s)

    old = split_records(old_p)
    fix = split_records(fix_p)

    keep = [r for r in old if rid(r) not in FIX_IDS and rid(r) not in DROP_IDS]
    dropped = Counter(rid(r) for r in old if rid(r) in FIX_IDS or rid(r) in DROP_IDS)
    merged = keep + fix

    with gzip.open(out_p, 'wb') as f:
        f.write(''.join(merged).encode('utf-8'))

    n_mol = len({rid(r) for r in merged})
    print('%-4s 原 %d 条 → 保留 %d（剔除 %d 个分子 %d 条）+ 修正版 %d 条 = %d 条 / %d 个分子'
          % (s, len(old), len(keep), len(dropped), sum(dropped.values()),
             len(fix), len(merged), n_mol))
    if len(merged) != 981 or n_mol != 109:
        print('   ⚠ 预期 981 条 / 109 分子，实际 %d 条 / %d 分子' % (len(merged), n_mol))
    total += len(merged)

print()
print('三文件合计 %d 条记录（预期 2943）' % total)
print('已写出: results/merged_{U1,U3,U10}.sdf.gz')
