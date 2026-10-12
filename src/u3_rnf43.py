#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""U3 的接触残基 vs LGR5–RNF43 界面残基的重合度

背景：LGR5 的 LRR 凹面是 RSPO1 / RNF43 / ZNRF3 的结合面。
      U3 位于该凹面；其 13 个接触残基与 RSPO1 / RNF43 界面均无重合，
      区域中心距 RSPO1 界面（A–M/P 界面集 36 残基）27.57 Å。
      本脚本回答：U3 是否与 **RNF43** 的结合界面重叠 → 判断
      「小分子占位可能影响 LGR5–RNF43 组装」这一假说在结构上是否可能。

口径：LGR5(A 链) 与 RNF43(E/F 链) 之间，任一重原子对距离 < 5.0 Å 即视为界面残基
      （与 README 中 RSPO1 界面残基的算法一致）。

用法: python3 u3_rnf43.py
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import csv
import os
import sys
from collections import defaultdict

import numpy as np

BASE = ROOT
FULL = os.path.join(BASE, 'data', '4KNG_full.pdb')
CONTACT_TSV = os.path.join(BASE, 'results', 'U3_接触残基与氢键.tsv')
OUT = os.path.join(BASE, 'results', 'U3_RNF43界面关系.tsv')
CUT = 5.0

if not os.path.exists(FULL):
    sys.exit('ERROR: 找不到 %s（完整复合物坐标）' % FULL)

chains = defaultdict(list)
for l in open(FULL, errors='replace'):
    if not l.startswith('ATOM'):
        continue
    el = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
    if el in ('H', 'D'):
        continue
    ch = l[21]
    chains[ch].append({'res': l[17:20].strip(), 'num': int(l[22:26]),
                       'xyz': np.array([float(l[30:38]), float(l[38:46]),
                                        float(l[46:54])])})
print('4KNG 各链重原子数: %s'
      % ', '.join('%s=%d' % (k, len(v)) for k, v in sorted(chains.items())))

A = chains.get('A', [])
if not A:
    sys.exit('ERROR: 4KNG_full.pdb 里没有 A 链')

AP = np.array([a['xyz'] for a in A])


def iface_with(other_chains):
    """返回 LGR5(A) 上与给定链任一重原子 < CUT 的残基集合"""
    O = []
    for c in other_chains:
        O.extend(chains.get(c, []))
    if not O:
        return set(), 0.0
    OP = np.array([o['xyz'] for o in O])
    res = set()
    dmin = 99.0
    for i, p in enumerate(AP):
        d = np.linalg.norm(OP - p, axis=1)
        m = float(d.min())
        if m < dmin:
            dmin = m
        if m < CUT:
            res.add('%s%d' % (A[i]['res'], A[i]['num']))
    return res, dmin


if_m, dmin_m = iface_with(['M', 'P'])
if_ef, dmin_ef = iface_with(['E', 'F'])
print()
print('LGR5(A) 与 RSPO1(M/P) 界面残基: %d 个（最近重原子距离 %.2f Å）'
      % (len(if_m), dmin_m))
print('LGR5(A) 与 RNF43(E/F) 界面残基: %d 个（最近重原子距离 %.2f Å）'
      % (len(if_ef), dmin_ef))

# U3 接触残基
u3 = []
with open(CONTACT_TSV, encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        u3.append('%s%s' % (r['残基'], r['链'] and ''))
u3set = set()
with open(CONTACT_TSV, encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        u3set.add('%s%s' % (r['残基'].strip(), r['链'].strip().replace('A', '')))
u3set = {x for x in u3set if x}
print()
print('U3 接触残基 %d 个: %s' % (len(u3set), ', '.join(sorted(u3set))))

ov_m = u3set & if_m
ov_ef = u3set & if_ef
print()
print('=== 重合度 ===')
print('  U3 ∩ RSPO1 界面 : %d 个 %s' % (len(ov_m), sorted(ov_m) if ov_m else ''))
print('  U3 ∩ RNF43 界面 : %d 个 %s' % (len(ov_ef), sorted(ov_ef) if ov_ef else ''))
print()
if len(ov_ef) == 0:
    print('判读：U3 与 RNF43 界面**完全无重合** → 「小分子占位直接阻断 LGR5–RNF43 组装」')
    print('      这一假说在结构上**不成立**（要影响组装需作用于界面本身或引起构象传递）。')
else:
    print('判读：U3 与 RNF43 界面有 %d 个残基重合 → 占位可能直接影响该界面，假说结构上可行。'
          % len(ov_ef))

with open(OUT, 'w', encoding='utf-8', newline='') as f:
    f.write('项目\t数量\t残基列表\n')
    f.write('U3 接触残基\t%d\t%s\n' % (len(u3set), ', '.join(sorted(u3set))))
    f.write('LGR5-RSPO1 界面(A链)\t%d\t%s\n' % (len(if_m), ', '.join(sorted(if_m))))
    f.write('LGR5-RNF43 界面(A链)\t%d\t%s\n' % (len(if_ef), ', '.join(sorted(if_ef))))
    f.write('U3 ∩ RSPO1 界面\t%d\t%s\n' % (len(ov_m), ', '.join(sorted(ov_m))))
    f.write('U3 ∩ RNF43 界面\t%d\t%s\n' % (len(ov_ef), ', '.join(sorted(ov_ef))))
    f.write('界面判据\t-\tLGR5(A) 与对应链任一重原子距离 < %.1f A\n' % CUT)
print('已写出: %s' % OUT)
