#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""由给定的 SMILES 清单建立配体库：
   1) 生成子库（指定编号的分子，使用显式指定立体构型的 SMILES）
   2) 生成归一化总表（剔除无法确定构型的条目）
   并自检：所用 SMILES 的未定义手性中心必须为 0。
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import csv
import os
import sys

from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

BASE = ROOT
SRC = HERE
FIXED = path('data', 'SMILES_修正版_110.txt')
sys.path.insert(0, SRC)
try:
    import stereo_gate as SG
except ImportError:
    SG = None
    print('ERROR: 找不到 stereo_gate.py，无法做构型自检')
    sys.exit(1)

FIX_IDS = ['016', '022', '041', '076', '077', '079', '080', '081', '082', '083', '084']
DROP_IDS = ['089']                      # 36 个手性中心全未定义，用户决定排除

lib_path = os.path.join(BASE, 'work', 'screen', 'ligands_source.tsv')
with open(lib_path, encoding='utf-8') as f:
    lib = list(csv.DictReader(f, delimiter='\t'))
cols = list(lib[0].keys())
fixed = [l.strip() for l in open(FIXED, encoding='utf-8') if l.strip()]
assert len(lib) == 110 and len(fixed) == 110, '行数不是 110'

# 修正版按行号对应序号 1..110
smi_new = {lib[i]['id']: fixed[i] for i in range(110)}

print('=== 修正版 11 个分子的构型自检 ===')
bad = 0
for rid in FIX_IDS:
    m = Chem.MolFromSmiles(smi_new[rid])
    st = SG.analyze(m)
    ik = Chem.MolToInchiKey(m)
    flag = 'OK' if not st['ambiguous'] else '⚠ 仍有未定义中心 %d' % st['n_undefined']
    if st['ambiguous']:
        bad += 1
    name = [r['name_cn'] for r in lib if r['id'] == rid][0]
    print('  %s %-12s 中心 %d 已定义 %d 未定义 %d  %s  %s'
          % (rid, name, st.get('n_real', 0), st['n_defined'], st['n_undefined'],
             ik, flag))
if bad:
    print('⚠ 有 %d 个分子修正后仍有未定义中心 —— 需复核' % bad)
else:
    print('全部 11 个分子的手性中心均已完整定义 ✓')

# ---- 写出 11 分子子库 ----
fixdir = os.path.join(BASE, 'work', 'screen', 'library_subset')
os.makedirs(fixdir, exist_ok=True)
p11 = os.path.join(fixdir, 'ligands_source.tsv')
with open(p11, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
    w.writeheader()
    for r in lib:
        if r['id'] in FIX_IDS:
            r2 = dict(r)
            r2['smiles'] = smi_new[r['id']]
            w.writerow(r2)
print()
print('写出 11 分子子库: %s' % p11)

# ---- 写出 109 行 v2 总表 ----
pv2 = os.path.join(BASE, 'work', 'screen', 'ligands_source_v2.tsv')
n_fix = n_drop = 0
with open(pv2, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
    w.writeheader()
    for r in lib:
        if r['id'] in DROP_IDS:
            n_drop += 1
            continue
        r2 = dict(r)
        if r['id'] in FIX_IDS:
            r2['smiles'] = smi_new[r['id']]
            n_fix += 1
        w.writerow(r2)
print('写出 v2 总表: %s （%d 行；换修正版 %d；剔除 %d）'
      % (pv2, 110 - n_drop, n_fix, n_drop))
