#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把中间 TSV 的 110 个 SMILES 展开成 RDKit 2D SDF 库（供 dock_prep.py 做 3D）。
   ID 用原表序号（三位），并保留中英文名/来源/原 Vina 分数作为 SD 属性。"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

TSV = sys.argv[1]
OUTDIR = sys.argv[2]
os.makedirs(OUTDIR, exist_ok=True)

ok, fail = 0, []
with open(TSV, encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < len(hdr):
            continue
        rec = dict(zip(hdr, p))
        smi = rec['smiles'].strip()
        m = Chem.MolFromSmiles(smi)
        if m is None:
            fail.append((rec['id'], rec['name_en'], smi[:40]))
            continue
        AllChem.Compute2DCoords(m)
        m.SetProp('_Name', rec['id'])
        m.SetProp('ID', rec['id'])
        m.SetProp('NAME_CN', rec['name_cn'])
        m.SetProp('NAME_EN', rec['name_en'])
        m.SetProp('HERB', rec['herb'])
        m.SetProp('VINA_ORIG', rec['vina_orig'])
        w = Chem.SDWriter(os.path.join(OUTDIR, '%s.sdf' % rec['id']))
        w.write(m)
        w.close()
        ok += 1

print('写出 %d 个分子到 %s' % (ok, OUTDIR))
if fail:
    print('SMILES 解析失败 %d 个: %s' % (len(fail), fail))
