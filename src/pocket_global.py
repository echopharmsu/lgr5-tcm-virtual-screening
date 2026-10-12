#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""同分辨率全局空腔对照

背景：pocket_scan.py 用 1.0 Å 网格，报告"胞外域最大空腔半径 4.69 Å"；
      u3_analysis.py 用 0.5 Å 网格在 U3 邻域实测到 5.74 Å —— 两者不可比
      （粗网格会低估距离场峰值，误差可达网格步长的半对角线 ≈0.87 Å）。
本脚本：用 **与 u3_analysis.py 相同的 0.5 Å 网格 + 相同溶剂过滤** 扫描整个胞外域，
        从而可以正确回答"U3 的腔是不是最大的"。

口径（与 pocket_scan.py 一致）：
  dmin  = 网格点到最近受体重原子的距离
  n8    = 8 Å 内的受体重原子数，n8 >= 20 才算"被受体包住"（否则是体相溶剂）

用法: python3 pocket_global_r05.py [receptor.pdb]
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import math
import os
import sys

import numpy as np

BASE = ROOT
REC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, 'data', 'receptor_apo_LGR5.pdb')
# 统一位点表**不随包提供**（其内容含作废的多候选位点框架）。
# 因此：U3 质心改从 config.yaml 取；「各区域对比」段仅在表存在时才运行。
UNI = os.path.join(BASE, 'results', 'unified_sites.tsv')   # 可选，存在才用

GRID = 0.5
MAXD = 8.0
LIN_MIN = 20

atoms = []
for l in open(REC, errors='replace'):
    if not l.startswith('ATOM'):
        continue
    el = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
    if el in ('H', 'D'):
        continue
    atoms.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
P = np.array(atoms, dtype=np.float32)
print('受体 %s : %d 重原子 ; 网格 %.1f Å' % (os.path.basename(REC), len(P), GRID))

lo = P.min(0) - MAXD - 2
hi = P.max(0) + MAXD + 2
nx = np.ceil((hi - lo) / GRID).astype(int)
print('网格 %d x %d x %d = %.1f 万点' % (nx[0], nx[1], nx[2], nx.prod() / 1e4))

dmin = np.full(tuple(nx), MAXD, dtype=np.float32)
n8 = np.zeros(tuple(nx), dtype=np.int16)
r = int(math.ceil(MAXD / GRID))
for p in P:
    c = np.floor((p - lo) / GRID).astype(int)
    i0 = np.maximum(c - r, 0)
    i1 = np.minimum(c + r + 1, nx)
    if np.any(i0 >= i1):
        continue
    gx = lo[0] + (np.arange(i0[0], i1[0]) + 0.5) * GRID
    gy = lo[1] + (np.arange(i0[1], i1[1]) + 0.5) * GRID
    gz = lo[2] + (np.arange(i0[2], i1[2]) + 0.5) * GRID
    d2 = ((gx[:, None, None] - p[0]) ** 2 + (gy[None, :, None] - p[1]) ** 2
          + (gz[None, None, :] - p[2]) ** 2)
    sl = (slice(i0[0], i1[0]), slice(i0[1], i1[1]), slice(i0[2], i1[2]))
    np.minimum(dmin[sl], np.sqrt(d2), out=dmin[sl])
    # 精确计数：窗口是 ±8 Å 的立方体，其对角可达 13.9 Å，必须用 d2 判定而不能直接 +1
    n8[sl] += (d2 <= MAXD * MAXD)

ok = n8 >= LIN_MIN
print('溶剂过滤：%d / %d 个体素通过（8 Å 内受体原子 ≥ %d）'
      % (int(ok.sum()), ok.size, LIN_MIN))
print('dmin 图完成（%.1f Å 网格），未过滤全局最大值 = %.2f Å（封顶值，体相溶剂）'
      % (GRID, dmin.max()))
print('对照：pocket_scan.py 同受体在 **1.0 Å** 网格下报 4.69 Å')
print()

# 通过过滤的全局最大
dm = np.where(ok, dmin, np.float32(-1))
gi = np.unravel_index(int(np.argmax(dm)), dm.shape)
GMAX = float(dm[gi])
GXYZ = lo + (np.array(gi) + 0.5) * GRID
gn8 = int(n8[gi])
print('=== 全局（已过滤）===')
print('  最大 dmin = %.2f Å @ (%.2f, %.2f, %.2f)   8Å内原子数 %d'
      % (GMAX, GXYZ[0], GXYZ[1], GXYZ[2], gn8))

# U3 区域（质心来自 config.yaml，不依赖外部位点表）
_u3cfg = (CFG.get('sites') or {}).get('U3') or {}
if 'center' not in _u3cfg:
    sys.exit('ERROR: config.yaml 缺少 sites.U3.center')
CEN = np.array([float(v) for v in _u3cfg['center']], dtype=np.float32)

gx, gy, gz = np.meshgrid(np.arange(nx[0]), np.arange(nx[1]), np.arange(nx[2]),
                         indexing='ij')
xyz_all = np.stack([lo[0] + (gx.ravel() + 0.5) * GRID,
                    lo[1] + (gy.ravel() + 0.5) * GRID,
                    lo[2] + (gz.ravel() + 0.5) * GRID], axis=1).astype(np.float32)
CENf = CEN.reshape(1, 3)
d3 = np.linalg.norm(xyz_all - CENf, axis=1) <= 3.0
m3 = (d3 & ok.ravel())
val3 = np.where(m3, dmin.ravel(), np.float32(-1))
i3 = int(np.argmax(val3))
U3MAX = float(val3[i3])

print()
print('=== U3 区域（区域中心 3 Å 内，同网格同过滤）===')
print('  U3 质心 (%.2f, %.2f, %.2f)' % (CEN[0], CEN[1], CEN[2]))
print('  最大 dmin = %.2f Å @ (%.2f, %.2f, %.2f)'
      % (U3MAX, xyz_all[i3][0], xyz_all[i3][1], xyz_all[i3][2]))

# 每个区域的局部最大值（同口径，3 Å 内）
print()
print('=== 各区域各自的局部最大 dmin（同口径，区域中心 3 Å 内）===')
sites = []
if os.path.exists(UNI):
    with open(UNI, encoding='utf-8') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < len(hdr):
                continue
            sites.append(dict(zip(hdr, p)))
else:
    print('  （未提供统一位点表 %s，本段跳过；U3 的局部最大 dmin 见上一节）'
          % os.path.relpath(UNI, BASE))
res = []
for s in sites:
    c = np.array([float(s['x']), float(s['y']), float(s['z'])], dtype=np.float32)
    m = (np.linalg.norm(xyz_all - c.reshape(1, 3), axis=1) <= 3.0) & ok.ravel()
    v = float(np.where(m, dmin.ravel(), np.float32(-1)).max()) if m.any() else float('nan')
    res.append((v, s['位点'], s['姿态数'], s['埋藏比']))
for v, name, npose, bur in sorted(res, reverse=True):
    print('  %-4s 局部最大 dmin %5.2f Å   姿态 %-3s 埋藏比 %s' % (name, v, npose, bur))

print()
d_gg = float(np.linalg.norm(GXYZ - CEN))
print('结论（同口径 0.5 Å 网格 + 溶剂过滤）：全局最大 %.2f Å，U3 局部 %.2f Å，相差 %.2f Å'
      % (GMAX, U3MAX, GMAX - U3MAX))
print('      全局最大点距 U3 质心 %.2f Å' % d_gg)
if abs(GMAX - U3MAX) < 0.05:
    print('      → 全胞外域的最大空腔就在 U3 区域')
elif U3MAX >= GMAX - 0.3:
    print('      → U3 的腔与全局最大同级')
else:
    print('      → U3 的腔比全局最大小 %.2f Å' % (GMAX - U3MAX))
