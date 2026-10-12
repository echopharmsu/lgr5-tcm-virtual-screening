#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
几何口袋检测（fpocket 思路的轻量实现）：
  1) 1 Å 网格，计算每点到最近受体重原子的距离 dmin
  2) 候选口袋中心 = dmin 的局部极大值且 dmin >= MIN_R（能容纳一个探针球）
  3) 包埋度 enclosure = 从中心向 64 个方向射射线，被受体原子挡住的比例
  4) 口袋体积 = 从中心泛洪 dmin >= PROBE 的连通自由体素，并检查是否与盒边界相连（=向溶剂敞开）
  5) 打分 = 综合 dmin、enclosure、体积

用法: python3 pocket_scan.py <receptor.pdb> [out_prefix] [--min-r 3.2] [--probe 1.4]
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import math
import os
import sys
from collections import deque

import numpy as np

REC = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith('--') else \
    os.path.splitext(os.path.basename(REC))[0] + '_pockets'
MIN_R = 3.2
PROBE = 1.4
GRID = 1.0
MAXD = 8.0

argv = sys.argv[3:] if (len(sys.argv) > 2 and not sys.argv[2].startswith('--')) else sys.argv[2:]
for i, a in enumerate(argv):
    if a == '--min-r':
        MIN_R = float(argv[i + 1])
    elif a == '--probe':
        PROBE = float(argv[i + 1])


def load_atoms(path):
    xs = []
    for l in open(path, errors='replace'):
        if not l.startswith('ATOM'):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
        if el in ('H', 'D'):
            continue
        xs.append((l[17:20].strip(), int(l[22:26]),
                   (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return xs


atoms = load_atoms(REC)
P = np.array([a[2] for a in atoms], dtype=np.float32)
print('受体 %s : %d 重原子' % (os.path.basename(REC), len(P)))

lo = P.min(0) - MAXD - 2
hi = P.max(0) + MAXD + 2
nx = np.ceil((hi - lo) / GRID).astype(int)
print('网格 %d x %d x %d = %.1f 万点 (%.1f Å)' % (nx[0], nx[1], nx[2], nx.prod() / 1e4, GRID))

# ---------- dmin 图 ----------
dmin = np.full(tuple(nx), MAXD, dtype=np.float32)
r = int(math.ceil(MAXD / GRID))
for p in P:
    i0 = np.floor((p - lo) / GRID).astype(int) - r
    i1 = np.floor((p - lo) / GRID).astype(int) + r + 1
    i0 = np.maximum(i0, 0); i1 = np.minimum(i1, nx)
    if np.any(i0 >= i1):
        continue
    gx = lo[0] + (np.arange(i0[0], i1[0]) + 0.5) * GRID
    gy = lo[1] + (np.arange(i0[1], i1[1]) + 0.5) * GRID
    gz = lo[2] + (np.arange(i0[2], i1[2]) + 0.5) * GRID
    d2 = ((gx[:, None, None] - p[0]) ** 2 + (gy[None, :, None] - p[1]) ** 2
          + (gz[None, None, :] - p[2]) ** 2)
    sub = dmin[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]]
    np.minimum(sub, np.sqrt(d2), out=sub)
print('dmin 图完成: 最大空腔半径 %.2f Å' % dmin.max())

# ---------- 局部极大值（须在空腔内，排除体相溶剂）----------
cand = []
LIN_MIN = 20          # 中心 8 Å 内至少要有这么多受体重原子（否则是溶剂区）
for i in range(1, nx[0] - 1):
    for j in range(1, nx[1] - 1):
        for k in range(1, nx[2] - 1):
            v = dmin[i, j, k]
            if v < MIN_R or v >= MAXD - 0.05:      # 排除体相（被 MAXD 截断的平坦区）
                continue
            nb = dmin[i - 1:i + 2, j - 1:j + 2, k - 1:k + 2]
            if v >= nb.max() - 1e-6:
                cand.append((i, j, k, float(v)))
cand.sort(key=lambda c: -c[3])
# 近邻去重（同一空腔只留最大者）
keep = []
for c in cand:
    pos = lo + (np.array(c[:3]) + 0.5) * GRID
    if all(np.linalg.norm(pos - kp) > 4.0 for kp in keep):
        keep.append(pos)
    if len(keep) >= 80:
        break
print('候选空腔中心: %d 个（dmin>=%.1f Å 且非体相，4 Å 内去重）' % (len(keep), MIN_R))

# ---------- 原子格（射线与泛洪用）----------
CELL = 4.0
grid_atoms = {}
for idx, p in enumerate(P):
    k = tuple(((p - lo) // CELL).astype(int))
    grid_atoms.setdefault(k, []).append(idx)

VDW = 1.7          # 重原子范德华半径近似
RAY_N = 64
RAY_LEN = 14.0
RAY_STEP = 0.5


def neighbors(p, rad):
    """返回 p 附近 rad 内的原子索引"""
    k0 = np.floor((p - rad - lo) / CELL).astype(int)
    k1 = np.floor((p + rad - lo) / CELL).astype(int) + 1
    out = []
    for i in range(k0[0], k1[0] + 1):
        for j in range(k0[1], k1[1] + 1):
            for k in range(k0[2], k1[2] + 1):
                out.extend(grid_atoms.get((i, j, k), ()))
    return out


# Fibonacci 球面方向
dirs = []
ga = math.pi * (3 - math.sqrt(5))
for i in range(RAY_N):
    y = 1 - (i / (RAY_N - 1)) * 2
    rad = math.sqrt(max(0.0, 1 - y * y))
    th = ga * i
    dirs.append(np.array([math.cos(th) * rad, y, math.sin(th) * rad], dtype=np.float32))
dirs = np.array(dirs)


def enclosure(center):
    blocked = 0
    for d in dirs:
        hit = False
        for t in np.arange(PROBE + VDW, RAY_LEN, RAY_STEP):
            pt = center + d * t
            for ai in neighbors(pt, VDW + 0.2):
                if np.linalg.norm(P[ai] - pt) < VDW:
                    hit = True
                    break
            if hit:
                break
        blocked += 1 if hit else 0
    return blocked / len(dirs)


def label_components():
    """对 dmin >= PROBE 的自由体素做一次全局连通标记，返回 (标签数组, 每个分量的体积与是否触边界)"""
    lab = np.full(tuple(nx), -1, dtype=np.int32)
    stats = {}
    free = dmin >= PROBE
    cur = 0
    for i0 in range(nx[0]):
        for j0 in range(nx[1]):
            for k0 in range(nx[2]):
                if not free[i0, j0, k0] or lab[i0, j0, k0] >= 0:
                    continue
                q = deque([(i0, j0, k0)])
                lab[i0, j0, k0] = cur
                n = 0
                touches = False
                while q:
                    t = q.popleft()
                    n += 1
                    if t[0] in (0, nx[0] - 1) or t[1] in (0, nx[1] - 1) or t[2] in (0, nx[2] - 1):
                        touches = True
                    for dd in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                        a, b, c = t[0] + dd[0], t[1] + dd[1], t[2] + dd[2]
                        if 0 <= a < nx[0] and 0 <= b < nx[1] and 0 <= c < nx[2] \
                                and free[a, b, c] and lab[a, b, c] < 0:
                            lab[a, b, c] = cur
                            q.append((a, b, c))
                stats[cur] = (n * GRID ** 3, touches)
                cur += 1
    return lab, stats


print('连通自由空间标记中 ...')
LAB, CSTATS = label_components()
print('自由空间分量数: %d（含体相溶剂）' % len(CSTATS))


def pocket_volume_local(center, rad=12.0):
    """有界泛洪：只统计中心 rad 半径内、dmin>=PROBE 的连通自由体素体积（不受与溶剂连通影响）"""
    c = np.floor((center - lo) / GRID).astype(int)
    q = deque([tuple(c)])
    seen = {tuple(c)}
    while q:
        t = q.popleft()
        for dd in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            a, b, cc = t[0] + dd[0], t[1] + dd[1], t[2] + dd[2]
            if not (0 <= a < nx[0] and 0 <= b < nx[1] and 0 <= cc < nx[2]):
                continue
            if (a, b, cc) in seen or dmin[a, b, cc] < PROBE:
                continue
            pos = lo + (np.array([a, b, cc]) + 0.5) * GRID
            if np.linalg.norm(pos - center) > rad:
                continue
            seen.add((a, b, cc))
            q.append((a, b, cc))
    return len(seen) * GRID ** 3


rows = []
for pos in keep:
    nlin_all = len(neighbors(pos, 8.0))
    if nlin_all < LIN_MIN:
        continue
    enc = enclosure(pos)
    c = tuple(np.floor((pos - lo) / GRID).astype(int))
    l = LAB[c]
    _, openb = CSTATS.get(int(l), (0.0, True)) if l >= 0 else (0.0, True)
    vol = pocket_volume_local(pos)
    lin = {}
    for ai in neighbors(pos, 8.0):
        d = np.linalg.norm(P[ai] - pos)
        key = '%s%d' % (atoms[ai][0], atoms[ai][1])
        lin[key] = min(lin.get(key, 99), d)
    score = enc * min(1.0, vol / 400.0)
    rows.append({'pos': pos, 'r': float(dmin[c]), 'enc': enc, 'vol': vol, 'open': openb,
                 'score': score, 'nlin': len(lin),
                 'lin': ','.join(k for k, _ in sorted(lin.items(), key=lambda kv: kv[1])[:8])})

# 分类：内部空腔（体积异常大且几乎全封闭）vs 表面口袋
surface = [r for r in rows if r['vol'] <= 3000 and r['enc'] <= 0.93]
interior = [r for r in rows if r not in surface]
surface.sort(key=lambda r: -r['score'])
interior.sort(key=lambda r: -r['score'])
rows = surface

print()
print('=== 表面口袋候选（已排除蛋白内部空腔）===')
print('%-5s %22s %7s %8s %10s %7s %6s %s'
      % ('#', '中心 (x, y, z)', '空腔半径', '包埋度', '体积Å³', '敞口', '衬里', '衬里残基'))
print('-' * 120)
for i, r in enumerate(rows[:20], 1):
    print('%-5d (%7.2f,%7.2f,%7.2f) %7.2f %8.2f %10.0f %7s %6d %s'
          % (i, r['pos'][0], r['pos'][1], r['pos'][2], r['r'], r['enc'], r['vol'],
             '是' if r['open'] else '否', r['nlin'], r['lin'][:60]))
print()
print('=== 蛋白内部空腔（参考，不作为小分子结合区域）===')
for i, r in enumerate(interior[:5], 1):
    print('%-5d (%7.2f,%7.2f,%7.2f) %7.2f %8.2f %10.0f'
          % (i, r['pos'][0], r['pos'][1], r['pos'][2], r['r'], r['enc'], r['vol']))
print()
print('类药口袋参考线（fpocket 经验）: 空腔半径 >= 3.5 Å 且 包埋度 0.55-0.90 且 体积 200-1500 Å³')
n_drug = sum(1 for r in rows if r['r'] >= 3.5 and 0.55 <= r['enc'] <= 0.90 and 200 <= r['vol'] <= 1500)
print('本受体满足该判据的区域数: %d' % n_drug)

with open(OUT + '.tsv', 'w', encoding='utf-8') as f:
    f.write('类别\t排名\tx\ty\tz\t空腔半径A\t包埋度\t体积A3\t向溶剂敞开\t衬里残基数\t衬里残基\n')
    for i, r in enumerate(rows, 1):
        f.write('表面口袋\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\t%.0f\t%s\t%d\t%s\n'
                % (i, r['pos'][0], r['pos'][1], r['pos'][2], r['r'], r['enc'], r['vol'],
                   'yes' if r['open'] else 'no', r['nlin'], r['lin']))
    for i, r in enumerate(interior, 1):
        f.write('内部空腔\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\t%.0f\t%s\t%d\t%s\n'
                % (i, r['pos'][0], r['pos'][1], r['pos'][2], r['r'], r['enc'], r['vol'],
                   'yes' if r['open'] else 'no', r['nlin'], r['lin']))
print()
print('已写出: %s.tsv' % OUT)
print('判读提示: 真正的类药口袋通常满足 空腔半径>=3.5 Å + 包埋度>=0.55 + 体积>=300 Å³ + 不向溶剂敞开')
