#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""U3 区域专项补算（一次跑完三件事）

  1) 区域成员重判：按 unified_sites.tsv 的 U3 质心，从各轮结果里挑出属于 U3 的姿态
     （判据：配体重原子质心距 U3 质心 ≤ 5.0 Å）
  2) 残基级接触表 + 氢键明细
     接触 = 重原子间 ≤ 4.0 Å ；氢键 = N/O/S 之间 ≤ 3.5 Å（与 analyze_blind.py 同口径）
  3) U3 空腔几何：局部 dmin 最大值 / 包埋度(64 射线) / 衬里残基
     dmin 口径与 pocket_scan.py 完全一致（1 Å 网格、到最近受体重原子的距离）
     参照：LGR5 胞外域全局 dmin 最大 4.69 Å（落在 U10）；9UVZ 已验证口袋 6.01 Å

用法: python3 u3_analysis.py [receptor.pdb]
      （U3 质心取自 config.yaml 的 sites.U3.center，本脚本**不需要**随包外的统一位点表）
输出: results/U3_接触残基与氢键.tsv / results/U3_空腔性质.txt
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import glob
import gzip
import math
import os
import sys
from collections import defaultdict

import numpy as np
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

BASE = ROOT
REC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, 'data', 'receptor_apo_LGR5.pdb')
# 说明：本脚本不再读取 unified_sites.tsv（该表不随包提供），质心改从 config.yaml 取

CONTACT = 4.0
HB = 3.5
CLUSTER = 5.0
POLAR = {'N', 'O', 'S'}

# ---------------- 受体 ----------------
rec = []
for l in open(REC, errors='replace'):
    if not l.startswith('ATOM'):
        continue
    el = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
    if el in ('H', 'D'):
        continue
    rec.append({'res': l[17:20].strip(), 'num': int(l[22:26]), 'ch': l[21],
                'name': l[12:16].strip(), 'el': el,
                'xyz': np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])})
print('受体 %s : %d 重原子' % (os.path.basename(REC), len(rec)))
RP = np.array([a['xyz'] for a in rec])

# ---------------- U3 质心（来自 config.yaml，不依赖随包外的位点表）----------------
_u3cfg = (CFG.get('sites') or {}).get('U3') or {}
if 'center' not in _u3cfg:
    sys.exit('ERROR: config.yaml 缺少 sites.U3.center')
CEN = np.array([float(v) for v in _u3cfg['center']])
print('U3 质心 (%.2f, %.2f, %.2f)（取自 config.yaml 的 sites.U3.center）'
      % (CEN[0], CEN[1], CEN[2]))


def load_poses(path):
    """返回 [(来源轮次, CNNscore, affinity, {元素:xyz 列表})]"""
    out = []
    tag = os.path.basename(path).replace('.sdf.gz', '')
    with gzip.open(path, 'rb') as fh:
        supp = Chem.ForwardSDMolSupplier(fh, removeHs=False, sanitize=False)
        for m in supp:
            if m is None:
                continue
            conf = m.GetConformer()
            atoms = []
            for a in m.GetAtoms():
                if a.GetSymbol().upper() in ('H', 'D'):
                    continue
                p = conf.GetAtomPosition(a.GetIdx())
                deg = sum(1 for nb in a.GetNeighbors()
                          if nb.GetSymbol().upper() not in ('H', 'D'))
                dbl = any(b.GetBondTypeAsDouble() > 1.5 for b in a.GetBonds()
                          if b.GetOtherAtom(a).GetSymbol().upper() not in ('H', 'D'))
                atoms.append({'el': a.GetSymbol().upper(),
                              'xyz': np.array([p.x, p.y, p.z]),
                              'idx': a.GetIdx(), 'deg': deg, 'dbl': dbl})
            if not atoms:
                continue
            cnn = float(m.GetProp('CNNscore')) if m.HasProp('CNNscore') else float('nan')
            aff = (float(m.GetProp('minimizedAffinity'))
                   if m.HasProp('minimizedAffinity') else float('nan'))
            out.append((tag, cnn, aff, atoms))
    return out


def contacts(atoms):
    """返回 (逐残基字典, 最近接触, 氢键明细)"""
    per = defaultdict(lambda: [99.0, 0, 0, None])   # key -> [最近距, 原子数, 氢键数, 氢键明细]
    dmin = 99.0
    hbonds = []
    for at in atoms:
        el, xyz = at['el'], at['xyz']
        d = np.linalg.norm(RP - xyz, axis=1)
        i = int(d.argmin())
        if d[i] < dmin:
            dmin = float(d[i])
        for j in np.where(d <= CONTACT)[0]:
            a = rec[j]
            k = '%s%d%s' % (a['res'], a['num'], a['ch'])
            c = per[k]
            if d[j] < c[0]:
                c[0] = float(d[j])
            c[1] += 1
            if el in POLAR and a['el'] in POLAR and d[j] <= HB:
                c[2] += 1
                hbonds.append((k, a['name'], el, float(d[j]), at['idx'],
                               at['deg'], at['dbl']))
    return per, dmin, hbonds


# ---------------- 收集 U3 姿态 ----------------
files = [os.path.join(BASE, 'results', 'blind_apo_e128.sdf.gz')] + \
        sorted(glob.glob(os.path.join(BASE, 'results', 'focused_S*.sdf.gz')))
members = []
for f in files:
    if not os.path.exists(f):
        continue
    for tag, cnn, aff, atoms in load_poses(f):
        cen = np.mean([a['xyz'] for a in atoms], axis=0)
        dd = float(np.linalg.norm(cen - CEN))
        if dd <= CLUSTER:
            members.append({'tag': tag, 'cnn': cnn, 'aff': aff, 'atoms': atoms, 'dc': dd})

print()
print('=== 1) U3 区域成员重判（质心 ≤%.1f Å）===' % CLUSTER)
byrun = defaultdict(int)
for m in members:
    byrun[m['tag']] += 1
for k in sorted(byrun):
    print('  %-22s %d 个姿态' % (k, byrun[k]))
print('  合计 %d 个姿态（区域表记载 22，含 blind_apo_e128 与 focused 轮）' % len(members))

if not members:
    sys.exit('ERROR: 没找到属于 U3 的姿态')

best_cnn = max(members, key=lambda m: m['cnn'])
best_aff = min(members, key=lambda m: m['aff'])
print('  最佳 CNN : %s  CNN %.4f  aff %.2f  (距质心 %.2f Å)'
      % (best_cnn['tag'], best_cnn['cnn'], best_cnn['aff'], best_cnn['dc']))
print('  最佳 aff : %s  CNN %.4f  aff %.2f  (距质心 %.2f Å)'
      % (best_aff['tag'], best_aff['cnn'], best_aff['aff'], best_aff['dc']))

# ---------------- 2) 接触表 ----------------
per_c, dmin_c, hb_c = contacts(best_cnn['atoms'])
per_a, dmin_a, hb_a = contacts(best_aff['atoms'])

print()
print('=== 2) 残基级接触（最佳 CNN 姿态，%s）===' % best_cnn['tag'])
print('  接触残基 %d 个 ; 最近接触 %.2f Å ; 氢键残基 %d 个'
      % (len(per_c), dmin_c, sum(1 for v in per_c.values() if v[2])))
print('  %-10s %8s %6s %6s' % ('残基', '最近距Å', '原子数', '氢键'))
for k, v in sorted(per_c.items(), key=lambda kv: kv[1][0]):
    print('  %-10s %8.2f %6d %6d' % (k, v[0], v[1], v[2]))

print()
print('=== 极性接触明细（最佳 CNN 姿态）===')
print('  判据：N/O/S 之间 ≤ %.1f Å（几何判据，不显式判定供受体）' % HB)
print('  %-10s %-6s %-8s %6s %8s  %-22s %s'
      % ('受体残基', '受体原子', '配体原子', '重邻数', '距离Å', '配体侧性质', '供受体判读'))
REC_DON = {'N', 'NZ', 'NH1', 'NH2', 'NE', 'ND2', 'NE2', 'NE1',
           'OG', 'OG1', 'OH', 'ND1'}
REC_ACC = {'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH',
           'ND1', 'NE2', 'SD'}
for h in sorted(hb_c, key=lambda x: x[3]):
    k, aname, el, dd, idx, deg, dbl = h
    if deg == 0 and not dbl:
        lig_kind, lig_don = '羟基 O（可作供体）', True
    elif dbl:
        lig_kind, lig_don = '羰基 O（仅受体）', False
    else:
        lig_kind, lig_don = '醚/甲氧基 O（仅受体）', False
    rec_don = aname in REC_DON
    rec_acc = aname in REC_ACC
    ok = (rec_don and True) or (lig_don and rec_acc)
    verdict = '几何上合理的氢键' if ok else '两侧均为受体 → 非氢键（静电不利）'
    print('  %-10s %-6s %-8s %6d %8.2f  %-22s %s'
          % (k, aname, '%s%d' % (el, idx), deg, dd, lig_kind, verdict))

# 两轮重合度（按残基集合）
def rset(per):
    return {k[:3] + k[3:-1] for k in per}      # 去掉链字母 -> 'SER402'

s_c, s_a = rset(per_c), rset(per_a)
inter = s_c & s_a
jac = len(inter) / max(1, len(s_c | s_a))
print()
print('=== 双姿态一致性（最佳CNN vs 最佳aff）===')
print('  接触残基集合 Jaccard = %.3f  （交集 %d / 并集 %d）' % (jac, len(inter), len(s_c | s_a)))

# ---------------- 3) 空腔几何 ----------------
print()
print('=== 3) U3 空腔几何 ===')
R8 = 8.0
STEP = 0.5
LIN_MIN = 20          # 与 pocket_scan.py 完全一致：8 Å 内至少 20 个受体重原子，否则视为体相溶剂
rng = np.arange(-R8, R8 + 1e-9, STEP)
gx, gy, gz = np.meshgrid(CEN[0] + rng, CEN[1] + rng, CEN[2] + rng, indexing='ij')
G = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
keep = np.linalg.norm(G - CEN, axis=1) <= R8
G = G[keep]
dmin_local = np.full(len(G), 99.0)
n8 = np.zeros(len(G), dtype=np.int32)
for p in RP:
    dd = np.linalg.norm(G - p, axis=1)
    np.minimum(dmin_local, dd, out=dmin_local)
    n8 += (dd <= R8)
# 溶剂过滤：只在"被受体包住"的点里取最大值（否则会取到蛋白外的溶剂空隙）
ok = n8 >= LIN_MIN
n_solv = int((~ok).sum())
d3 = (np.linalg.norm(G - CEN, axis=1) <= 3.0) & ok
if not ok.any():
    print('  ⚠ 8 Å 域内没有任何点通过溶剂过滤（%d 个点被排除）' % n_solv)
    r_local = float('nan')
    PMAX = CEN
else:
    r_local = float(dmin_local[ok].max())
    PMAX = G[int(np.where(ok)[0][int(dmin_local[ok].argmax())])]
r_site = float(dmin_local[d3].max()) if d3.any() else float('nan')
print('  溶剂过滤：%d 个网格点 / 共 %d 个（8 Å 内受体原子 < %d 者剔除，与 pocket_scan.py 同规则）'
      % (int(ok.sum()), len(G), LIN_MIN))
print('  局部 dmin 最大（8 Å 域内，已过滤）: %.2f Å  @ (%.2f, %.2f, %.2f)'
      % (r_local, PMAX[0], PMAX[1], PMAX[2]))
print('  区域中心 3 Å 内 dmin 最大（已过滤）: %.2f Å' % r_site)
print('  参照：LGR5 胞外域全局最大 4.69 Å（落在 U10）；9UVZ 已验证口袋 6.01 Å')

# 包埋度：64 射线
dirs = []
for i in range(8):
    for j in range(8):
        th = math.pi * (i + 0.5) / 8
        ph = 2 * math.pi * (j + 0.5) / 8
        dirs.append(np.array([math.sin(th) * math.cos(ph),
                              math.sin(th) * math.sin(ph), math.cos(th)]))
blocked = 0
for d in dirs:
    hit = False
    for t in np.arange(0.5, 12.0, 0.5):
        q = PMAX + d * t
        if np.linalg.norm(RP - q, axis=1).min() < 1.8:
            hit = True
            break
    blocked += hit
enc = blocked / len(dirs)
print('  包埋度（64 射线，%s）      : %.2f' % ('中心最大空腔点', enc))

# 衬里残基：距最大空腔点 (r_local + 1.4) Å 内的残基
cut = r_local + 1.4
d = np.linalg.norm(RP - PMAX, axis=1)
lining = sorted({('%s%d' % (rec[i]['res'], rec[i]['num'])) for i in np.where(d <= cut)[0]})
print('  衬里残基（≤%.1f Å，%d 个）   : %s' % (cut, len(lining), ', '.join(lining)))

# ---------------- 写出 ----------------
out_tsv = os.path.join(BASE, 'results', 'U3_接触残基与氢键.tsv')
with open(out_tsv, 'w', encoding='utf-8', newline='') as f:
    f.write('残基\t链\t最近距A\t接触原子数\t氢键数\t'
            '最佳CNN姿态距离A\t最佳aff姿态距离A\t最佳aff姿态氢键数\n')
    for k, v in sorted(per_c.items(), key=lambda kv: kv[1][0]):
        va = per_a.get(k, [float('nan'), 0, 0, None])
        f.write('%s\t%s\t%.2f\t%d\t%d\t%.2f\t%.2f\t%d\n'
                % (k[:-1], k[-1], v[0], v[1], v[2], v[0], va[0], va[2]))
print()
print('已写出: %s' % out_tsv)

out_txt = os.path.join(BASE, 'results', 'U3_空腔性质.txt')
with open(out_txt, 'w', encoding='utf-8') as f:
    f.write('U3 位点空腔几何（口径与 pocket_scan.py 一致：1 Å 网格 dmin = 到最近受体重原子距离）\n')
    f.write('受体: %s\n' % os.path.basename(REC))
    f.write('U3 质心: (%.2f, %.2f, %.2f)\n' % tuple(CEN))
    f.write('局部 dmin 最大(8A域内): %.2f A @ (%.2f, %.2f, %.2f)\n'
            % (r_local, PMAX[0], PMAX[1], PMAX[2]))
    f.write('位点中心 3A 内 dmin 最大: %.2f A\n' % r_site)
    f.write('包埋度(64射线): %.2f\n' % enc)
    f.write('衬里残基(%d 个): %s\n' % (len(lining), ', '.join(lining)))
    f.write('参照: LGR5 胞外域全局最大 4.69 A（落在 U10）; 9UVZ 已验证口袋 6.01 A\n')
    f.write('接触残基(最佳CNN姿态) %d 个; 氢键残基 %d 个\n'
            % (len(per_c), sum(1 for v in per_c.values() if v[2])))
print('已写出: %s' % out_txt)
