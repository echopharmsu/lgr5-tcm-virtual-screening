#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
交付物生成：候选清单（四道过滤）+ 质检门禁 + Excel 汇总。
判据：最近接触 >= 2.5 Å、affinity < 0、CNNscore >= 0.60、距区域中心 <= 6 Å。

用法: python3 make_deliverables.py <ligands_source.tsv> <receptor.pdb> <out_prefix> <screen_U1.sdf.gz> [screen_U3...]
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import gzip
import math
import os
import sys
from collections import defaultdict

SRC_TSV, REC, OUT = sys.argv[1:4]
RESULTS = sys.argv[4:]
SITES = [os.path.basename(r).replace('screen_', '').replace('.sdf.gz', '') for r in RESULTS]
TARGET = '001'
MIN_CONTACT, MIN_CNN, SITE_CUT = 2.5, 0.60, 6.0
CONTACT, HB = 4.0, 3.5
POLAR = {'N', 'O', 'S'}

# 区域中心 = 对接盒中心（写死以保证判据独立于结果）
SITE_CENTER = {
    'U1': (23.55, -3.49, 21.44),
    'U3': (25.17, -17.39, 41.21),
    'U10': (29.97, -15.32, 62.40),
}

meta = {}
with open(SRC_TSV, encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= len(hdr):
            meta[p[0]] = dict(zip(hdr, p))

rec = []
for l in open(REC, errors='replace'):
    if not l.startswith('ATOM'):
        continue
    el = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
    if el in ('H', 'D'):
        continue
    rec.append((l[17:20].strip(), l[22:26].strip(), el,
                (float(l[30:38]), float(l[38:46]), float(l[46:54]))))


def parse(path):
    op = gzip.open if path.endswith('.gz') else open
    txt = op(path, 'rt', errors='replace').read()
    mols = defaultdict(list)
    for b in txt.split('$$$$'):
        if 'CNNscore' not in b:
            continue
        lines = b.split('\n')
        coords = []
        for i, l in enumerate(lines):
            if 'V2000' in l:
                na = int(l[0:3])
                for j in range(i + 1, i + 1 + na):
                    a = lines[j]
                    coords.append(((a[31:34].strip() or a[34:36].strip()).upper(),
                                   float(a[0:10]), float(a[10:20]), float(a[20:30])))
                break
        props = {}
        for i, l in enumerate(lines):
            if l.strip().startswith('> <'):
                k = l.split('<')[1].split('>')[0]
                if i + 1 < len(lines):
                    props[k] = lines[i + 1].strip()
        if not coords:
            continue
        mid = props.get('ID') or props.get('CAS') or lines[0].strip()
        lig = [(e, (x, y, z)) for e, x, y, z in coords if e not in ('H', 'D')]
        contacts, dmin = defaultdict(lambda: [99.0, 0, 0]), 99.0
        for el, xyz in lig:
            for res, num, rel, rxyz in rec:
                d = math.dist(xyz, rxyz)
                dmin = min(dmin, d)
                if d <= CONTACT:
                    c = contacts[res + num]
                    c[0] = min(c[0], d); c[1] += 1
                    if el in POLAR and rel in POLAR and d <= HB:
                        c[2] += 1
        cen = (sum(x for _, (x, _, _) in lig) / len(lig),
               sum(y for _, (_, y, _) in lig) / len(lig),
               sum(z for _, (_, _, z) in lig) / len(lig))
        mols[mid].append({
            'cnn': float(props.get('CNNscore', 'nan')),
            'caff': float(props.get('CNNaffinity', 'nan')),
            'aff': float(props.get('minimizedAffinity', 'nan')),
            'dmin': dmin, 'ncont': len(contacts), 'cen': cen,
            'nhb': sum(1 for c in contacts.values() if c[2]),
            'contacts': sorted(contacts.items(), key=lambda kv: kv[1][0]),
        })
    return mols


sites_data, site_centers = {}, {}
all_rows = {}
for path, site in zip(RESULTS, SITES):
    mols = parse(path)
    # 区域中心：优先用预设盒中心；未预设时退回"该轮最优姿态质心的中位数"
    if site in SITE_CENTER:
        cen = SITE_CENTER[site]
    else:
        tops = [v[0] for v in mols.values()]
        cen = tuple(sorted(t['cen'][k] for t in tops)[len(tops) // 2] for k in range(3))
    site_centers[site] = cen
    sites_data[site] = mols
    rows = []
    for mid, poses in mols.items():
        p = poses[0]
        d_site = math.dist(p['cen'], cen)
        reasons = []
        if not (p['dmin'] == p['dmin'] and p['dmin'] >= MIN_CONTACT):
            reasons.append('穿透(%.2fÅ)' % p['dmin'])
        if not (p['aff'] == p['aff'] and p['aff'] < 0):
            reasons.append('亲和力≥0')
        if not (p['cnn'] == p['cnn'] and p['cnn'] >= MIN_CNN):
            reasons.append('CNN<%.2f' % MIN_CNN)
        if not (d_site <= SITE_CUT):
            reasons.append('偏离位点(%.1fÅ)' % d_site)
        rows.append({'id': mid, 'cnn': p['cnn'], 'caff': p['caff'], 'aff': p['aff'],
                     'dmin': p['dmin'], 'd_site': d_site, 'ncont': p['ncont'], 'nhb': p['nhb'],
                     'pass': not reasons, 'reason': ';'.join(reasons) or 'OK',
                     'contacts': p['contacts']})
    rows.sort(key=lambda r: (0 if r['pass'] else 1, -(r['cnn'] if r['cnn'] == r['cnn'] else -1)))
    all_rows[site] = rows

L = []
W = L.append
W('=' * 108)
W('候选清单与质检门禁（GNINA 筛选 × LGR5 胞外域；各路规模见下方逐区域表）')
W('  过滤条件：最近接触 >= %.1f Å  且 affinity < 0  且 CNNscore >= %.2f  且 距区域中心 <= %.1f Å'
  % (MIN_CONTACT, MIN_CNN, SITE_CUT))
W('=' * 108)
W('')
W('%-6s %8s %8s %10s %10s %12s' % ('区域', '分子数', '通过', '通过率', '最佳CNN', '芒柄花素名次'))
for site in SITES:
    rows = all_rows[site]
    passed = [r for r in rows if r['pass']]
    tgt_rank = next((i for i, r in enumerate(rows, 1) if r['id'] == TARGET), None)
    W('%-6s %8d %8d %9.0f%% %10.4f %12s'
      % (site, len(rows), len(passed), 100.0 * len(passed) / len(rows),
         passed[0]['cnn'] if passed else float('nan'), tgt_rank))

for site in SITES:
    W('')
    W('─' * 108)
    W('[区域 %s] 候选清单（前 15）  区域中心 (%.2f, %.2f, %.2f)'
      % (site, *site_centers[site]))
    W('%-5s %-12s %-15s %9s %10s %8s %8s %7s %6s'
      % ('名次', '中文名', '英文名', 'CNNscore', 'affinity', '接触Å', '距区域中心Å', '氢键', '接触'))
    rows = all_rows[site]
    shown = 0
    for r in rows:
        if not r['pass']:
            continue
        mm = meta.get(r['id'], {})
        W('%-5d %-12s %-15s %9.4f %10.2f %8.2f %8.2f %7d %6d'
          % (shown + 1, (mm.get('name_cn') or '')[:10], (mm.get('name_en') or '')[:14],
             r['cnn'], r['aff'], r['dmin'], r['d_site'], r['nhb'], r['ncont']))
        shown += 1
        if shown >= 15:
            break
    # 推荐候选（config.yaml 的 recommended_molecule）那一行的汇总
    tgt = next((r for r in rows if r['id'] == TARGET), None)
    if tgt:
        rk = next((i for i, r in enumerate(rows, 1) if r['id'] == TARGET), None)
        W('  → 芒柄花素：总名次 %s / %d ；CNNscore %.4f ；affinity %.2f ；接触 %.2f Å ；距区域中心 %.2f Å ；判定 %s（%s）'
          % (rk, len(rows), tgt['cnn'], tgt['aff'], tgt['dmin'], tgt['d_site'],
             'PASS' if tgt['pass'] else 'FAIL', tgt['reason']))

# ---- 质检门禁 ----
W('')
W('=' * 108)
W('质检门禁 · 各区域的几何统计（本靶点无共晶小分子，故以几何门禁替代重对接黄金对照）')
W('=' * 108)
W('%-6s %14s %12s %12s %16s' % ('区域', '盒内姿态占比(均值)', '首姿态穿透数', 'CNN-接触相关r', '首姿态入盒分子数'))
for site in SITES:
    mols = sites_data[site]
    cen = site_centers[site]
    occ = []
    for mid, poses in mols.items():
        n = sum(1 for p in poses if math.dist(p['cen'], cen) <= SITE_CUT)
        occ.append(n / len(poses))
    mean_occ = sum(occ) / len(occ)
    clash = sum(1 for mid, poses in mols.items() if poses[0]['dmin'] < MIN_CONTACT)
    xs = [(p[0]['cnn'], p[0]['dmin']) for p in mols.values() if p[0]['cnn'] == p[0]['cnn']]
    if len(xs) >= 10:
        mx = sum(a for a, _ in xs) / len(xs)
        my = sum(b for _, b in xs) / len(xs)
        num = sum((a - mx) * (b - my) for a, b in xs)
        dx = math.sqrt(sum((a - mx) ** 2 for a, _ in xs))
        dy = math.sqrt(sum((b - my) ** 2 for _, b in xs))
        r = num / (dx * dy) if dx and dy else float('nan')
    else:
        r = float('nan')
    correct = sum(1 for mid, poses in mols.items() if math.dist(poses[0]['cen'], cen) <= SITE_CUT)
    W('%-6s %13.0f%% %12d %12.3f %10d/%-3d'
      % (site, 100 * mean_occ, clash, r, correct, len(mols)))
W('')
W('注：上表四项均为【几何统计量】，不含任何实验对照 —— ')
W('    · 盒内姿态占比 / 首姿态入盒分子数：随随机种子与配体起跑位置明显变化，【不作性能指标】；')
W('    · 首姿态穿透数：首姿态与受体最近接触 < 2.5 A 的分子数（几何合理性检查）；')
W('    · CNN-接触相关r：首姿态 CNNscore 与最小接触距离的 Pearson 相关系数（打分-几何一致性检查）。')

open(OUT + '_report.txt', 'w', encoding='utf-8', newline='').write('\n'.join(L) + '\n')

for site in SITES:
    with open('%s_候选清单_%s.tsv' % (OUT, site), 'w', encoding='utf-8-sig', newline='') as f:
        f.write('名次\t编号\t中文名\t英文名\t来源中药\tCNN打分(0-1越高越好)\t结合能(kcal/mol越负越好)\t'
                'CNN预测亲和力\t最近受体接触(A，需≥2.5)\t距位点中心(A，需≤6)\t接触残基数\t氢键残基数\t'
                '判定\t未通过原因\n')
        for i, r in enumerate(all_rows[site], 1):
            mm = meta.get(r['id'], {})
            f.write('%d\t%s\t%s\t%s\t%s\t%.4f\t%.2f\t%.3f\t%.2f\t%.2f\t%d\t%d\t%s\t%s\n'
                    % (i, r['id'], mm.get('name_cn', ''), mm.get('name_en', ''), mm.get('herb', ''),
                       r['cnn'], r['aff'], r['caff'], r['dmin'], r['d_site'], r['ncont'], r['nhb'],
                       'PASS' if r['pass'] else 'FAIL', r['reason']))
    print('写出候选清单: %s_候选清单_%s.tsv（通过 %d / %d）'
          % (OUT, site, sum(1 for r in all_rows[site] if r['pass']), len(all_rows[site])))

print('\n'.join(L[:40]))
print('\n完整报告: %s_report.txt' % OUT)
