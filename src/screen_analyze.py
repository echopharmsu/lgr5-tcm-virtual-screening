#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
全库筛选结果分析：排名、跨盒共识、与原 Vina 分数对照。

用法: python3 screen_analyze.py <ligands_source.tsv> <receptor.pdb> <out_prefix> <screen_U1.sdf.gz> <screen_U3.sdf.gz> <screen_U10.sdf.gz>
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
TARGET = '001'          # 刺芒柄花素
CONTACT, CLASH, HB = 4.0, 2.5, 3.5
POLAR = {'N', 'O', 'S'}

meta = {}
with open(SRC_TSV, encoding='utf-8') as f:
    hdr = f.readline().rstrip('\n').split('\t')
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= len(hdr):
            meta[p[0]] = dict(zip(hdr, p))

# 分子尺寸（用于偏倚检查：浅表面上大分子接触面积大、分数易虚高）
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    for m, v in meta.items():
        mol = Chem.MolFromSmiles(v.get('smiles', ''))
        if mol is not None:
            v['mw'] = Descriptors.MolWt(mol)
            v['ha'] = mol.GetNumHeavyAtoms()
        else:
            v['mw'], v['ha'] = float('nan'), 0
except ImportError:
    for v in meta.values():
        v['mw'], v['ha'] = float('nan'), 0

rec = []
for l in open(REC, errors='replace'):
    if not l.startswith('ATOM'):
        continue
    el = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
    if el in ('H', 'D'):
        continue
    rec.append((l[17:20].strip(), l[22:26].strip(), el,
                (float(l[30:38]), float(l[38:46]), float(l[46:54]))))


def parse_result(path):
    op = gzip.open if path.endswith('.gz') else open
    txt = op(path, 'rt', errors='replace').read()
    mols = {}
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
                    c[0] = min(c[0], d)
                    c[1] += 1
                    if el in POLAR and rel in POLAR and d <= HB:
                        c[2] += 1
        cen = (sum(x for _, (x, _, _) in lig) / len(lig),
               sum(y for _, (_, y, _) in lig) / len(lig),
               sum(z for _, (_, _, z) in lig) / len(lig))
        rec_ = {
            'id': mid,
            'cnn': float(props.get('CNNscore', 'nan')),
            'caff': float(props.get('CNNaffinity', 'nan')),
            'aff': float(props.get('minimizedAffinity', 'nan')),
            'ncont': len(contacts), 'dmin': dmin,
            'density': sum(c[1] for c in contacts.values()) / len(lig),
            'nhb': sum(1 for c in contacts.values() if c[2]),
            'clash': dmin < CLASH,
            'contacts': sorted(contacts.items(), key=lambda kv: kv[1][0]),
            'cen': cen, 'natoms': len(lig),
        }
        # gnina 按 CNNscore 降序输出：只保留每分子首个（最优）姿态
        if mid not in mols:
            mols[mid] = rec_
    return mols


data = {}
for path, site in zip(RESULTS, SITES):
    data[site] = parse_result(path)
    print('%-5s 分子数 %d' % (site, len(data[site])))

# ---------------- 每位点排名 ----------------
def rank_by(site, key, reverse):
    items = [(m, v[key]) for m, v in data[site].items() if v[key] == v[key]]
    items.sort(key=lambda kv: -kv[1] if reverse else kv[1])
    return {m: i for i, (m, _) in enumerate(items, 1)}, len(items)


L = []
W = L.append
W('=' * 112)
W('全库筛选：%d 个中药小分子 × %d 个对接盒' % (len(meta), len(SITES)))
W('  协议：24 Å 盒 · exhaustiveness 8 · num_modes 9 · seed 42 · 受体为 PDB 4KNG 三元复合物抽链 A（非游离态）')
W('=' * 112)
W('')
W('--- 芒柄花素（001 刺芒柄花素 Formononetin，原 Vina −7.3 / 综合排名 9）在各轮的表现 ---')
W('%-6s %10s %10s %12s %10s %8s %8s %10s %10s'
  % ('区域', 'CNNscore', 'CNN名次', 'affinity', 'aff名次', '接触', '氢键', '接触密度', '判定'))
summary_rows = []
for site in SITES:
    r_cnn, n = rank_by(site, 'cnn', True)
    r_aff, _ = rank_by(site, 'aff', False)
    t = data[site].get(TARGET)
    if t is None:
        W('%-6s %s' % (site, '推荐候选无结果'))
        continue
    verdict = 'PASS' if (not t['clash'] and t['aff'] < 0) else 'FAIL'
    W('%-6s %10.4f %6s/%-3d %12.2f %6s/%-3d %8d %8d %10.2f %10s'
      % (site, t['cnn'], r_cnn.get(TARGET), n, t['aff'], r_aff.get(TARGET), n,
         t['ncont'], t['nhb'], t['density'], verdict))
    summary_rows.append((site, t, r_cnn.get(TARGET), n, r_aff.get(TARGET)))

W('')
for site in SITES:
    W('─' * 112)
    W('[区域 %s] Top 12（按 CNNscore）' % site)
    W('%-5s %-14s %-16s %9s %10s %7s %6s %8s %7s'
      % ('ID', '中文名', '英文名', 'CNNscore', 'affinity', 'MW', '接触', '氢键', '穿透'))
    order = sorted([m for m in data[site] if data[site][m]['cnn'] == data[site][m]['cnn']],
                   key=lambda m: -data[site][m]['cnn'])[:12]
    for m in order:
        v = data[site][m]
        mm = meta.get(m, {})
        mark = '  ← 芒柄花素' if m == TARGET else ''
        W('%-5s %-14s %-16s %9.4f %10.2f %7.1f %6d %8d %7s%s'
          % (m, (mm.get('name_cn') or '')[:12], (mm.get('name_en') or '')[:15],
             v['cnn'], v['aff'], mm.get('mw', float('nan')), v['ncont'], v['nhb'],
             'Y' if v['clash'] else '-', mark))

# ---------------- 尺寸偏倚与汇总 ----------------
W('')
W('=' * 112)
W('跨盒共识（各盒都进入前 20 的分子）')
W('=' * 112)
ranks = {}
for site in SITES:
    r_cnn, n = rank_by(site, 'cnn', True)
    for m, rk in r_cnn.items():
        ranks.setdefault(m, {})[site] = rk
consensus = [(m, d) for m, d in ranks.items() if all(d.get(s, 999) <= 20 for s in SITES)]
consensus.sort(key=lambda kv: sum(kv[1].values()))
W('%-5s %-14s %-16s %s' % ('ID', '中文名', '英文名',
                            '  '.join('%s名次' % s for s in SITES) + '   名次和'))
for m, d in consensus[:20]:
    mm = meta.get(m, {})
    mark = '  ← 芒柄花素' if m == TARGET else ''
    W('%-5s %-14s %-16s %s    %d%s'
      % (m, (mm.get('name_cn') or '')[:12], (mm.get('name_en') or '')[:15],
         '  '.join('%7d' % d[s] for s in SITES), sum(d.values()), mark))
if not consensus:
    W('  无：各盒没有共同的高分分子')

# ---------------- 与原 Vina 分数对照 ----------------
W('')
W('=' * 112)
W('与背景分数的对照（仅作参照，不参与判定；仅针对各盒都跑到的分子）')
W('=' * 112)


def spearman(a, b):
    if len(a) < 3 or len(a) != len(b):
        return float('nan')

    def rk(xs):
        order = sorted(range(len(xs)), key=lambda i: xs[i])
        r = [0.0] * len(xs)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    ra, rb = rk(a), rk(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da = math.sqrt(sum((x - ma) ** 2 for x in ra))
    db = math.sqrt(sum((y - mb) ** 2 for y in rb))
    return num / (da * db) if da and db else float('nan')


for site in SITES:
    xs, ys_cnn, ys_aff, ids = [], [], [], []
    for m, v in data[site].items():
        vo = meta.get(m, {}).get('vina_orig', '')
        try:
            vo = float(vo)
        except (TypeError, ValueError):
            continue
        if v['aff'] != v['aff']:
            continue
        xs.append(vo)
        ys_cnn.append(v['cnn'])
        ys_aff.append(v['aff'])
        ids.append(m)
    if len(xs) < 10:
        W('  %-5s 有效配对数不足（%d）' % (site, len(xs)))
        continue
    W('  %-5s n=%3d   Spearman(Vina原分 vs gnina affinity) = %+.3f   vs CNNscore = %+.3f'
      % (site, len(xs), spearman(xs, ys_aff), spearman(xs, ys_cnn)))

W('')
W('ρ 为正值 = 两套打分的排序方向一致（Vina 越负 ↔ gnina affinity 越负 / CNN 越高）；')
W('两套打分都是“越负/越低越好”的同一朝向，因此 ρ 越接近 +1 越一致；|ρ|<0.2 才算无关。')

# ---------------- 尺寸偏倚检查 ----------------
W('')
W('=' * 112)
W('尺寸偏倚检查（浅表面上大分子接触面积大，打分容易虚高）')
W('=' * 112)
for site in SITES:
    mws, cnns, affs = [], [], []
    for m, v in data[site].items():
        mw = meta.get(m, {}).get('mw', float('nan'))
        if mw != mw or v['cnn'] != v['cnn'] or v['aff'] != v['aff']:
            continue
        mws.append(mw); cnns.append(v['cnn']); affs.append(v['aff'])
    W('  %-5s n=%3d  Spearman(MW vs CNNscore) = %+.3f   Spearman(MW vs affinity) = %+.3f'
      % (site, len(mws), spearman(mws, cnns), spearman(mws, affs)))
    sub = [(m, v) for m, v in data[site].items()
           if 200 <= meta.get(m, {}).get('mw', 0) <= 450 and v['cnn'] == v['cnn']]
    sub.sort(key=lambda kv: -kv[1]['cnn'])
    rk = next((i for i, (m, _) in enumerate(sub, 1) if m == TARGET), None)
    W('       尺寸匹配子集（MW 200–450，n=%d）内芒柄花素 CNN 名次: %s / %d'
      % (len(sub), rk, len(sub)))

open(OUT + '.txt', 'w', encoding='utf-8', newline='').write('\n'.join(L) + '\n')

with open(OUT + '_per_site.tsv', 'w', encoding='utf-8', newline='') as f:
    f.write('位点\tID\t中文名\t英文名\t来源中药\tMW\t重原子数\tCNNscore\tCNN名次\tCNNaffinity\taffinity\t'
            'aff名次\t接触残基数\t氢键残基数\t接触密度\t最近接触A\t穿透\t原Vina分数\t'
            '原综合评分\t原综合排名\tx\ty\tz\n')
    for site in SITES:
        r_cnn, n = rank_by(site, 'cnn', True)
        r_aff, _ = rank_by(site, 'aff', False)
        for m, v in sorted(data[site].items(), key=lambda kv: kv[1]['cnn'], reverse=True):
            mm = meta.get(m, {})
            f.write('%s\t%s\t%s\t%s\t%s\t%.1f\t%d\t%.4f\t%d\t%.3f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f\t%s\t%s\t%s\t%s\t'
                    '%.2f\t%.2f\t%.2f\n'
                    % (site, m, mm.get('name_cn', ''), mm.get('name_en', ''), mm.get('herb', ''),
                       mm.get('mw', float('nan')), mm.get('ha', 0),
                       v['cnn'], r_cnn.get(m, 0), v['caff'], v['aff'], r_aff.get(m, 0),
                       v['ncont'], v['nhb'], v['density'], v['dmin'],
                       'YES' if v['clash'] else 'no', mm.get('vina_orig', ''),
                       mm.get('score_total', ''), mm.get('rank_orig', ''),
                       v['cen'][0], v['cen'][1], v['cen'][2]))

# 汇总表：每分子在各盒的名次
with open(OUT + '_consensus.tsv', 'w', encoding='utf-8', newline='') as f:
    f.write('ID\t中文名\t英文名\t来源中药\t' + '\t'.join('%s_CNN名次' % s for s in SITES)
            + '\t名次和\t' + '\t'.join('%s_CNNscore' % s for s in SITES)
            + '\t原Vina分数\t原综合排名\n')
    rows = []
    for m in meta:
        d = ranks.get(m, {})
        if not d:
            continue
        rows.append((m, sum(d.get(s, 999) for s in SITES), d))
    rows.sort(key=lambda r: r[1])
    for m, s, d in rows:
        mm = meta.get(m, {})
        f.write('%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n'
                % (m, mm.get('name_cn', ''), mm.get('name_en', ''), mm.get('herb', ''),
                   '\t'.join(str(d.get(s, '')) for s in SITES), s,
                   '\t'.join('%.4f' % data[s][m]['cnn'] if m in data[s] else ''
                             for s in SITES),
                   mm.get('vina_orig', ''), mm.get('rank_orig', '')))

print('\n'.join(L))
print('已写出: %s.txt / %s_per_site.tsv / %s_consensus.tsv' % (OUT, OUT, OUT))
