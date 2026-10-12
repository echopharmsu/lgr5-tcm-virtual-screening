#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""为盲对接准备受体：4KNG 抽链 + altLoc 去重 + 计算覆盖全蛋白的盒子。"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import os
from collections import Counter, defaultdict

SRC = path('data', '4KNG_full.pdb')
OUT = path('data')
os.makedirs(OUT, exist_ok=True)
PAD = 6.0


def build(chains, out_name):
    kept, dropped_alt, het_skip, alt_seen = [], 0, Counter(), defaultdict(set)
    for line in open(SRC, errors='replace'):
        rec = line[:6].strip()
        if rec not in ('ATOM', 'HETATM'):
            continue
        if line[21] not in chains:
            continue
        if rec == 'HETATM':
            het_skip[line[17:20].strip()] += 1
            continue
        alt = line[16]
        if alt not in (' ', 'A'):
            dropped_alt += 1
            alt_seen['%s%s' % (line[17:20].strip(), line[22:26].strip())].add(alt)
            continue
        kept.append(line.rstrip('\n'))

    xs = [float(l[30:38]) for l in kept]
    ys = [float(l[38:46]) for l in kept]
    zs = [float(l[46:54]) for l in kept]
    cx, cy, cz = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2
    sx, sy, sz = max(xs) - min(xs) + PAD, max(ys) - min(ys) + PAD, max(zs) - min(zs) + PAD

    p = os.path.join(OUT, out_name)
    with open(p, 'w') as f:
        f.write('\n'.join(kept) + '\nEND\n')

    res = {int(l[22:26]) for l in kept}
    print('--- %s ---' % out_name)
    print('   链 %s | ATOM %d | 残基 %d | 剔除 altLoc %d（%s）| 剔除 HETATM %s'
          % (''.join(sorted(chains)), len(kept), len(res), dropped_alt,
             dict((k, ''.join(sorted(v))) for k, v in alt_seen.items()) or '无', dict(het_skip) or '无'))
    print('   蛋白跨度 %.1f × %.1f × %.1f Å' % (sx - PAD, sy - PAD, sz - PAD))
    print('   盒子中心 (%.3f, %.3f, %.3f)  尺寸 %.1f × %.1f × %.1f Å'
          % (cx, cy, cz, sx, sy, sz))
    print('   体积 %.0f Å³' % (sx * sy * sz))
    return (cx, cy, cz, sx, sy, sz)


b_apo = build({'A'}, 'receptor_apo_LGR5.pdb')
b_holo = build({'A', 'M'}, 'receptor_holo_LGR5_RSPO1.pdb')

with open(os.path.join(OUT, 'boxes.txt'), 'w') as f:
    # 注：文件名沿用历史惯例；此处的「apo」仅表示「只保留 LGR5 A 链」，
    #     该构象取自 4KNG（LGR5–RSPO1–RNF43 三元复合物），并非游离态结构。
    f.write('# A = 仅 LGR5 A 链；A+M = LGR5 A + RSPO1 M（均取自 4KNG 三元复合物）\n')
    f.write('apo  %.3f %.3f %.3f %.1f %.1f %.1f\n' % b_apo)
    f.write('holo %.3f %.3f %.3f %.1f %.1f %.1f\n' % b_holo)
print('\n盒子参数已写出: %s' % os.path.join(OUT, 'boxes.txt'))
