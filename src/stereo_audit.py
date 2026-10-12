#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stereo_audit.py —— 配体库立体化学体检（只读，不改任何文件）

用法:
    python3 stereo_audit.py <sdf_dir_or_file> [out_tsv]

输出（逐分子）：
    文件名 / 重原子数 / 真手性中心数 / 已定义 / 未定义 / 未指定双键立体 /
    可能立体异构体数 / 判定（完全定义 | 无歧义(假不对称) | ⚠ 存在歧义）
    以及文件自带 PubChem 元数据是否可用于修复

退出码：0 = 全部无歧义；1 = 存在立体歧义分子（便于在流水线里做门禁）

背景：源头 PubChem 记录可能本身就是"立体未指定"（InChIKey 第二段为 UHFFFAOY），
      此时 RDKit ETKDG 会**任选一个构型**。本工具用于在对接前识别这类分子。
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import glob
import os
import re
import sys

from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import stereo_gate as SG
except ImportError:
    SG = None


def read_props(path):
    props = {}
    try:
        lines = open(path, errors='replace').read().split('\n')
    except Exception:
        return props
    for i, l in enumerate(lines):
        if l.strip().startswith('> <'):
            try:
                k = l.split('<')[1].split('>')[0]
            except IndexError:
                continue
            if i + 1 < len(lines):
                props[k] = lines[i + 1].strip()
    return props


def inchikey_stereo_flag(ik):
    """InChIKey 第二段前 8 位 == UHFFFAOY 表示该记录无立体层"""
    if not ik or ik.count('-') < 2:
        return None
    return '无立体' if ik.split('-')[1][:8] == 'UHFFFAOY' else '有立体'


target = sys.argv[1] if len(sys.argv) > 1 else 'ligands_2d'
out_tsv = sys.argv[2] if len(sys.argv) > 2 else 'stereo_audit.tsv'

files = sorted(glob.glob(os.path.join(target, '*.sdf'))) if os.path.isdir(target) else [target]
if not files:
    print('ERROR: 未找到 sdf 文件: %s' % target)
    sys.exit(2)

rows = []
for p in files:
    name = os.path.splitext(os.path.basename(p))[0]
    m = Chem.MolFromMolFile(p, removeHs=False, sanitize=True)
    if m is None:
        rows.append({'name': name, 'err': 'parse-fail'})
        continue
    st = SG.analyze(m) if SG else {'n_defined': 0, 'n_undefined': 0, 'n_db_any': 0,
                                   'n_isomers': 1, 'ambiguous': False, 'n_real': 0}
    props = read_props(p)
    ik = props.get('PUBCHEM_IUPAC_INCHIKEY', '')
    fixable = False
    if st['ambiguous'] and SG is not None:
        fixed, _ = SG.fix_from_props(m, props)
        if fixed is not None and not SG.analyze(fixed)['ambiguous']:
            fixable = True
    verdict = ('完全定义' if not st['ambiguous'] and not (st['n_undefined'] or st['n_db_any'])
               else '无歧义(假不对称中心)' if not st['ambiguous']
               else '⚠ 存在立体歧义')
    rows.append({'name': name, 'heavy': m.GetNumHeavyAtoms(),
                 'n_real': st.get('n_real', 0), 'n_def': st['n_defined'],
                 'n_und': st['n_undefined'], 'n_db': st['n_db_any'],
                 'n_iso': st['n_isomers'], 'verdict': verdict,
                 'ik_flag': inchikey_stereo_flag(ik) or '-',
                 'fixable': '是' if fixable else ('否' if st['ambiguous'] else '-')})

with open(out_tsv, 'w', encoding='utf-8') as f:
    f.write('文件\t重原子数\t真手性中心\t已定义中心\t未定义中心\t未指定双键\t可能异构体数\t判定\t'
            'PubChem记录立体层\t元数据可修复\n')
    for r in rows:
        if 'err' in r:
            f.write('%s\t%s\n' % (r['name'], r['err']))
            continue
        f.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n'
                % (r['name'], r['heavy'], r['n_real'], r['n_def'], r['n_und'], r['n_db'],
                   r['n_iso'], r['verdict'], r['ik_flag'], r['fixable']))

amb = [r for r in rows if 'err' not in r and r['verdict'].startswith('⚠')]
fake = [r for r in rows if 'err' not in r and r['verdict'].startswith('无歧义')]
print('=== 立体化学体检: %s ===' % target)
print('  分子数 %d ；解析失败 %d' % (len(rows), sum(1 for r in rows if 'err' in r)))
print('  ✅ 完全定义     : %d' % sum(1 for r in rows if r.get('verdict') == '完全定义'))
print('  ○ 无歧义(假不对称): %d' % len(fake))
print('  ⚠ 存在立体歧义  : %d' % len(amb))
if amb:
    print()
    print('  %-14s %-8s %-10s %-12s %-10s %s'
          % ('文件', '未定义中心', '可能异构体', 'PubChem立体层', '可修复', '备注'))
    for r in sorted(amb, key=lambda x: -x['n_iso']):
        note = '异构体数过多，建议人工指定' if r['n_iso'] > 8 else \
               ('可由元数据修复' if r['fixable'] == '是' else '需枚举或人工指定')
        print('  %-14s %-8d %-10d %-12s %-10s %s'
              % (r['name'], r['n_und'], r['n_iso'], r['ik_flag'], r['fixable'], note))
print()
print('已写出: %s' % out_tsv)
print('建议: 存在歧义的分子在对接前应经 dock_prep.py 的立体闸门处理'
      '（默认 --stereo-policy auto：先修元数据、再枚举、超限则拦截）')
sys.exit(1 if amb else 0)
