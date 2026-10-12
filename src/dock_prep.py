#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
逐文件配体预处理：PubChem 2D SDF -> RDKit 3D SDF（供 gnina 对接）

用法:
    python3 dock_prep.py <src_dir> <dst_dir> <merged_out> <report_tsv> [cx cy cz]
                         [--stereo-policy auto|fix|enumerate|fail|warn|asis]
                         [--stereo-max-isomers N]

说明:
  - src_dir 下每个 *.sdf 视为一个分子（文件名即 CAS 号，作为唯一标识）。
  - 输出 dst_dir/<cas>.sdf（3D、加氢、MMFF 优化），并把 CAS 写入 <CAS> 属性、
    分子名(_Name)也设为 CAS，便于 gnina 输出中回溯。
  - **若给出 cx cy cz**：把每个分子整体平移到该坐标（质心对齐）。
    RDKit 的 EmbedMolecule 生成的构象恒定以原点为中心，不平移的话
    配体会离对接盒子几百埃远（触发 "initial pose not within box"），
    并且 9 个姿态在盒内随机散布、收敛性变差。
  - 同时合并所有 3D 分子为 merged_out（纯文本拼接，保留全部属性）。
  - 生成 report_tsv；给出质心时额外记录质心距口袋的距离并做自检。
  - 随机种子固定为 42，保证可复现。

【立体化学闸门】（2026-09 新增，见 src/stereo_gate.py）
  问题：若分子含**未定义**手性中心（源头记录本身没写立体，如 PubChem 的
        樟脑 CID 2537、α-蒎烯 CID 6654，其 InChIKey 立体层为 UHFFFAOY），
        ETKDG 会**任选一个构型**：固定 seed 下结果确定，但化学上任意。
  判定：不数"潜在中心"，而是**枚举后按唯一 InChI 去重**的异构体数 ——
        自动排除假不对称中心（如 1,8-桉叶油素实为非手性）。
  策略（--stereo-policy）：
        auto(默认) 先尝试用元数据修复（PUBCHEM_SMILES/PUBCHEM_IUPAC_INCHI）；
                   修不了则枚举全部立体异构体（各自成一个配体，文件名 <cas>__isoK）；
                   异构体数超过 --stereo-max-isomers（默认 8）则**明确失败**，交人工决定。
        fix        只尝试用元数据修复，修不了按 fail 处理。
        enumerate  直接枚举（同上限）。
        fail       只要存在歧义就跳过该分子并在报告中标 FAIL。
        warn       照样嵌 3D，但在日志与报告中显式告警。
        asis       嵌 3D 且不告警（仅在需复现历史结果时使用）。
  报告新增列（追加在原 12 列之后，不影响既有解析）：
        lig_id / stereo_policy / n_undefined_centers / n_stereoisomers /
        stereo_source / stereo_note
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import os
import sys
import glob

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Geometry import Point3D
from rdkit.Chem import rdMolTransforms
from rdkit import RDLogger

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import stereo_gate as SG
except ImportError:                     # 找不到模块时退化为旧行为
    SG = None

RDLogger.DisableLog('rdApp.*')

SEED = 42

# ---------------- 解析可选开关（不影响既有位置参数）----------------
STDIO = list(sys.argv[1:])
POLICY = 'auto'
MAX_ISO = 8
rest = []
i = 0
while i < len(STDIO):
    a = STDIO[i]
    if a == '--stereo-policy' and i + 1 < len(STDIO):
        POLICY = STDIO[i + 1]
        i += 2
        continue
    if a.startswith('--stereo-policy='):
        POLICY = a.split('=', 1)[1]
        i += 1
        continue
    if a == '--stereo-max-isomers' and i + 1 < len(STDIO):
        MAX_ISO = int(STDIO[i + 1])
        i += 2
        continue
    if a.startswith('--stereo-max-isomers='):
        MAX_ISO = int(a.split('=', 1)[1])
        i += 1
        continue
    rest.append(a)
    i += 1
if POLICY not in ('auto', 'fix', 'enumerate', 'fail', 'warn', 'asis'):
    print('ERROR: --stereo-policy 只能是 auto/fix/enumerate/fail/warn/asis')
    sys.exit(1)
if SG is None and POLICY != 'asis':
    print('WARN: 未找到 stereo_gate.py，立体闸门关闭（等同 asis）')
    POLICY = 'asis'

src_dir = rest[0] if len(rest) > 0 else 'ligands_2d'
dst_dir = rest[1] if len(rest) > 1 else 'ligands_3d'
merged_out = rest[2] if len(rest) > 2 else 'ligands_3d.sdf'
report_tsv = rest[3] if len(rest) > 3 else 'prep_report.tsv'

POCKET = None
if len(rest) >= 7:
    try:
        POCKET = Point3D(float(rest[4]), float(rest[5]), float(rest[6]))
    except ValueError:
        print('ERROR: 口袋坐标无法解析: %s' % rest[4:7])
        sys.exit(1)

os.makedirs(dst_dir, exist_ok=True)

files = sorted(glob.glob(os.path.join(src_dir, '*.sdf')))
if not files:
    print('ERROR: no *.sdf found in %s' % src_dir)
    sys.exit(1)

print('stereo policy : %s（异构体上限 %d）' % (POLICY, MAX_ISO))

rows = []
ok_list = []
fail_list = []
stereo_flagged = []


def read_props(path):
    """从 SDF 文本里取 SD 属性（闸门修复时要用）"""
    props = {}
    try:
        with open(path, errors='replace') as f:
            lines = f.read().split('\n')
    except Exception:
        return props
    i = 0
    while i < len(lines):
        l = lines[i]
        if l.strip().startswith('> <') or l.startswith('>  <'):
            try:
                k = l.split('<')[1].split('>')[0]
            except IndexError:
                i += 1
                continue
            if i + 1 < len(lines):
                props[k] = lines[i + 1].strip()
        i += 1
    return props


def embed_and_write(mol, lig_id, cas):
    """加氢 -> ETKDG 3D -> MMFF -> 平移 -> 写盘；返回 (out_path, info) 或抛异常"""
    molh = Chem.AddHs(mol)
    # ⚠ 枚举/复制出来的分子可能带着输入的 2D 构象，会与手性约束冲突导致 ETKDG 失败
    molh.RemoveAllConformers()
    note = ''
    ps = AllChem.ETKDGv3()
    ps.randomSeed = SEED
    if AllChem.EmbedMolecule(molh, ps) != 0:
        ps2 = AllChem.ETKDGv3()
        ps2.randomSeed = SEED
        ps2.useRandomCoords = True
        if AllChem.EmbedMolecule(molh, ps2) != 0:
            raise ValueError('embed-failed')
        note = 'random-coords'
    try:
        AllChem.MMFFOptimizeMolecule(molh, maxIters=500)
    except Exception:
        note = (note + ';mmff-skipped').strip(';')

    conf = molh.GetConformer()
    if POCKET is not None:
        cent = rdMolTransforms.ComputeCentroid(conf)
        shift = Point3D(POCKET.x - cent.x, POCKET.y - cent.y, POCKET.z - cent.z)
        for i in range(molh.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Point3D(p.x + shift.x, p.y + shift.y, p.z + shift.z))

    conf = molh.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(molh.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(molh.GetNumAtoms())]
    zs = [conf.GetAtomPosition(i).z for i in range(molh.GetNumAtoms())]
    span = (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
    max_span = max(span)
    natoms = molh.GetNumAtoms()

    d_pocket = float('nan')
    if POCKET is not None:
        cent = rdMolTransforms.ComputeCentroid(conf)
        d_pocket = ((cent.x - POCKET.x) ** 2 + (cent.y - POCKET.y) ** 2
                    + (cent.z - POCKET.z) ** 2) ** 0.5

    if max_span < 0.5:
        raise ValueError('degenerate-3d')

    molh.SetProp('_Name', lig_id)
    molh.SetProp('CAS', cas)
    molh.SetProp('LIG_ID', lig_id)
    out_path = os.path.join(dst_dir, lig_id + '.sdf')
    w = Chem.SDWriter(out_path)
    w.write(molh)
    w.close()
    return out_path, {'note': note, 'span': span, 'max_span': max_span,
                      'natoms': natoms, 'd_pocket': d_pocket}


for path in files:
    cas = os.path.splitext(os.path.basename(path))[0]
    status = 'ok'
    note = ''
    heavy = 0
    formula = ''
    mw = 0.0
    span = (0.0, 0.0, 0.0)
    max_span = 0.0
    natoms = 0
    d_pocket = float('nan')
    n_undef = 0
    n_iso = 1
    stereo_source = 'input'
    s_note = ''

    try:
        supp = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
        mol = supp[0] if len(supp) else None
        if mol is None:
            supp = Chem.SDMolSupplier(path, removeHs=False, sanitize=False)
            mol = supp[0] if len(supp) else None
            if mol is not None:
                try:
                    mol.UpdatePropertyCache(strict=False)
                    Chem.SanitizeMol(mol)
                    note = 'sanitize-recovered'
                except Exception as e:
                    mol = None
                    note = 'sanitize-failed:%s' % str(e)[:60]
            else:
                note = 'unreadable'

        if mol is None:
            raise ValueError(note or 'mol is None')

        heavy = mol.GetNumHeavyAtoms()
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)

        # ---------------- 立体化学闸门 ----------------
        mols_to_write = [mol]
        if SG is not None and POLICY != 'asis':
            st = SG.analyze(mol)
            n_undef = st['n_undefined']
            n_iso = st['n_isomers']
            if st['ambiguous']:
                stereo_flagged.append((cas, st['n_undefined'], st['n_isomers']))
                if POLICY in ('auto', 'fix', 'enumerate'):
                    fixed, how = SG.fix_from_props(mol, read_props(path), title=cas)
                    if fixed is not None:
                        st2 = SG.analyze(fixed)
                        if not st2['ambiguous']:
                            mols_to_write = [fixed]
                            stereo_source = 'metadata:' + how
                            s_note = '由元数据修复立体'
                        elif POLICY == 'fix':
                            mols_to_write = []
                            stereo_source = 'unfixable'
                            s_note = '元数据仍不明确: %s' % how
                        else:
                            iso_list, why = SG.enumerate_isomers(mol, MAX_ISO, title=cas)
                            if not iso_list:
                                mols_to_write = []
                                stereo_source = 'unfixable'
                                s_note = why
                            else:
                                mols_to_write = iso_list
                                stereo_source = 'enumerated'
                                n_iso = len(iso_list)
                                s_note = why
                    elif POLICY == 'fix':
                        mols_to_write = []
                        stereo_source = 'unfixable'
                        s_note = how
                    else:
                        iso_list, why = SG.enumerate_isomers(mol, MAX_ISO, title=cas)
                        if not iso_list:
                            mols_to_write = []
                            stereo_source = 'unfixable'
                            s_note = why
                        else:
                            mols_to_write = iso_list
                            stereo_source = 'enumerated'
                            n_iso = len(iso_list)
                            s_note = why
                elif POLICY == 'fail':
                    mols_to_write = []
                    stereo_source = 'undefined'
                    s_note = '立体未定义且策略=fail'
                else:                       # warn
                    stereo_source = 'picked-by-ETKDG'
                    s_note = ('⚠ 立体未定义（%d 个中心 / %d 种异构体），'
                              'ETKDG 将任选一个构型' % (st['n_undefined'], st['n_isomers']))
            else:
                stereo_source = 'input'
                if st['n_undefined'] or st['n_db_any']:
                    s_note = '潜在中心经对称性去重后仅 1 种异构体（无歧义）'

        if not mols_to_write:
            raise ValueError('stereo-blocked: %s' % (s_note or '立体未定义'))

        # ---------------- 逐异构体嵌 3D 写盘 ----------------
        # 注意：部分枚举出的异构体在几何上不可能存在（如樟脑的双环[2.2.1]桥头构型冲突），
        #       ETKDG 会拒绝嵌入 —— 这是正确行为，应**跳过该异构体**而不是判定整个分子失败。
        n_written = 0
        skipped_iso = []
        for k, m in enumerate(mols_to_write, 1):
            lig_id = cas if len(mols_to_write) == 1 else '%s__iso%d' % (cas, k)
            try:
                out_path, info = embed_and_write(m, lig_id, cas)
            except Exception as e:
                skipped_iso.append('%s(%s)' % (lig_id, str(e)[:24]))
                continue
            n_written += 1
            ok_list.append((lig_id, out_path))
            extra = s_note if len(mols_to_write) == 1 else \
                ('异构体 %d/%d；%s' % (k, len(mols_to_write), s_note))
            note_k = ';'.join(x for x in (note, info['note'], extra) if x)
            rows.append((cas, 'ok', heavy, info['natoms'], formula, '%.1f' % mw,
                         '%.1f' % info['span'][0], '%.1f' % info['span'][1],
                         '%.1f' % info['span'][2], '%.1f' % info['max_span'],
                         '%.2f' % info['d_pocket'] if info['d_pocket'] == info['d_pocket'] else '',
                         note_k, lig_id, POLICY, n_undef, n_iso, stereo_source,
                         s_note if len(mols_to_write) == 1 else
                         ('异构体 %d/%d' % (k, len(mols_to_write)))))
        if n_written == 0:
            raise ValueError('all-isomers-failed: %s' % ';'.join(skipped_iso)[:90])
        if skipped_iso:
            print('  注：%s 有 %d 个枚举异构体无法嵌入（几何不可能），已跳过：%s'
                  % (cas, len(skipped_iso), ', '.join(skipped_iso)))
        continue
    except Exception as e:
        status = 'FAIL'
        note = (note + '|' + str(e))[:120]
        fail_list.append((cas, str(e)[:80]))

    rows.append((cas, status, heavy, natoms, formula, '%.1f' % mw,
                 '%.1f' % span[0], '%.1f' % span[1], '%.1f' % span[2],
                 '%.1f' % max_span, '%.2f' % d_pocket if d_pocket == d_pocket else '',
                 note, cas, POLICY, n_undef, n_iso, stereo_source, s_note))

# ---- 合并 ----
merged_ok = True
try:
    with open(merged_out, 'w') as fo:
        for lig_id, p in ok_list:
            with open(p, 'r') as fi:
                txt = fi.read()
            if not txt.rstrip().endswith('$$$$'):
                txt = txt.rstrip() + '\n$$$$\n'
            fo.write(txt)
except Exception as e:
    merged_ok = False
    print('ERROR merging: %s' % e)

# ---- 报告 ----
with open(report_tsv, 'w') as f:
    f.write('cas\tstatus\theavy_atoms\ttotal_atoms_with_H\tformula\tMW\t'
            'span_x\tspan_y\tspan_z\tmax_span\td_centroid_to_pocket\tnote\t'
            'lig_id\tstereo_policy\tn_undefined_centers\tn_stereoisomers\t'
            'stereo_source\tstereo_note\n')
    for r in rows:
        f.write('\t'.join(str(x) for x in r) + '\n')

print('total files : %d' % len(files))
print('3D ok       : %d  输出配体数: %d' % (len({r[0] for r in rows if r[1] == 'ok'}),
                                           len(ok_list)))
print('failed      : %d' % len(fail_list))
for cas, err in fail_list:
    print('  FAIL %s : %s' % (cas, err))

# ---- 立体化学小结 ----
if POLICY != 'asis':
    # 行结构: 0 cas,1 status,2 heavy,3 natoms,4 formula,5 MW,6-8 span,9 max_span,
    #         10 d_pocket,11 note,12 lig_id,13 policy,14 n_undef,15 n_iso,
    #         16 stereo_source,17 stereo_note
    IDX_SRC = 16
    n_amb = len({c for c, _, _ in stereo_flagged})
    n_enum = len({r[0] for r in rows if r[1] == 'ok' and r[IDX_SRC] == 'enumerated'})
    n_fix = len({r[0] for r in rows if r[1] == 'ok' and str(r[IDX_SRC]).startswith('metadata')})
    n_warn = len({r[0] for r in rows if r[1] == 'ok' and r[IDX_SRC] == 'picked-by-ETKDG'})
    n_blocked = len({r[0] for r in rows if r[1] != 'ok'
                     and r[IDX_SRC] in ('unfixable', 'undefined')})
    n_iso_out = len([r for r in rows if r[1] == 'ok' and '__iso' in str(r[12])])
    print()
    print('--- 立体化学小结（策略 %s）---' % POLICY)
    print('检出立体歧义 : %d 个分子' % n_amb)
    if n_amb:
        for cas, nu, ni in sorted(stereo_flagged, key=lambda x: -x[2])[:12]:
            print('    %-14s 未定义中心 %2d → %d 种立体异构体' % (cas, nu, ni))
    print('由元数据修复 : %d ；枚举展开成 %d 个分子（额外产出 %d 个异构体配体）；'
          '仍任选(告警) : %d ；被拦截 : %d'
          % (n_fix, n_enum, n_iso_out, n_warn, n_blocked))
    if n_blocked:
        print('  ⚠ 被拦截的分子未参与对接（请人工指定构型后单独处理）')

if merged_ok and ok_list:
    spans = [float(r[9]) for r in rows if r[1] == 'ok']
    print('merged file : %s (%d mols)' % (merged_out, len(ok_list)))
    print('max span    : %.1f A  (box must exceed this)' % max(spans))

if POCKET is not None:
    ds = [float(r[10]) for r in rows if r[1] == 'ok' and str(r[10]) != '']
    if ds:
        print('pocket      : (%.3f, %.3f, %.3f)' % (POCKET.x, POCKET.y, POCKET.z))
        print('SELF-CHECK  : 质心距口袋 最近 %.2f / 中位 %.2f / 最远 %.2f A'
              % (min(ds), sorted(ds)[len(ds) // 2], max(ds)))
        if max(ds) < 2.0:
            print('SELF-CHECK  : PASS (全部 < 2 A，配体已就位)')
        else:
            print('SELF-CHECK  : FAIL (仍有配体未进入盒子)')
            sys.exit(3)
else:
    print('pocket      : 未指定（未做平移）')

print('report      : %s' % report_tsv)
