#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stereo_gate.py —— 配体立体化学闸门（供 dock_prep.py 调用）

背景：
    若分子含**未定义**手性中心（源头记录未写立体，InChIKey 立体层为 UHFFFAOY），
    RDKit ETKDG 三维嵌入会**任选一个构型**（固定 seed 下确定，但化学上是任意的），
    且不报错、不告警 —— 对接结果因而建立在任意构型上。

本模块提供：
    analyze(mol)                  -> 立体化学状态（含"是否真有歧义"判定）
    status_text(st)               -> 一行可读描述
    fix_from_props(mol, props)    -> 用权威元数据（PubChem 异构 SMILES / InChI）重建
    enumerate_isomers(mol, max_n) -> 枚举全部立体异构体（去重按 InChI）

设计原则：
    · "有歧义"的判定不靠中心计数，而靠**枚举后唯一 InChI 的个数** ——
      这样自动排除"假不对称中心"（如 1,8-桉叶油素，分子本身有对称面、实为非手性）。
    · 修复优先用**权威来源**（文件自带的 PubChem 异构 SMILES / InChI），
      因为 SDF 里可能只是转换时丢了立体，而元数据仍然完整。
    · 一律不"猜"：修不了就枚举；枚举超出上限就明确失败，交给人工决定。
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
from rdkit import Chem
from rdkit.Chem import EnumerateStereoisomers
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

_ISO_OPTS = EnumerateStereoisomers.StereoEnumerationOptions(
    onlyUnassigned=True, unique=True, tryEmbedding=False)


def _inchi(mol):
    """返回 InChI（用于按对称性正确去重）；失败返回 None"""
    try:
        return Chem.MolToInchi(mol)
    except Exception:
        return None


def _flat_smiles(mol):
    return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=False)),
                            isomericSmiles=False)


def count_defined(mol):
    """已定义 / 未定义 手性中心数、未指定双键立体数"""
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                       useLegacyImplementation=False)
    n_def = sum(1 for _, c in centers if c != '?')
    n_und = sum(1 for _, c in centers if c == '?')
    db = sum(1 for b in mol.GetBonds() if b.GetStereo() == Chem.BondStereo.STEREOANY)
    return n_def, n_und, db


def real_centers(mol):
    """
    **真正**的潜在手性中心（旧实现 force=True）。
    与 includeUnassigned 的"潜在中心"列表不同，它会剔除**假不对称中心**：
    例：1,8-桉叶油素（470-82-6）有对称面、实为非手性，新实现会报 2 个潜在中心，
    本函数返回 []（实测：其两个"异构体"是同一分子的镜像，InChI 仅差 /m 层，
    不能靠 InChI 去重区分）。因此"是否存在立体歧义"必须以本函数为准。
    """
    try:
        return Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                        useLegacyImplementation=True, force=True)
    except Exception:
        return Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                        useLegacyImplementation=False)


def analyze(mol):
    """
    返回 dict:
        n_defined / n_undefined / n_db_any : 计数（潜在中心口径）
        n_real      : 真手性中心（已定义+未定义）
        n_isomers   : 真正不同的立体异构体数（按唯一 InChI 去重）
        ambiguous   : True 表示"存在真实立体歧义"（ETKDG 会任选一个构型）
    """
    n_def, n_und, db = count_defined(mol)
    rc = real_centers(mol)
    st = {'n_defined': n_def, 'n_undefined': n_und, 'n_db_any': db,
          'n_real': len(rc), 'n_isomers': 1, 'ambiguous': False, 'isomers': [mol]}
    # 没有真手性中心、也没有未指定双键 → 无歧义（假不对称中心不算）
    if st['n_real'] == 0 and db == 0:
        return st
    if n_und == 0 and db == 0:
        return st                      # 全部已定义
    try:
        isos = list(EnumerateStereoisomers.EnumerateStereoisomers(mol, _ISO_OPTS))
    except Exception:
        isos = [mol]
    seen, uniq = set(), []
    for m in isos:
        key = _inchi(m) or Chem.MolToSmiles(m, isomericSmiles=True)
        if key not in seen:
            seen.add(key)
            uniq.append(m)
    st['n_isomers'] = max(1, len(uniq))
    st['isomers'] = uniq
    st['ambiguous'] = st['n_isomers'] > 1
    return st


def status_text(st):
    if not st['ambiguous']:
        if st['n_undefined'] or st['n_db_any']:
            return '无不明确立体（%d 个"潜在中心"经对称性去重后只有 1 种异构体）' % st['n_undefined']
        return '立体已完全定义'
    return '⚠ 立体不明确：%d 个未定义手性中心 → 存在 %d 种立体异构体' % (
        st['n_undefined'], st['n_isomers'])


def fix_from_props(mol, props, title=None):
    """
    用权威元数据重建分子（保留原连接性）。props 为 SD 属性字典。
    优先 PUBCHEM_SMILES（异构 SMILES），其次 PUBCHEM_IUPAC_INCHI。
    返回 (new_mol 或 None, 来源说明)
    """
    cands = []
    smi = (props or {}).get('PUBCHEM_SMILES') or (props or {}).get('PUBCHEM_ISOMERIC_SMILES')
    if smi:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            cands.append((m, 'PUBCHEM_SMILES'))
    inchi = (props or {}).get('PUBCHEM_IUPAC_INCHI')
    if inchi:
        try:
            m = Chem.MolFromInchi(inchi)
            if m is not None:
                cands.append((m, 'PUBCHEM_IUPAC_INCHI'))
        except Exception:
            pass

    flat_ref = _flat_smiles(mol)
    best = None
    for m, src in cands:
        if _flat_smiles(m) != flat_ref:
            continue                      # 连接性必须一致，否则不用
        st = analyze(m)
        if not st['ambiguous']:
            best = (m, src, st)
            break
        if best is None:
            best = (m, src, st)           # 至少比原来好（记下来）
    if best is None:
        return None, '元数据无可用立体信息' if not cands else '元数据连接性与输入不一致'
    m, src, st = best
    if title:
        m.SetProp('_Name', title)
    return m, '%s（%d 种异构体）' % (src, st['n_isomers'])


def enumerate_isomers(mol, max_n=8, title=None):
    """
    枚举全部立体异构体。返回 (list_of_mol, note)
    若需要的异构体数超过 max_n，返回 ([], 原因) —— 交由调用方决定失败/跳过。
    """
    st = analyze(mol)
    if not st['ambiguous']:
        return [mol], '立体已定义'
    if st['n_isomers'] > max_n:
        return [], ('需要的立体异构体数 %d 超过上限 %d（未定义中心 %d 个），'
                    '请人工指定构型' % (st['n_isomers'], max_n, st['n_undefined']))
    out = []
    for k, m in enumerate(st['isomers'], 1):
        if title:
            m.SetProp('_Name', '%s__iso%d' % (title, k))
            m.SetProp('STEREO_ISOMER', str(k))
            m.SetProp('STEREO_ISOMER_TOTAL', str(len(st['isomers'])))
        out.append(m)
    return out, '已枚举 %d 种立体异构体' % len(out)
