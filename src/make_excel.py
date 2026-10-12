#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
把筛选结果汇总成多工作表 Excel。
用法: python -X utf8 make_excel.py <results_dir> <out_xlsx>
"""
from config import ROOT, HERE, path, CFG   # 路径全部来自 config.yaml，禁止硬编码
import csv
import os
import sys

import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

RES = sys.argv[1]
OUT = sys.argv[2]
SITES = ['U1', 'U3', 'U10']

HDR = PatternFill('solid', fgColor='2F5597')
HDRF = Font(bold=True, color='FFFFFF', size=10)
GOLD = PatternFill('solid', fgColor='FFE699')
OKF = PatternFill('solid', fgColor='E2EFDA')
BADF = PatternFill('solid', fgColor='FCE4E4')
INFO = PatternFill('solid', fgColor='DDEBF7')
THIN = Side(style='thin', color='BFBFBF')
BD = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)


def read_tsv(p, enc='utf-8-sig'):
    if not os.path.exists(p):
        return []
    return list(csv.DictReader(open(p, encoding=enc, errors='replace'), delimiter='\t'))


def fn(s, d=None):
    try:
        return float(s)
    except (TypeError, ValueError):
        return d


def head(ws, widths, left=(), h=32):
    n = len(widths)
    for c in range(1, n + 1):
        cell = ws.cell(row=1, column=c)
        cell.fill = HDR
        cell.font = HDRF
        cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
        cell.border = BD
    for i, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = 'A2'
    ws.row_dimensions[1].height = h
    for r in range(2, ws.max_row + 1):
        for c in range(1, n + 1):
            cell = ws.cell(row=r, column=c)
            cell.border = BD
            cell.font = Font(size=10)
            cell.alignment = Alignment(horizontal='left' if c in left else 'center', wrap_text=False)
    ws.auto_filter.ref = 'A1:%s%d' % (get_column_letter(n), ws.max_row)


wb = openpyxl.Workbook()

# ---- 1) 说明 ----
ws = wb.active
ws.title = '说明与参数'
notes = [
    ('LGR5(GPR49) gnina 筛选 —— 109 个中药小分子（v2 · 立体构型修正版）', ''),
    ('', ''),
    ('靶点', 'PDB 4KNG = human LGR5–RSPO1–RNF43；受体仅取 LGR5 A 链（胞外域，残基 30–543，3616 重原子）'),
    ('靶点核对', 'RCSB: entity1 = LGR5, UniProt O75473（别名 GPR49）；本地序列与 O75473 逐位 0 处不一致'),
    ('配体库', '原 110 个中药小分子（综合评分表(1).xlsx）。v2 的两处变更：'
             '① 11 个分子的立体构型改用「修正版 SMILES」并重跑对接'
             '（016/022/041/076/077/079/080/081/082/083/084）；'
             '② 剔除 #89 桔梗皂苷D3（36 个手性中心全部未定义，公开库无法补全）→ 本表为 109 个分子'),
    ('配体 3D', 'RDKit ETKDGv3（seed 42）+ MMFF，逐区域平移到盒中心；立体构型由 SMILES 显式指定'),
    ('立体化学', '⚠ 原 110 个分子中有 9 个在 PubChem 源头记录里立体未指定（按英文名命中了"父记录"），'
             'ETKDG 会任选构型并静默固化；v2 已用修正版构型重跑其中 11 个。'
             '重跑日志见 logs/screen_U*.log'),
    ('筛选协议', '24 Å 盒 · exhaustiveness 8 · num_modes 9 · seed 42 · 受体为 PDB 4KNG 三元复合物抽链 A（非游离态）'),
    ('过滤条件', '最近接触 ≥ 2.5 Å 且 affinity < 0 且 CNNscore ≥ 0.60 且 距区域中心 ≤ 6 Å'),
    ('', ''),
    ('筛选区域', 'LRR13–16 凹面（残基 353–446；LRR 边界取自 UniProt O75473 的 LRR 注释）'),
    ('', ''),
    ('黄金对照', '本靶点无共晶小分子（RCSB 全库 LGR5 结构均无类药配体）→ 无法做重对接 RMSD 验证，已改为【盒内姿态占比 / 首姿态穿透数 / CNN-接触相关r / 首姿态入盒分子数】四项几何统计（均不作性能指标）'),
    ('局限', '① 受体刚性且仅胞外域（不含 7TM）；② 目标区域埋藏比约 0.30，属浅表面结合，非真实口袋；'
             '③ CNNscore 不可跨轮比较；④ 结果性质为假说生成'),
    ('生成脚本', 'predict.py（主运行入口，见 README 第三节）'),
]
ws.append(['项目', '说明'])
for a, b in notes:
    ws.append([a, b])
head(ws, [26, 130], h=28)
for r in range(2, ws.max_row + 1):
    a = ws.cell(row=r, column=1)
    b = ws.cell(row=r, column=2)
    a.font = Font(size=10, bold=str(a.value).startswith(('选中的', 'LGR5')))
    b.font = Font(size=10)
    a.alignment = Alignment(vertical='top')
    b.alignment = Alignment(wrap_text=True, vertical='top')
    if str(a.value).startswith(('筛选协议', '过滤条件', '黄金对照', '局限')):
        a.fill = INFO
        b.fill = INFO
ws.auto_filter.ref = None

# ---- 2) 候选清单（每区域一表）----
for site in SITES:
    rows = read_tsv(os.path.join(RES, 'screen_候选清单_%s.tsv' % site))
    if not rows:
        continue
    w = wb.create_sheet('候选-%s' % site)
    cols = list(rows[0].keys())
    w.append(cols)
    for r in rows:
        w.append([fn(r[c]) if c not in ('编号', '中文名', '英文名', '来源中药', '判定', '未通过原因')
                  else r[c] for c in cols])
    head(w, [6, 6, 14, 17, 10, 11, 12, 11, 11, 11, 9, 9, 8, 26, 10], left=(3, 4, 14))
    for r in range(2, w.max_row + 1):
        is_pass = w.cell(row=r, column=cols.index('判定') + 1).value == 'PASS'
        w.cell(row=r, column=c).fill = OKF if is_pass else BADF

# ---- 3) 跨盒对比 ----
cons = read_tsv(os.path.join(RES, 'screen_unified_consensus.tsv'))
if cons:
    w = wb.create_sheet('跨盒对比')
    cols = list(cons[0].keys())
    w.append(cols)
    for r in cons:
        w.append([r[c] if c in ('ID', '中文名', '英文名', '来源中药') else fn(r[c]) for c in cols])
    head(w, [6, 14, 17, 10, 12, 12, 12, 9, 12, 12, 12, 10, 10], left=(2, 3, 4))
    for r in range(2, w.max_row + 1):
        if str(w.cell(row=r, column=1).value) == '001':
            for c in range(1, len(cols) + 1):
                w.cell(row=r, column=c).fill = GOLD
        elif r <= 11:
            for c in range(1, len(cols) + 1):
                w.cell(row=r, column=c).fill = OKF

# ---- 4) 全部明细（每区域合并）----
w = wb.create_sheet('全部明细')
w.append(['位点', 'ID', '中文名', '英文名', '来源中药', 'CNNscore', 'CNN名次', 'affinity',
          'aff名次', '接触残基数', '氢键残基数', '接触密度', '最近接触A', '穿透',
          '原Vina分数', '原综合评分', '原综合排名'])
for site in SITES:
    for r in read_tsv(os.path.join(RES, 'screen_unified_per_site.tsv')):
        if r.get('位点') != site:
            continue
        w.append([site, r['ID'], r['中文名'], r['英文名'], r['来源中药'],
                  fn(r['CNNscore']), fn(r['CNN名次']), fn(r['affinity']), fn(r['aff名次']),
                  fn(r['接触残基数']), fn(r['氢键残基数']), fn(r['接触密度']), fn(r['最近接触A']),
                  r['穿透'], fn(r['原Vina分数']), fn(r['原综合评分']), fn(r['原综合排名'])])
head(w, [6, 6, 14, 17, 10, 10, 9, 10, 9, 10, 10, 9, 11, 7, 10, 10, 10], left=(3, 4, 5))

# ---- 5) 芒柄花素专项 ----
w = wb.create_sheet('芒柄花素专项')
w.append(['位点', 'CNNscore', 'CNN名次', 'affinity', 'aff名次', '接触残基数', '氢键残基数', '接触密度',
          '最近接触A', '原Vina分数(kcal/mol)', '原综合评分', '原综合排名'])
w.append(['（背景参照，不参与判定）', '', '', '', '', '', '', '', '', -7.3, 77.33, 9])
for r in read_tsv(os.path.join(RES, 'screen_unified_per_site.tsv')):
    if r.get('ID') != '001':
        continue
    w.append([r['位点'], fn(r['CNNscore']), fn(r['CNN名次']), fn(r['affinity']), fn(r['aff名次']),
              fn(r['接触残基数']), fn(r['氢键残基数']), fn(r['接触密度']), fn(r['最近接触A']),
              fn(r['原Vina分数']), fn(r['原综合评分']), fn(r['原综合排名'])])
head(w, [22, 10, 9, 10, 9, 10, 10, 9, 11, 14, 11, 11])
for r in range(2, w.max_row + 1):
    for c in range(1, 13):
        w.cell(row=r, column=c).fill = GOLD if r == 2 else OKF

wb.save(OUT)
print('已生成: %s' % OUT)
print('工作表: %s' % wb.sheetnames)
