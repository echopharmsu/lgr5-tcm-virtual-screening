#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成最终结果文件（附件5 第四章）

必填字段：候选编号 / 所属赛道 / SMILES 结构 / 结构文件 / 关键预测指标 /
          对应模型与运行版本 / 备注
本项目在其上补充：中文名、来源中药、各区域指标与判定、推荐区域。

输出（附件5 第四章：以模板为准；未另行规定时默认 results.xlsx 或 UTF-8 results.csv）：
    results/results.xlsx    —— 若环境装有 openpyxl
    results/results.csv     —— 总是写出（UTF-8 with BOM，Excel 可直接打开）

用法: python3 make_results.py [结果目录]
"""
import csv
import os
import sys

try:
    import openpyxl
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter
    HAVE_XLSX = True
except ImportError:
    HAVE_XLSX = False

from config import ROOT, path, cfg

RES = sys.argv[1] if len(sys.argv) > 1 else path(cfg('paths', 'results'))
SITES = list((cfg('sites', default={}) or {}).keys())
REC_SITE = cfg('project', 'recommended_site', default='U3')          # 内部键：查 cand / 结构文件名
REC_LABEL = cfg('project', 'recommended_site_label', default=REC_SITE)   # 对外显示名：交付表「推荐位点」列
TARGET = cfg('project', 'recommended_molecule', default='001')
TRACK = cfg('project', 'track', default='')
MODEL_VER = '%s（CNN 模型 %s，--cnn_scoring rescore）' % (
    cfg('docking', 'software', default='gnina'), cfg('docking', 'cnn_model', default='dense'))

COLS = ['候选编号', '中文名', '英文名', '来源中药', '所属赛道', 'SMILES',
        '结构文件', '推荐位点', '推荐位点判定', 'CNNscore', 'CNN预测亲和力',
        '结合能(kcal/mol)', '距位点中心(A)', '接触残基数', '氢键残基数',
        '对应模型与运行版本', '随机种子']

# ADME 列（若 results/ 下有 ADME 明细则追加；来源为第三方工具 ADMET-AI，不参与任何判定）
ADME_COLS = ['血脑屏障安全度', 'hERG安全度', '无毒性评分', '溶解度评分', '生物利用度评分']
ADME_HEAD = {'血脑屏障安全度': 'Blood-Brain Barrier Safe', 'hERG安全度': 'hERG Safe',
             '无毒性评分': 'Non-Toxic', '溶解度评分': 'Soluble',
             '生物利用度评分': 'Bioavailable'}


def read_tsv(p):
    if not os.path.exists(p):
        return []
    with open(p, encoding='utf-8-sig', errors='replace') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def fn(s, d=None):
    try:
        return float(s)
    except (TypeError, ValueError):
        return d


def intn(s, d=None):
    """计数型字段：写成**整数**而非浮点。

    2026-09-29 修复：接触残基数 / 氢键残基数是计数，用 fn() 会变成 float，
    csv.writer 于是写出 "6.0" 而 xlsx 显示 "6" —— 两份"同内容"文件在逐字对表时
    会出现 220 个单元格的不一致。改用本函数写整数。
    读不出数值时返回 d（与 fn 的缺省语义一致，csv 会写成空）。
    """
    v = fn(s)
    if v is None:
        return d
    return int(round(v))


# ---------------------------------------------------------------- 输入
meta = {r['id']: r for r in read_tsv(path(cfg('paths', 'ligand_source')))}
cand = {s: {r['编号']: r for r in
            read_tsv(os.path.join(RES, 'screen_候选清单_%s.tsv' % s))} for s in SITES}
if not meta:
    sys.exit('ERROR: 读不到配体表 %s' % path(cfg('paths', 'ligand_source')))
if not any(cand.values()):
    sys.exit('ERROR: 结果目录 %s 下没有候选清单，先跑 predict.py' % RES)

# 排序（2026-09-29 修正）：
#   ① **不再把 TARGET(001) 强行置顶** —— 固定在第 1 位属自我印证；现与其它分子同一把尺子。
#   ② **推荐位点 PASS 的排在 FAIL 之前**——候选清单不该以"未通过"开头。
#   ③ 组内按推荐位点的 CNNscore 降序。
order = sorted([i for i in meta if any(i in cand.get(s, {}) for s in SITES)],
               key=lambda i: (
    0 if (cand.get(REC_SITE, {}).get(i) or {}).get('判定') == 'PASS' else 1,
    -(fn((cand.get(REC_SITE, {}).get(i) or {}).get('CNN打分(0-1越高越好)'), 0.0) or 0.0)))
if not order:
    sys.exit('ERROR: 候选清单里没有任何配体，先跑 predict.py')

# ---------------------------------------------------------------- ADME（可选）
ADME = {}
ADME_SRC = ''
adme_path = path(cfg('paths', 'adme', default='')) if cfg('paths', 'adme', default='') else ''
if adme_path and os.path.exists(adme_path):
    if not HAVE_XLSX:
        print('提示: 无 openpyxl，读不了 ADME 明细 %s → 不追加 ADME 列' % adme_path)
    else:
        _wb = openpyxl.load_workbook(adme_path, data_only=True, read_only=True)
        _ws = _wb['总览']
        _col, _hdr = {}, None
        for _rr in _ws.iter_rows(values_only=True):
            if _hdr is None:
                _hdr = ['' if x is None else str(x) for x in _rr]
                _col = {n: _hdr.index(h) for n, h in ADME_HEAD.items() if h in _hdr}
                continue
            if not _rr or _rr[0] is None:
                continue
            ADME[str(_rr[0]).strip()] = {n: _rr[_col[n]] for n in _col}
        _wb.close()
        ADME_SRC = os.path.relpath(adme_path, ROOT).replace(os.sep, '/')
        print('ADME 明细已载入: %s（%d 个分子，%d 列）' % (ADME_SRC, len(ADME), len(_col)))
if not ADME:
    print('★ 警告: 未载入 ADME 明细 → 结果表不追加 ADME 列（将少 5 列，23 → 18）。')
    print('        本次输出会改用带「无ADME_非交付件」标记的文件名，【不会覆盖交付件】。')

COLS = COLS + (ADME_COLS if ADME else []) + ['备注']

rows = []
for i in order:
    m = meta.get(i, {})
    c = cand.get(REC_SITE, {}).get(i, {})
    struct = ('results/结构文件/%s_%s_%s复合物.pdb'
              % (i, m.get('name_en', '').replace(' ', '_'), REC_SITE)
              if i == TARGET else '')
    note = []
    if i == TARGET:
        note.append('★ 推荐候选（库内选择性检查：同尺寸子集第 1 名）')
    if c.get('判定') == 'FAIL':
        note.append('推荐位点未通过：%s' % (c.get('未通过原因') or ''))
    if not c:
        note.append('该位点无结果')
    row = [i, m.get('name_cn', ''), m.get('name_en', ''), m.get('herb', ''), TRACK,
           m.get('smiles', ''), struct, REC_LABEL, c.get('判定', '-'),
           fn(c.get('CNN打分(0-1越高越好)')), fn(c.get('CNN预测亲和力')),
           fn(c.get('结合能(kcal/mol越负越好)')), fn(c.get('距位点中心(A，需≤6)')),
           intn(c.get('接触残基数')), intn(c.get('氢键残基数')),
           MODEL_VER, cfg('docking', 'seed', default=42)]
    if ADME:
        row += [ADME.get(i, {}).get(n, '') for n in ADME_COLS]
    row.append('；'.join(x for x in note if x))
    rows.append(row)

NOTES = [
    ('项目', cfg('project', 'name', default='')),
    ('赛道', TRACK),
    ('靶点', '%s；结构来源 RCSB PDB %s'
     % (cfg('project', 'target', default=''), cfg('project', 'pdb_id', default=''))),
    ('推荐分子 / 区域', '%s（%s） / %s'
     % (TARGET, meta.get(TARGET, {}).get('name_cn', ''), REC_LABEL)),
    ('', ''),
    ('对接软件', '%s；CNN 模型 %s' % (cfg('docking', 'software', default=''),
                                 cfg('docking', 'cnn_model', default=''))),
    ('对接参数', '盒 %s Å 立方 · exhaustiveness %s · num_modes %s · --seed %s'
     % (cfg('docking', 'box_size'), cfg('docking', 'exhaustiveness'),
        cfg('docking', 'num_modes'), cfg('docking', 'seed'))),
    ('受体', 'LGR5 胞外域 A 链（残基 30-543，3616 重原子）；来源 4KNG，原始下载件见 data/4KNG_full.pdb'),
    ('配体库', '%d 个中药小分子；立体构型经立体化学闸门校验' % len(meta)),
    ('过滤判据', '最近接触 ≥ %s Å 且 affinity < %s 且 CNNscore ≥ %s 且 距区域中心 ≤ %s Å'
     % (cfg('filters', 'min_contact'), cfg('filters', 'max_affinity'),
        cfg('filters', 'min_cnn'), cfg('filters', 'site_cut'))),
    ('结果生成入口', 'python3 predict.py'),
    ('', ''),
    ('重要限制', '① 受体刚性、仅胞外域（不含 7TM 与胞内区），484-536 共 53 残基未建模；'
              '② 本靶点无任何实验解析的小分子配体，无法做重对接黄金对照；'
              '③ 目标区域埋藏比约 0.30，属浅表面结合，非真实口袋；'
              '④ 结果为计算假说，不代表已验证的结合区域；'
              '⑤ CNNscore 不可跨轮/跨受体网格比较'),
]
if ADME:
    NOTES.insert([a for a, _ in NOTES].index('结果生成入口'),
                 ('ADME 列来源', '%s —— 第三方工具 ADMET-AI（Chemprop 系）'
                  '在本机单独运行产出，取其中 52 个端点的 5 个安全度轴。'
                  '该预测不在对接流程内，不参与任何排序、筛选或判定，也不可与 gnina 分数'
                  '比较或加权（量纲与含义不同）；数值为计算结果而非实验数据。' % ADME_SRC))

# ---------------------------------------------------------------- 输出 CSV
# ★ 只有 ADME 明细成功载入时才写「交付件」文件名；否则改用带标记的名字，
#   避免把 18 列的不完整表当成 23 列交付件覆盖掉（2026-10-11 事故的直接教训）。
csv_out = os.path.join(RES, 'results.csv' if ADME else 'results_无ADME_非交付件.csv')
with open(csv_out, 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f)
    for a, b in NOTES:
        if a:
            w.writerow(['# %s' % a, b])
    w.writerow([])
    w.writerow(COLS)
    w.writerows(rows)
print('已写出（UTF-8 CSV）: %s' % csv_out)

# ---------------------------------------------------------------- 输出 XLSX
if HAVE_XLSX:
    HDR = PatternFill('solid', fgColor='2F5597')
    HDRF = Font(bold=True, color='FFFFFF', size=10)
    OKF = PatternFill('solid', fgColor='E2EFDA')
    BADF = PatternFill('solid', fgColor='FCE4E4')
    GOLD = PatternFill('solid', fgColor='FFE699')
    THIN = Side(style='thin', color='BFBFBF')
    BD = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = '候选清单'
    ws.append(COLS)
    for r in rows:
        ws.append(r)
    for c in range(1, len(COLS) + 1):
        cell = ws.cell(row=1, column=c)
        cell.fill, cell.font, cell.border = HDR, HDRF, BD
        cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
    for r in range(2, ws.max_row + 1):
        is_tgt = ws.cell(row=r, column=1).value == TARGET
        is_ok = ws.cell(row=r, column=9).value == 'PASS'
        for c in range(1, len(COLS) + 1):
            cell = ws.cell(row=r, column=c)
            cell.border = BD
            cell.font = Font(size=9)
            cell.alignment = Alignment(horizontal='left', vertical='center')
            if is_tgt:
                cell.fill = GOLD
            elif not is_ok:
                cell.fill = BADF
            elif c == 9:
                cell.fill = OKF
    WIDTHS = [10, 14, 22, 10, 22, 46, 40, 9, 12, 10, 13, 13, 13, 10, 10,
              12, 12, 12, 26, 9]
    WIDTHS += [11] * (len(ADME_COLS) if ADME else 0)
    WIDTHS += [26]                      # 备注
    for i, wd in enumerate(WIDTHS, 1):
        ws.column_dimensions[get_column_letter(i)].width = wd
    ws.freeze_panes = 'A2'
    ws.row_dimensions[1].height = 30
    ws.auto_filter.ref = 'A1:%s%d' % (get_column_letter(len(COLS)), ws.max_row)

    ws2 = wb.create_sheet('说明')
    ws2.append(['项目', '说明'])
    for a, b in NOTES:
        ws2.append([a, b])
    for c in (1, 2):
        ws2.cell(row=1, column=c).fill = HDR
        ws2.cell(row=1, column=c).font = HDRF
    ws2.column_dimensions['A'].width = 20
    ws2.column_dimensions['B'].width = 120
    for r in range(2, ws2.max_row + 1):
        ws2.cell(row=r, column=1).font = Font(size=10, bold=True)
        ws2.cell(row=r, column=2).font = Font(size=10)
        ws2.cell(row=r, column=2).alignment = Alignment(wrap_text=True, vertical='top')

    xlsx_out = os.path.join(RES, 'results.xlsx' if ADME
                            else 'results_无ADME_非交付件.xlsx')
    wb.save(xlsx_out)
    print('已写出（Excel）: %s   工作表 %s' % (xlsx_out, wb.sheetnames))
else:
    print('★ 警告: 本环境没有 openpyxl —— 既未生成 results.xlsx，也读不了 ADME 明细，')
    print('        因此【没有写出 23 列的交付件】。要生成完整交付件：')
    print('            pip install openpyxl 后重跑本脚本。')

print('候选行数 %d ；字段数 %d' % (len(rows), len(COLS)))
npass = sum(1 for r in rows if r[8] == 'PASS')
print('推荐区域 %s：通过 %d / 未通过 %d' % (REC_SITE, npass, len(rows) - npass))
