#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""★ 主运行入口（附件5 第二章 predict.py）——一键生成最终结果文件

用法：
    python3 predict.py                    # 全流程（109 分子 × 全部区域 → results/results.xlsx）
    python3 predict.py --demo             # 演示案例（小规模，约 1 分钟；输出到 results/demo/）
    python3 predict.py --build-library    # 只建配体库
    python3 predict.py --report-only      # 复用已有对接结果，只重算清单与 Excel
    python3 predict.py --sites <区域键>    # 只跑指定区域（键名见 config.yaml 的 sites）

流程（每一步都会打印进度并落日志到 logs/）：
    ① 建配体 2D 库（由 data/ligand_source.tsv 的 SMILES 生成）
    ② 逐区域：配体 3D 化（ETKDGv3, seed 42，含立体化学闸门）→ gnina 对接
    ③ 区域分析：逐分子指标、尺寸偏倚
    ④ 四道过滤 → 候选清单
    ⑤ 输出 results/results.xlsx（附件5 第四章要求的字段）

环境：见 README「环境配置」；所有参数来自 config.yaml。
"""
import argparse
import os
import subprocess
import sys
import time

# 本文件位于 project/ 根目录，公共模块在 src/ 下 —— 先注入路径再导入
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))

from config import ROOT, HERE, path, cfg   # noqa: E402


def run(step, cmd, logname=None):
    print('\n' + '=' * 72)
    print('▶ %s' % step)
    print('  %s' % ' '.join(str(c) for c in cmd))
    print('=' * 72)
    t0 = time.time()
    p = subprocess.run([str(c) for c in cmd])
    dt = time.time() - t0
    if p.returncode != 0:
        sys.exit('✗ %s 失败（退出码 %d）%s' % (step, p.returncode,
                                           '  日志: %s' % logname if logname else ''))
    print('✓ %s 完成（%.1f s）' % (step, dt))


def main():
    ap = argparse.ArgumentParser(description='LGR5 中药小分子虚拟筛选 · 主运行入口')
    ap.add_argument('--demo', action='store_true', help='小规模演示案例')
    ap.add_argument('--sites', default=None, help='逗号分隔的区域（默认全部）')
    ap.add_argument('--build-library', action='store_true', help='只建配体库')
    ap.add_argument('--report-only', action='store_true', help='跳过对接，只重算清单与 Excel')
    args = ap.parse_args()

    all_sites = list((cfg('sites', default={}) or {}).keys())
    if args.demo:
        sites = cfg('demo', 'sites', default=['U3'])
        limit = cfg('demo', 'n_ligands', default=5)
        exh = cfg('demo', 'exhaustiveness', default=4)
        outdir = path(cfg('paths', 'results'), 'demo')
        print('★ 演示模式：%d 个配体 × %s（exhaustiveness %d）→ %s'
              % (limit, ','.join(sites), exh, outdir))
    else:
        sites = [s.strip() for s in args.sites.split(',')] if args.sites else all_sites
        limit = 0
        exh = None
        outdir = path(cfg('paths', 'results'))
    os.makedirs(outdir, exist_ok=True)

    src = HERE
    py = sys.executable

    # ---- ① 配体 2D 库 ----
    ligdir = path(cfg('paths', 'ligands_2d'))
    n_sdf = len([f for f in os.listdir(ligdir) if f.endswith('.sdf')]) \
        if os.path.isdir(ligdir) else 0
    if n_sdf == 0 or args.build_library:
        run('① 建配体 2D 库',
            [py, os.path.join(src, 'build_library.py'),
             path(cfg('paths', 'ligand_source')), ligdir],
            'logs/build_library.log')
    else:
        print('① 配体库已存在：%s（%d 个分子）' % (ligdir, n_sdf))
    if args.build_library:
        return

    # ---- ② 逐区域对接 ----
    if not args.report_only:
        for s in sites:
            # 演示运行与全量运行必须写【不同】的日志文件，否则 demo 会覆盖全量轮日志
            logname = 'logs/%sscreen_%s.log' % ('demo_' if args.demo else '', s)
            cmd = [py, os.path.join(src, 'screen.py'), '--site', s,
                   '--out', os.path.join(outdir, 'screen_%s.sdf.gz' % s),
                   '--log', logname]
            if limit:
                cmd += ['--limit', str(limit)]
            if exh:
                cmd += ['--exh', str(exh)]
            run('② 区域 %s 对接' % s, cmd, logname)

    # ---- ③ 区域分析 ----
    sdfs = [os.path.join(outdir, 'screen_%s.sdf.gz' % s) for s in sites]
    missing = [f for f in sdfs if not os.path.exists(f)]
    if missing:
        sys.exit('✗ 缺少对接结果：%s\n  先去掉 --report-only 重跑' % ', '.join(missing))
    run('③ 区域分析（逐分子指标 / 尺寸偏倚）',
        [py, os.path.join(src, 'screen_analyze.py'),
         path(cfg('paths', 'ligand_source')), path(cfg('paths', 'receptor')),
         os.path.join(outdir, 'screen_unified')] + sdfs,
        'logs/screen_analyze.log')

    # ---- ④ 四道过滤 → 候选清单 ----
    run('④ 四道过滤 → 候选清单',
        [py, os.path.join(src, 'make_deliverables.py'),
         path(cfg('paths', 'ligand_source')), path(cfg('paths', 'receptor')),
         os.path.join(outdir, 'screen')] + sdfs,
        'logs/make_deliverables.log')

    # ---- ⑤ 最终结果文件 ----
    run('⑤ 生成 results.xlsx（附件5 第四章字段）',
        [py, os.path.join(src, 'make_results.py'), outdir],
        'logs/make_results.log')

    print('\n' + '=' * 72)
    print('全部完成。最终结果文件：')
    print('  %s' % os.path.join(outdir, 'results.xlsx'))
    print('候选清单：')
    for s in sites:
        print('  %s' % os.path.join(outdir, 'screen_候选清单_%s.tsv' % s))
    print('日志目录：%s' % path(cfg('paths', 'logs')))


if __name__ == '__main__':
    main()
