#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""虚拟筛选入口：把配体做 3D 化并投到指定区域做 gnina 对接（附件5 第二章 screen.py）

用法：
    python3 screen.py --site U3                       # 全库 × U3
    python3 screen.py --site U3 --limit 5 --exh 4     # 演示：5 个分子
    python3 screen.py --site U3 --skip-prep           # 复用已 3D 化的配体

流程：
    1) 读 config.yaml 取受体 / 配体库 / 盒参数 / 随机种子
    2) 调 src/dock_prep.py 把 2D SDF 库做 ETKDGv3(seed) 3D 化 + 平移到盒中心
       （内含立体化学闸门：未指定手性的分子会被修复/枚举/拦截，见 stereo_gate.py）
    3) 调 gnina 对接，原始输出与日志分别落到 results/ 与 logs/
"""
import argparse
import os
import subprocess
import sys
import time

from config import ROOT, HERE, path, cfg

SEP = '$$$$\n'


def sh(cmd, logfile):
    """跑外部命令并把 stdout+stderr 写进日志；返回退出码（附件5 第三条：调用日志须留存）"""
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    with open(logfile, 'w', encoding='utf-8') as lf:
        lf.write('# %s\n# cwd: %s\n\n' % (' '.join(cmd), os.getcwd()))
        lf.flush()
        p = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    return p.returncode


def count_records(p):
    import gzip
    if not os.path.exists(p):
        return 0
    with gzip.open(p, 'rb') as f:
        return f.read().decode('utf-8', 'replace').count(SEP)


def main():
    ap = argparse.ArgumentParser(description='gnina 虚拟筛选（单一区域）')
    ap.add_argument('--site', required=True, help='区域键名（config.yaml 的 sites 键，如 U3）')
    ap.add_argument('--limit', type=int, default=0, help='只跑前 N 个配体（0=全部）')
    ap.add_argument('--exh', type=int, default=None, help='覆盖 exhaustiveness')
    ap.add_argument('--skip-prep', action='store_true', help='复用已 3D 化的配体')
    ap.add_argument('--out', default=None, help='输出 sdf.gz 路径（默认 results/screen_<site>.sdf.gz）')
    ap.add_argument('--log', default=None,
                    help='gnina 日志路径（默认 logs/screen_<site>.log）——'
                         '演示运行请传 logs/demo_screen_<site>.log，避免覆盖全量轮日志')
    args = ap.parse_args()

    sites = cfg('sites', default={}) or {}
    if args.site not in sites:
        sys.exit('ERROR: config.yaml 里没有区域 %s（可选：%s）'
                 % (args.site, ', '.join(sites)))
    cx, cy, cz = sites[args.site]['center']

    box = cfg('docking', 'box_size', default=24)
    exh = args.exh if args.exh else cfg('docking', 'exhaustiveness', default=8)
    modes = cfg('docking', 'num_modes', default=9)
    seed = cfg('docking', 'seed', default=42)
    cpu = cfg('docking', 'cpu', default=0)

    rec = path(cfg('paths', 'receptor'))
    lig2d = path(cfg('paths', 'ligands_2d'))
    if not os.path.exists(rec):
        sys.exit('ERROR: 受体不存在 %s' % rec)
    if not os.path.isdir(lig2d):
        sys.exit('ERROR: 配体库不存在 %s\n先运行: python3 predict.py --build-library' % lig2d)

    work = path(cfg('paths', 'work'), 'screen_%s' % args.site)
    os.makedirs(work, exist_ok=True)
    lig3d = os.path.join(work, 'ligands_3d.sdf')
    logdir = path(cfg('paths', 'logs'))
    out = args.out or path(cfg('paths', 'results'), 'screen_%s.sdf.gz' % args.site)

    print('=== 区域 %s  中心 (%.2f, %.2f, %.2f)  盒 %.0f Å ===' % (args.site, cx, cy, cz, box))
    print('受体 %s' % rec)

    # ---- 1) 配体 3D 化（含立体化学闸门）----
    if not args.skip_prep or not os.path.exists(lig3d):
        cmd = [sys.executable, os.path.join(HERE, 'dock_prep.py'),
               lig2d, os.path.join(work, 'ligands_3d'),
               lig3d, os.path.join(work, 'prep_report.tsv'),
               str(cx), str(cy), str(cz),
               '--stereo-policy', cfg('ligand_prep', 'stereo_policy', default='auto'),
               '--stereo-max-isomers', str(cfg('ligand_prep', 'stereo_max_isomers', default=8))]
        print('--- 配体 3D 化 ---')
        rc = sh(cmd, os.path.join(logdir, 'prep_%s.log' % args.site))
        if rc != 0:
            sys.exit('ERROR: dock_prep 失败，见 logs/prep_%s.log' % args.site)
    else:
        print('--- 复用已有 3D 配体 %s ---' % lig3d)

    # 演示/限量模式：截取前 N 个分子
    if args.limit and args.limit > 0:
        import gzip
        with open(lig3d, 'r', errors='replace') as f:
            txt = f.read()
        recs = [r for r in txt.split(SEP) if r.strip()]
        recs = recs[:args.limit]
        lig3d_use = os.path.join(work, 'ligands_3d_top%d.sdf' % args.limit)
        with open(lig3d_use, 'w') as f:
            f.write(SEP.join(recs) + SEP)
        print('--- 限量模式：只跑前 %d 个分子 ---' % args.limit)
    else:
        lig3d_use = lig3d

    # ---- 2) gnina 对接 ----
    exe = 'gnina'
    if subprocess.run(['which', exe], capture_output=True).returncode != 0:
        sys.exit('ERROR: 找不到 gnina。安装说明见 README「环境配置」；'
                 '或把可执行文件放进 PATH。')
    cmd = [exe, '-r', rec, '-l', lig3d_use,
           '--center_x', str(cx), '--center_y', str(cy), '--center_z', str(cz),
           '--size_x', str(box), '--size_y', str(box), '--size_z', str(box),
           '--exhaustiveness', str(exh), '--num_modes', str(modes),
           '--seed', str(seed), '-o', out]
    if cpu and int(cpu) > 0:
        cmd += ['--cpu', str(cpu)]
    print('--- gnina 对接（exh %d / modes %d / seed %d）---' % (exh, modes, seed))
    t0 = time.time()
    logfile = args.log or os.path.join(logdir, 'screen_%s.log' % args.site)
    rc = sh(cmd, logfile)
    dt = time.time() - t0
    if rc != 0:
        sys.exit('ERROR: gnina 退出码 %d，见 %s' % (rc, logfile))
    print('完成：%d 条姿态，用时 %.1f s -> %s' % (count_records(out), dt, out))
    print('日志：%s' % logfile)


if __name__ == '__main__':
    main()
