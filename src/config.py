#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""统一路径与参数解析（提交包内部模块）

本项目所有脚本都通过本模块取路径，路径与参数集中在根目录 config.yaml。

    from config import ROOT, HERE, path, CFG
    rec = path(*CFG['paths']['receptor'].split('/'))
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))       # project/src
ROOT = os.path.dirname(HERE)                            # project/

if HERE not in sys.path:
    sys.path.insert(0, HERE)

CFG_PATH = os.path.join(ROOT, 'config.yaml')
CFG = {}
try:
    import yaml
    with open(CFG_PATH, encoding='utf-8') as _f:
        CFG = yaml.safe_load(_f) or {}
except ImportError:
    # 没有 PyYAML 时退化为极简解析（只支持本项目的两层结构）
    _sec = None
    for _line in open(CFG_PATH, encoding='utf-8'):
        _s = _line.rstrip('\n')
        if not _s.strip() or _s.lstrip().startswith('#'):
            continue
        if not _s.startswith((' ', '\t')) and _s.rstrip().endswith(':'):
            _sec = _s.strip()[:-1]
            CFG[_sec] = {}
            continue
        if _sec and ':' in _s:
            _k, _v = _s.strip().split(':', 1)
            _v = _v.split('#')[0].strip()
            if _v == '':
                CFG[_sec][_k] = {}
            else:
                try:
                    CFG[_sec][_k] = eval(_v, {'__builtins__': {}}, {})
                except Exception:
                    CFG[_sec][_k] = _v.strip('"\'')


def path(*parts, **kw):
    """拼出 project/ 下的相对路径；kw: mkdir=True 时自动建父目录"""
    p = os.path.join(ROOT, *[str(x) for x in parts])
    if kw.get('mkdir'):
        os.makedirs(os.path.dirname(p) if os.path.splitext(p)[1] else p,
                    exist_ok=True)
    return p


def cfg(*keys, default=None):
    """按层级取配置：cfg('docking', 'seed')"""
    cur = CFG
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def ensure_work():
    """创建运行中间目录"""
    w = path(cfg('paths', 'work', default='work'))
    os.makedirs(w, exist_ok=True)
    return w


if __name__ == '__main__':
    print('ROOT =', ROOT)
    print('HERE =', HERE)
    print('config:', CFG_PATH, '存在' if os.path.exists(CFG_PATH) else '缺失')
    for k, v in CFG.items():
        print('  [%s] %s' % (k, v if not isinstance(v, dict) else list(v.keys())))
    print('receptor ->', path(cfg('paths', 'receptor')))
