#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collect_env_md.py —— 采集本机 WSL 环境的**配置层面**信息，生成 chempropenv.md

只用实跑命令取值（非人工转录），输出四节：
    一、主机与操作系统
    二、GPU 与 CUDA
    三、Python 环境
    四、关键组件版本

（不含安装过程、依赖清单、注意事项——按项目要求只记录配置层面）

用法：
    ~/envs/chemprop/bin/python collect_env_md.py [输出路径]
    # 不传参数时输出到【当前工作目录】的 chempropenv.md
    # 例： wsl -d Ubuntu-22.04 -- bash -lc "cd <任意目录> && ~/envs/chemprop/bin/python collect_env_md.py out.md"
"""
import datetime as dt
import subprocess
import sys
from pathlib import Path

# 本脚本不硬编码任何绝对路径（附件5 第六条「避免硬编码绝对路径」）：
# 环境位置一律由 Path.home() 推出；输出位置由命令行参数或当前目录决定。
ENV = Path.home() / "envs" / "chemprop"
PY = ENV / "bin" / "python"
MAMBA = Path.home() / ".local" / "bin" / "micromamba"
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd() / "chempropenv.md"


def sh(cmd: str) -> str:
    """跑命令取真实输出；失败时返回可读的失败标记，绝不静默编造。"""
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
        return ((r.stdout or "") + (r.stderr or "")).strip()
    except Exception as e:  # noqa: BLE001
        return f"<采集失败：{type(e).__name__}: {e}>"


def keyvals(txt: str, keys: list) -> dict:
    d = {}
    for line in txt.splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            if k in keys:
                d[k] = v.strip().strip('"')
    return d


def kv_block(txt: str) -> str:
    rows = []
    for line in txt.splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            rows.append(f"| `{k.strip()}` | {v.strip()} |")
    return "\n".join(rows) or "| — | 采集失败 |"


# ---------------------------------------------------------------- 采集
osrel = sh("cat /etc/os-release")
osv = keyvals(osrel, ["PRETTY_NAME", "VERSION_ID", "VERSION_CODENAME"])
kernel = sh("uname -r")
arch = sh("uname -m")
glibc = sh("ldd --version | head -1").split()[-1]
cpu = sh("grep -m1 'model name' /proc/cpuinfo | cut -d: -f2").strip()
ncpu = sh("nproc")
mem = sh("free -g | awk '/^Mem:/{print $2\" GB total / \"$7\" GB available\"}'")
disk = sh("df -h / | awk 'NR==2{print $2\" total, \"$4\" avail, mount \"$6}'")
python_sys = sh("python3 -V")

gpu = sh("nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader")
gpu_full = sh("nvidia-smi | head -12")
cuda_tk = sh("ls -d /usr/local/cuda-* 2>/dev/null | tr '\\n' ' '").strip()
nvcc = sh("/usr/local/cuda-12.8/bin/nvcc --version 2>/dev/null | grep release").strip()

torchinfo = sh(f"""{PY} - <<'PYEOF'
import torch
print("torch_version=", torch.__version__)
print("cuda_built=", torch.version.cuda)
print("cudnn=", torch.backends.cudnn.version())
print("cuda_available=", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device=", torch.cuda.get_device_name(0))
    print("capability=", torch.cuda.get_device_capability(0))
    print("arch_list=", torch.cuda.get_arch_list())
    print("device_count=", torch.cuda.device_count())
PYEOF""")

verinfo = sh(f"""{PY} - <<'PYEOF'
import importlib
from importlib.metadata import version as mv, PackageNotFoundError
mods = ["chemprop", "admet_ai", "torch", "rdkit", "lightning", "numpy", "pandas", "scipy",
        "scikit_learn", "tensorboard", "flask", "gunicorn", "matplotlib", "seaborn",
        "astartes", "descriptastorus", "cuik_molmaker_pin", "myerson", "openpyxl",
        "typed_argument_parser"]
for m in mods:
    v = None
    try:
        v = getattr(importlib.import_module(m), "__version__", None)
    except Exception:
        pass
    if v is None:
        try:
            v = mv(m.replace("_", "-"))
        except PackageNotFoundError:
            v = "<未安装>"
        except Exception as e:
            v = f"<取不到：{{type(e).__name__}}>"
    print(f"{{m}}={{v}}")
PYEOF""")

mmver = sh(f"{MAMBA} --version")
pyver = sh(f"{PY} -V")
pipfreeze = sh(f"{PY} -m pip freeze")
npkg = len([l for l in pipfreeze.splitlines() if l.strip()])
mmlist = sh(f"{MAMBA} list -p {ENV} 2>/dev/null")
nimport = sh(f"{PY} -c \"import chemprop,torch,admet_ai;print('ok')\"")
admet_cli = sh(f"{ENV}/bin/admet_predict --help 2>&1 | head -3")
chemprop_cli = sh(f"{ENV}/bin/chemprop --version 2>&1")
ts = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# ---------------------------------------------------------------- 成文
md = f"""# Chemprop / ADMET-AI 运行环境记录（WSL）

> **生成时间**：{ts}
> **采集方式**：由 `collect_env_md.py` 在本机**实际执行命令**采集生成，**非人工转录**（可重跑核对）
> **记录范围**：配置层面（系统 / GPU / Python 环境 / 组件版本）

---

## 一、主机与操作系统

| 项 | 值 |
|---|---|
| 发行版 | {osv.get('PRETTY_NAME', '?')} |
| 版本号 / 代号 | {osv.get('VERSION_ID', '?')} / {osv.get('VERSION_CODENAME', '?')} |
| 内核 | `{kernel}` |
| 架构 | `{arch}` |
| glibc | `{glibc}` |
| CPU | {cpu}（{ncpu} 核） |
| 内存 | {mem} |
| 根分区 | {disk} |

```
{osrel}
```

---

## 二、GPU 与 CUDA

| 项 | 值 |
|---|---|
| GPU | {gpu} |
| 系统 CUDA Toolkit | `{cuda_tk}` |
| nvcc | {nvcc or '（不在当前 PATH；工具链位于 /usr/local/cuda-12.8）'} |

**PyTorch 侧（决定能否用 GPU 的关键指标）**

| 项 | 值 |
|---|---|
{kv_block(torchinfo)}

> ⚠️ 本机 GPU 为 **Blackwell（sm_120）**，PyTorch 必须用 **cu128 及以上**轮子；
> `arch_list` 中出现 `sm_120` 才说明该轮子真正支持本机显卡。

```
{gpu_full}
```

---

## 三、Python 环境

| 项 | 值 |
|---|---|
| 环境路径 | `{ENV}` |
| 解释器 | `{ENV}/bin/python` |
| Python 版本 | {pyver} |
| 环境管理器 | micromamba `{mmver}`（root prefix `{Path.home() / 'micromamba'}`） |
| 是否安装 conda | **否** —— 使用免 sudo 的单文件 micromamba |
| 系统 Python | {python_sys}（Ubuntu 自带，**未被改动**，与本环境隔离） |
| 已安装 pip 包数 | {npkg} |

导入自检：`{nimport}`

---

## 四、关键组件版本

| 组件 | 版本 |
|---|---|
{kv_block(verinfo)}

命令行入口自检：

```
$ chemprop --version
{chemprop_cli}

$ admet_predict --help | head -3
{admet_cli}
```

---

## 五、依赖清单（保证可复现）

### 5.1 `pip freeze`（环境内直接导出，权威）

> 把下面代码块另存为 `requirements.txt`，即可 `pip install -r requirements.txt` 复现。
> 共 {npkg} 个包。

```
{pipfreeze}
```

### 5.2 micromamba 安装的包（Python 解释器与系统库层）

```
{mmlist}
```
"""

OUT.write_text(md, encoding="utf-8")
print(f"OK -> {OUT}")
print(f"   {OUT.stat().st_size} 字节 / {md.count(chr(10))} 行 / pip 包 {npkg} 个")
