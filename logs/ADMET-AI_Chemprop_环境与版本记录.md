<!-- 本文件由项目根目录的 chempropenv.md 原件复制而来（2026-09-30）。 -->
<!-- 采集脚本随包：logs/collect_env_md.py（已去除硬编码绝对路径，可直接重跑复核）。 -->

# Chemprop / ADMET-AI 运行环境记录（WSL）

> **生成时间**：2026-09-28 22:18:21
> **采集方式**：由 **`logs/collect_env_md.py`** 在本机**实际执行命令**采集生成，**非人工转录**
> —— **该脚本随包提交**，可从任意目录重跑核对（不硬编码绝对路径，输出路径由命令行参数或当前目录决定）：
>
> ```bash
> ~/envs/chemprop/bin/python logs/collect_env_md.py <输出路径>
> ```
> **记录范围**：配置层面（系统 / GPU / Python 环境 / 组件版本）

---

## 一、主机与操作系统

| 项 | 值 |
|---|---|
| 发行版 | Ubuntu 22.04.5 LTS |
| 版本号 / 代号 | 22.04 / jammy |
| 内核 | `6.18.33.2-microsoft-standard-WSL2` |
| 架构 | `x86_64` |
| glibc | `2.35` |
| CPU | Intel(R) Core(TM) Ultra 5 225H（14 核） |
| 内存 | 15 GB total / 14 GB available |
| 根分区 | 1007G total, 919G avail, mount / |

```
PRETTY_NAME="Ubuntu 22.04.5 LTS"
NAME="Ubuntu"
VERSION_ID="22.04"
VERSION="22.04.5 LTS (Jammy Jellyfish)"
VERSION_CODENAME=jammy
ID=ubuntu
ID_LIKE=debian
HOME_URL="https://www.ubuntu.com/"
SUPPORT_URL="https://help.ubuntu.com/"
BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
UBUNTU_CODENAME=jammy
```

---

## 二、GPU 与 CUDA

| 项 | 值 |
|---|---|
| GPU | NVIDIA GeForce RTX 5060 Laptop GPU, 591.74, 8151 MiB |
| 系统 CUDA Toolkit | `/usr/local/cuda-12 /usr/local/cuda-12.8` |
| nvcc | Cuda compilation tools, release 12.8, V12.8.93 |

**PyTorch 侧（决定能否用 GPU 的关键指标）**

| 项 | 值 |
|---|---|
| `torch_version` | 2.11.0+cu128 |
| `cuda_built` | 12.8 |
| `cudnn` | 91900 |
| `cuda_available` | True |
| `device` | NVIDIA GeForce RTX 5060 Laptop GPU |
| `capability` | (12, 0) |
| `arch_list` | ['sm_75', 'sm_80', 'sm_86', 'sm_90', 'sm_100', 'sm_120'] |
| `device_count` | 1 |

> ⚠️ 本机 GPU 为 **Blackwell（sm_120）**，PyTorch 必须用 **cu128 及以上**轮子；
> `arch_list` 中出现 `sm_120` 才说明该轮子真正支持本机显卡。

```
Mon Sep 28 22:17:58 2026       
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 590.52.01              Driver Version: 591.74         CUDA Version: 13.1     |
+-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA GeForce RTX 5060 ...    On  |   00000000:01:00.0  On |                  N/A |
| N/A   44C    P5              7W /  115W |    1352MiB /   8151MiB |      1%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
```

---

## 三、Python 环境

| 项 | 值 |
|---|---|
| 环境路径 | `/home/echo/envs/chemprop` |
| 解释器 | `/home/echo/envs/chemprop/bin/python` |
| Python 版本 | Python 3.12.14 |
| 环境管理器 | micromamba `2.9.0`（root prefix `/home/echo/micromamba`） |
| 是否安装 conda | **否** —— 使用免 sudo 的单文件 micromamba |
| 系统 Python | Python 3.10.12（Ubuntu 自带，**未被改动**，与本环境隔离） |
| 已安装 pip 包数 | 102 |

导入自检：`ok`

---

## 四、关键组件版本

| 组件 | 版本 |
|---|---|
| `chemprop` | 2.3.1 |
| `admet_ai` | 2.0.1 |
| `torch` | 2.11.0+cu128 |
| `rdkit` | 2026.03.6 |
| `lightning` | 2.6.6 |
| `numpy` | 2.5.3 |
| `pandas` | 3.0.6 |
| `scipy` | 1.18.1 |
| `scikit_learn` | 1.9.1 |
| `tensorboard` | 2.21.0 |
| `flask` | 3.1.3 |
| `gunicorn` | 26.2.0 |
| `matplotlib` | 3.11.2 |
| `seaborn` | 0.13.2 |
| `astartes` | 1.3.3 |
| `descriptastorus` | 2.8.0 |
| `cuik_molmaker_pin` | 2026.3.6 |
| `myerson` | 1.0.1 |
| `openpyxl` | 3.1.5 |
| `typed_argument_parser` | 1.12.0 |

命令行入口自检：

```
$ chemprop --version
chemprop 2.3.1

$ admet_predict --help | head -3
usage: admet_predict --data_path DATA_PATH [--models_dir MODELS_DIR]
                     [--save_path SAVE_PATH] [--include_physchem]
                     [--drugbank_path DRUGBANK_PATH] [--atc_code ATC_CODE]
```

---

## 五、依赖清单（保证可复现）

### 5.1 `pip freeze`（环境内直接导出，权威）

> 把下面代码块另存为 `requirements.txt`，即可 `pip install -r requirements.txt` 复现。
> 共 102 个包。

```
absl-py==2.5.0
admet_ai==2.0.1
aimsim_core==2.2.3
aiohappyeyeballs==2.7.1
aiohttp==3.14.3
aiosignal==1.4.0
astartes==1.3.3
attrs==26.1.0
blinker==1.9.0
chemprop==2.3.1
click==8.5.0
cloudpickle==3.1.2
ConfigArgParse==1.7.7
contourpy==1.4.0
cuda-bindings==12.9.7
cuda-pathfinder==1.6.0
cuda-toolkit==12.8.1
cuik_molmaker_pin==2026.3.6
cycler==0.12.1
descriptastorus==2.8.0
dill==0.4.1
docstring_parser==0.18.0
et_xmlfile==2.0.0
filelock==3.32.3
Flask==3.1.3
fonttools==4.66.0
frozenlist==1.8.0
fsspec==2026.7.0
grpcio==1.84.0
gunicorn==26.2.0
idna==3.20
itsdangerous==2.2.0
Jinja2==3.1.6
joblib==1.6.0
kiwisolver==1.5.1
lightning==2.6.6
lightning-utilities==0.15.3
Markdown==3.11
markdown-it-py==4.2.0
MarkupSafe==3.0.3
matplotlib==3.11.2
mdurl==0.1.2
mhfp==1.9.6
mordredcommunity==2.0.7
mpmath==1.3.0
multidict==6.9.1
multiprocess==0.70.19
myerson==1.0.1
narwhals==2.26.0
networkx==3.6.1
numpy==2.5.3
nvidia-cublas-cu12==12.8.4.1
nvidia-cuda-cupti-cu12==12.8.90
nvidia-cuda-nvrtc-cu12==12.8.93
nvidia-cuda-runtime-cu12==12.8.90
nvidia-cudnn-cu12==9.19.0.56
nvidia-cufft-cu12==11.3.3.83
nvidia-cufile-cu12==1.13.1.3
nvidia-curand-cu12==10.3.9.90
nvidia-cusolver-cu12==11.7.3.90
nvidia-cusparse-cu12==12.5.8.93
nvidia-cusparselt-cu12==0.7.1
nvidia-nccl-cu12==2.28.9
nvidia-nvjitlink-cu12==12.8.93
nvidia-nvshmem-cu12==3.4.5
nvidia-nvtx-cu12==12.8.90
openpyxl==3.1.5
packaging @ file:///home/conda/feedstock_root/build_artifacts/bld/rattler-build_packaging_1785888127/work
padelpy==0.1.17
pandas==3.0.6
pandas-flavor==0.8.1
pillow==12.3.0
propcache==0.5.4
protobuf==7.36.2
psutil==7.2.2
Pygments==2.21.0
pyparsing==3.3.3
python-dateutil==2.9.0.post0
pytorch-lightning==2.6.6
PyYAML==6.0.3
rdkit==2026.3.6
rich==15.0.0
scikit-learn==1.9.1
scipy==1.18.1
seaborn==0.13.2
setuptools==78.1.0
six==1.17.0
sympy==1.14.0
tabulate==0.10.0
tensorboard==2.21.0
tensorboard-data-server==0.7.2
threadpoolctl==3.7.0
torch==2.11.0+cu128
torchmetrics==1.9.0
tqdm==4.70.1
triton==3.6.0
typed-argument-parser==1.12.0
typing_extensions==4.16.0
Werkzeug==3.1.8
wheel==0.48.0
xarray==2026.7.0
yarl==1.25.1
```

### 5.2 micromamba 安装的包（Python 解释器与系统库层）

```
List of packages in environment: "/home/echo/envs/chemprop"

  Name                      Version       Build                 Channel    
─────────────────────────────────────────────────────────────────────────────
  ConfigArgParse            1.7.7         pypi_0                pypi       
  Flask                     3.1.3         pypi_0                pypi       
  Jinja2                    3.1.6         pypi_0                pypi       
  Markdown                  3.11          pypi_0                pypi       
  MarkupSafe                3.0.3         pypi_0                pypi       
  PyYAML                    6.0.3         pypi_0                pypi       
  Pygments                  2.21.0        pypi_0                pypi       
  Werkzeug                  3.1.8         pypi_0                pypi       
  _openmp_mutex             4.5           20_gnu                conda-forge
  absl-py                   2.5.0         pypi_0                pypi       
  admet_ai                  2.0.1         pypi_0                pypi       
  aimsim_core               2.2.3         pypi_0                pypi       
  aiohappyeyeballs          2.7.1         pypi_0                pypi       
  aiohttp                   3.14.3        pypi_0                pypi       
  aiosignal                 1.4.0         pypi_0                pypi       
  astartes                  1.3.3         pypi_0                pypi       
  attrs                     26.1.0        pypi_0                pypi       
  blinker                   1.9.0         pypi_0                pypi       
  bzip2                     1.0.8         hda65f42_10           conda-forge
  ca-certificates           2026.7.22     hbd8a1cb_0            conda-forge
  chemprop                  2.3.1         pypi_0                pypi       
  click                     8.5.0         pypi_0                pypi       
  cloudpickle               3.1.2         pypi_0                pypi       
  contourpy                 1.4.0         pypi_0                pypi       
  cuda-bindings             12.9.7        pypi_0                pypi       
  cuda-pathfinder           1.6.0         pypi_0                pypi       
  cuda-toolkit              12.8.1        pypi_0                pypi       
  cuik_molmaker_pin         2026.3.6      pypi_0                pypi       
  cycler                    0.12.1        pypi_0                pypi       
  descriptastorus           2.8.0         pypi_0                pypi       
  dill                      0.4.1         pypi_0                pypi       
  docstring_parser          0.18.0        pypi_0                pypi       
  et_xmlfile                2.0.0         pypi_0                pypi       
  filelock                  3.32.3        pypi_0                pypi       
  fonttools                 4.66.0        pypi_0                pypi       
  frozenlist                1.8.0         pypi_0                pypi       
  fsspec                    2026.7.0      pypi_0                pypi       
  grpcio                    1.84.0        pypi_0                pypi       
  gunicorn                  26.2.0        pypi_0                pypi       
  icu                       78.3          py310h44b86e0_2       conda-forge
  idna                      3.20          pypi_0                pypi       
  itsdangerous              2.2.0         pypi_0                pypi       
  joblib                    1.6.0         pypi_0                pypi       
  kiwisolver                1.5.1         pypi_0                pypi       
  ld_impl_linux-64          2.46.1        default_hbd61a6d_102  conda-forge
  libexpat                  2.8.5         hd2095e1_0            conda-forge
  libffi                    3.7.0         h81df57d_1            conda-forge
  libgcc                    16.2.0        ha9f2e26_7            conda-forge
  libgomp                   16.2.0        he0feb66_7            conda-forge
  liblzma                   5.8.3         hb03c661_1            conda-forge
  libnsl                    2.0.1         hb9d3cd8_1            conda-forge
  libpython                 3.12.14       h0c77377_3_cpython    conda-forge
  libsqlite                 3.53.4        h13e7031_1            conda-forge
  libstdcxx                 16.2.0        h934c35e_7            conda-forge
  libuuid                   2.42.4        hcfc3c73_0            conda-forge
  libxcrypt                 4.4.38        h280c20c_0            conda-forge
  libzlib                   1.3.2         h25fd6f3_3            conda-forge
  lightning                 2.6.6         pypi_0                pypi       
  lightning-utilities       0.15.3        pypi_0                pypi       
  markdown-it-py            4.2.0         pypi_0                pypi       
  matplotlib                3.11.2        pypi_0                pypi       
  mdurl                     0.1.2         pypi_0                pypi       
  mhfp                      1.9.6         pypi_0                pypi       
  mordredcommunity          2.0.7         pypi_0                pypi       
  mpmath                    1.3.0         pypi_0                pypi       
  multidict                 6.9.1         pypi_0                pypi       
  multiprocess              0.70.19       pypi_0                pypi       
  myerson                   1.0.1         pypi_0                pypi       
  narwhals                  2.26.0        pypi_0                pypi       
  ncurses                   6.6           hdb14827_1            conda-forge
  networkx                  3.6.1         pypi_0                pypi       
  numpy                     2.5.3         pypi_0                pypi       
  nvidia-cublas-cu12        12.8.4.1      pypi_0                pypi       
  nvidia-cuda-cupti-cu12    12.8.90       pypi_0                pypi       
  nvidia-cuda-nvrtc-cu12    12.8.93       pypi_0                pypi       
  nvidia-cuda-runtime-cu12  12.8.90       pypi_0                pypi       
  nvidia-cudnn-cu12         9.19.0.56     pypi_0                pypi       
  nvidia-cufft-cu12         11.3.3.83     pypi_0                pypi       
  nvidia-cufile-cu12        1.13.1.3      pypi_0                pypi       
  nvidia-curand-cu12        10.3.9.90     pypi_0                pypi       
  nvidia-cusolver-cu12      11.7.3.90     pypi_0                pypi       
  nvidia-cusparse-cu12      12.5.8.93     pypi_0                pypi       
  nvidia-cusparselt-cu12    0.7.1         pypi_0                pypi       
  nvidia-nccl-cu12          2.28.9        pypi_0                pypi       
  nvidia-nvjitlink-cu12     12.8.93       pypi_0                pypi       
  nvidia-nvshmem-cu12       3.4.5         pypi_0                pypi       
  nvidia-nvtx-cu12          12.8.90       pypi_0                pypi       
  openpyxl                  3.1.5         pypi_0                pypi       
  openssl                   3.6.4         h781a0a9_0            conda-forge
  packaging                 26.3          pyhc364b38_0          conda-forge
  padelpy                   0.1.17        pypi_0                pypi       
  pandas                    3.0.6         pypi_0                pypi       
  pandas-flavor             0.8.1         pypi_0                pypi       
  pillow                    12.3.0        pypi_0                pypi       
  pip                       26.2.1        pyh8b19718_0          conda-forge
  propcache                 0.5.4         pypi_0                pypi       
  protobuf                  7.36.2        pypi_0                pypi       
  psutil                    7.2.2         pypi_0                pypi       
  pyparsing                 3.3.3         pypi_0                pypi       
  python                    3.12.14       h5f976f7_3_cpython    conda-forge
  python-dateutil           2.9.0.post0   pypi_0                pypi       
  pytorch-lightning         2.6.6         pypi_0                pypi       
  rdkit                     2026.3.6      pypi_0                pypi       
  readline                  8.3           hd6e31c0_1            conda-forge
  rich                      15.0.0        pypi_0                pypi       
  scikit-learn              1.9.1         pypi_0                pypi       
  scipy                     1.18.1        pypi_0                pypi       
  seaborn                   0.13.2        pypi_0                pypi       
  setuptools                84.0.0        pyh332efcf_0          conda-forge
  six                       1.17.0        pypi_0                pypi       
  sympy                     1.14.0        pypi_0                pypi       
  tabulate                  0.10.0        pypi_0                pypi       
  tensorboard               2.21.0        pypi_0                pypi       
  tensorboard-data-server   0.7.2         pypi_0                pypi       
  threadpoolctl             3.7.0         pypi_0                pypi       
  tk                        8.6.13        noxft_h1df4ec4_4      conda-forge
  torch                     2.11.0+cu128  pypi_0                pypi       
  torchmetrics              1.9.0         pypi_0                pypi       
  tqdm                      4.70.1        pypi_0                pypi       
  triton                    3.6.0         pypi_0                pypi       
  typed-argument-parser     1.12.0        pypi_0                pypi       
  typing_extensions         4.16.0        pypi_0                pypi       
  tzdata                    2026c         h151e31d_0            conda-forge
  wheel                     0.48.0        pyhd8ed1ab_0          conda-forge
  xarray                    2026.7.0      pypi_0                pypi       
  yarl                      1.25.1        pypi_0                pypi       
  zstd                      1.5.7         hb78ec9c_7            conda-forge
```
