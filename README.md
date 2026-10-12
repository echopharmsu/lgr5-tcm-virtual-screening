# LGR5(GPR49) 中药小分子虚拟筛选 · gnina

> **在线仓库（公开）**：<https://github.com/echopharmsu/lgr5-tcm-virtual-screening>

> **赛道**：赛道三 · 离子通道、GPCR等药靶的小分子药物虚拟筛选
> **靶点**：人源 **LGR5 / GPR49**（UniProt **O75473**），实验结构 RCSB PDB **4KNG**（X-ray 2.50 Å）
> **分子库**：10 味中药、**109 个**天然产物单体
> **方法**：gnina（CNN 打分 + 蒙特卡洛构象搜索）分子对接
> **主运行入口**：`predict.py`

本项目从**公共结构 4KNG** 出发，可一键重跑出全部结果：
配体 3D 化（RDKit ETKDGv3）→ gnina 虚拟筛选 → 四道过滤 → 候选清单与 `results/results.xlsx`。
全程固定随机种子 `--seed 42`；路径与参数统一在 `config.yaml`。

---

## 一、目录结构（附件5 第二章）

> **关于 `U1` / `U3` / `U10` 这三个键名**：早期版本曾并行考察三处候选区域，**该三位点框架已作废**。
> 本项目最终**只推荐单一区域** —— **LRR13–16 凹面**（残基 353–446），其确定过程见 **§四**。
> 下列**文件名与配置键**中的 `U1/U3/U10`、以及 `apo/holo`，均为**历史命名，为保持溯源性不予改动**；
> **正文叙述中的“位点”一律指这单一区域。**

```
project/
├── README.md              本文件：环境配置、运行命令、输入输出与结果说明
├── requirements.txt       Python 依赖及必要版本
├── config.yaml            ★ 全部路径与参数（改这里，不用改代码）
├── predict.py             ★ 主运行入口：一键生成最终结果文件
├── screen.py              （在 src/ 内）虚拟筛选入口：配体 3D 化 + gnina 对接
├── train.py               说明性占位：本项目未训练模型（见文件内注释）
├── design.py              说明性占位：本项目不含生成式设计（见文件内注释）
├── data/                  输入数据与来源说明
│   ├── README.md          数据来源、获取方式与许可
│   ├── 4KNG_full.pdb      4KNG 原始下载件（全部 6 条链）
│   ├── receptor_apo_LGR5.pdb    实际对接用受体（LGR5 A 链胞外域）
│   ├── receptor_holo_LGR5_RSPO1.pdb  含 RSPO1 的对照受体
│   ├── ligand_source.tsv  ★ 109 个配体（编号/名称/来源/SMILES）
│   ├── SMILES_修正版_110.txt    立体构型修正后的 SMILES 清单
│   └── ligands/           109 个配体的 2D SDF
├── src/                   核心源代码（含注释）
│   ├── config.py          路径与参数解析（禁止硬编码绝对路径）
│   ├── dock_prep.py       配体 3D 化（ETKDGv3 + MMFF）+ 平移到盒中心
│   ├── stereo_gate.py     立体化学闸门（未指定手性的分子：修复/枚举/拦截）
│   ├── stereo_audit.py    立体化学只读体检工具
│   ├── build_library.py   SMILES → 2D SDF 库
│   ├── library_build.py   建库 + 构型自检 + 归一化总表
│   ├── screen.py          虚拟筛选入口（对接）
│   ├── screen_analyze.py  逐分子指标、尺寸偏倚（共识列沿用历史键名）
│   ├── make_deliverables.py  四道过滤 → 候选清单
│   ├── make_results.py    ★ 生成最终结果文件 results.xlsx / results.csv
│   ├── merge_runs.py      多轮结果合并
│   ├── u3_analysis.py     位点残基接触与空腔几何专项分析
│   ├── u3_rnf43.py        位点与 RSPO1 / RNF43 界面的重合度
│   ├── pocket_scan.py     几何口袋检测（表面口袋 vs 内部空腔）
│   ├── pocket_global.py   同分辨率全局空腔对照
│   ├── prep_receptor.py   从 4KNG 抽链、去水去糖、altLoc 去重
│   ├── export_complex.py  导出小分子–靶点复合物结构
│   └── make_excel.py      研究用多工作表汇总（非最终结果文件）
├── docs/
│   └── 流程说明.md          ★ **把这条流水线移植到别的靶点**：哪些通用、哪些必须重做、踩过的坑
├── models/
│   └── README.md          ★ Model Card：所用模型、版本、参数、适用范围与局限
├── notebooks/
│   └── README.md          说明：本项目未附 Notebook
├── results/               ★ 运行结果与候选清单
│   ├── results.xlsx       ★ 最终结果文件（2 表：候选清单 / 说明）
│   ├── results.csv        同内容 UTF-8 版本（无 openpyxl 环境可直读）
│   ├── screen_候选清单_U1/U3/U10.tsv
│   ├── screen_unified_per_site.tsv / _consensus.tsv / screen_unified.txt
│   ├── screen_report.txt  候选清单 + 四道质检门禁
│   ├── screen_U1/U3/U10.sdf.gz   原始对接结果（各 109 分子 × 9 姿态）
│   ├── ADME预测_109个分子.xlsx    ADME（成药性/安全性）预测明细，编号口径与 results.xlsx 对齐
│   └── 结构文件/          小分子–靶点复合物 PDB
└── logs/                  gnina 原始日志（含命令行、随机种子、分数表）+ 随机种子记录
                            + ADMET-AI_Chemprop_环境与版本记录.md（第三方工具的环境记录）
```

---

## 二、环境配置（附件5 第一条、第三章）

### 2.1 操作系统与硬件（本项目的实际构建环境）

| 项 | 值 |
| --- | --- |
| 操作系统 | **Ubuntu 22.04.5 LTS**（WSL2，内核 6.18.33.2-microsoft-standard-WSL2） |
| CPU | Intel Core Ultra 5 225H（14 逻辑核心） |
| 内存 | 15 GiB |
| GPU | **NVIDIA GeForce RTX 5060 Laptop**（8151 MiB，驱动 591.74，CUDA 13.1） |
| 计算方式 | 全部计算在 WSL 侧完成；Windows 侧仅用于生成 Excel（`openpyxl`） |

### 2.2 Python 依赖

见 `requirements.txt`。核心：**RDKit 2026.03.6**、**NumPy 1.21.5**；Python **3.10.12**。

```bash
pip install -r requirements.txt
```

### 2.3 gnina 的安装

gnina **不在 PyPI 上**，需单独安装（本项目用的是源码编译版）：

| 方式 | 命令 / 说明 |
| --- | --- |
| 预编译（推荐，快） | 从 <https://github.com/gnina/gnina/releases> 下载 Linux 版并解压，把 `gnina` 放进 PATH |
| 源码编译（**本项目采用**） | master 分支 archive（`https://github.com/gnina/gnina/archive/refs/heads/master.zip`）；**commit `6fe1ce2bb9c35c8067de9f49bb9169857dfbad70`**；构建于 2026-09-11，CUDA 12.8，sm_120（Blackwell）；安装于 `/usr/local/bin/gnina` |
| 校验 | `gnina --version` → `gnina N/A N/A:N/A   Built Sep 11 2026.`（版本串取自 git，而源码包不含 `.git`，故为空；版本以 commit 号为准） |



---

## 三、运行（附件5 第一章、第二章、第七章）

### 3.1 一键跑通（主运行入口）

```bash
cd project

# ① 演示案例：5 个分子 × LRR13–16 凹面（内部键名 U3），约 1 分钟，输出到 results/demo/
python3 predict.py --demo

# ② 完整重跑：全库 109 个分子，本机实测约 66 分钟（GPU；三个对接盒各约 20 分钟）
python3 predict.py

# ③ 复用已对接结果，只重算候选清单与结果文件（秒级）
python3 predict.py --report-only

# ④ 只跑指定区域（参数名沿用历史键名）
python3 predict.py --sites U3
```

**主运行入口的每一步都会打印命令与耗时，并把原始输出写进 `logs/`。**

### 3.2 分步运行（便于排查）

```bash
python3 src/build_library.py data/ligand_source.tsv data/ligands   # 建 2D 库
python3 src/screen.py --site U3                                    # 单点对接
python3 src/screen_analyze.py data/ligand_source.tsv data/receptor_apo_LGR5.pdb \
        results/screen_unified results/screen_U1.sdf.gz results/screen_U3.sdf.gz \
        results/screen_U10.sdf.gz
python3 src/make_deliverables.py data/ligand_source.tsv data/receptor_apo_LGR5.pdb \
        results/screen results/screen_U1.sdf.gz results/screen_U3.sdf.gz \
        results/screen_U10.sdf.gz
python3 src/make_results.py results                                 # 出 results.xlsx
```

### 3.3 预期运行时间与资源

| 任务 | 实测耗时 | 资源 |
| --- | --- | --- |
| 建配体 2D 库（109 个） | ~5 s | CPU |
| 配体 3D 化（109 个；含氢最多 **219** 原子 / 重原子最多 **107**） | ~30 s | CPU |
| 单点对接（109 分子 × 9 姿态，exh 8） | 约 2–10 分钟 | GPU 约 1 GB 显存 |
| 演示案例（5 分子 × 1 点，exh 4） | **63 s**（实测） | GPU |
| 分析 + 出结果文件 | <1 分钟 | CPU |

所有随机性均固定：**`--seed 42`**（gnina 与 RDKit 3D 化各自统一为 42）。

---

## 四、输入与输出

### 输入

| 文件 | 说明 |
| --- | --- |
| `data/ligand_source.tsv` | 109 个配体：编号 / 中文名 / 英文名 / 来源中药 / SMILES |
| `data/receptor_apo_LGR5.pdb` | 受体：LGR5 A 链胞外域（残基 30–543，3616 重原子）。**由 4KNG 三元复合物抽 LGR5 链得到，非游离态结构**（文件名沿用历史） |
| `config.yaml` | 各对接盒中心、盒尺寸、exhaustiveness、过滤阈值等全部参数 |

### 输出

| 文件 | 说明 |
| --- | --- |
| **`results/results.xlsx`** | **最终结果文件**（附件5 第四章）。字段（与 `results.csv` 表头逐字一致）：候选编号 / 中文名 / 英文名 / 来源中药 / 所属赛道 / SMILES / 结构文件 / 推荐位点 / 推荐位点判定 / CNNscore / CNN预测亲和力 / 结合能(kcal/mol) / 距位点中心(A) / 接触残基数 / 氢键残基数 / 对应模型与运行版本 / 随机种子 / 血脑屏障安全度 / hERG安全度 / 无毒性评分 / 溶解度评分 / 生物利用度评分 / 备注 |
| `results/results.csv` | 同内容 UTF-8 版本 |
| `results/screen_候选清单_U*.tsv` | 各区域候选清单（历史键名 `U1/U3/U10`；含判定与未通过原因） |
| `results/screen_report.txt` | 候选清单 + 各区域的**几何统计（质检门禁）**。⚠️ 其中「**盒内姿态占比**」「**首姿态入盒分子数**」随随机种子与配体起跑位置明显变化，**不作性能指标**；该文件为**过程记录**，其统计结论**不用于对外表述**（与 `pocket_scan.py` 同一处理原则） |
| `results/screen_U*.sdf.gz` | 原始对接结果（109 分子 × 9 姿态） |
| **`results/blind_apo_e128.sdf.gz`** | **盲对接轮**（无先验信息、exhaustiveness 128 全蛋白扫描，40 个姿态）的原始结果 |
| **`results/focused_S1..S8.sdf.gz`** | **8 个聚焦盒**轮的原始结果（各 20 个姿态）。与上一行一起是 `src/u3_analysis.py` 的输入 —— **补上后 `results/U3_接触残基与氢键.tsv` 与 `results/U3_空腔性质.txt` 可在包内完整复算且逐字节重现**（已实测） |
| **`results/复现_8次运行/*.sdf.gz`** | **8 次复现运行的原始对接结果**（2 种起跑口径 × 4 个随机种子 `42/123/777/2024`；每次 109 分子 × 9 姿态）。**供评测方自行复算「哪些分子在所有条件下都通过」** —— 说明见该目录的 `README.md` |
| **`results/U3_接触残基与氢键.tsv`** | **报告的支撑数据**：13 个接触残基的距离 / 接触原子数 / 氢键数（最佳 CNN 与最佳 affinity 双列） |
| **`results/U3_RNF43界面关系.tsv`** | **报告的支撑数据**：与 RSPO1 / RNF43 界面的残基重合判定 |
| **`results/U3_空腔性质.txt`** | **报告的支撑数据**：位点局部自由球半径、包埋度、衬里残基 |
| **`results/U3_接触距离图.png`** | **报告用图**：13 个残基的接触距离与接触原子数（5 个几何氢键标红、ASP358 标为「非氢键」）。<br>⚠️ **该图副标题里的 CNNscore 0.8382 / affinity −7.75 取自「位点精算轮」（`logs/blind_apo_e128.log`），不是交付轮数字** —— 交付轮 001 为 **0.8350 / −7.76**（见 §5.2） |
| **`results/U3_site_local.png`** / **`U3_site_global.png`** | **报告用图**：位点局部视图（凹面 + 配体 + 13 个接触残基标注）/ 胞外域整体（含位点中心） |
| `results/ADME预测_109个分子.xlsx` | ADME（成药性/安全性）预测明细：52 个端点原始值与 DrugBank 百分位、风险标记 427 条。由本机以 **ADMET-AI（Chemprop 系，v2.0.1）** 实跑；**仅为成药性参考，不与对接分数合并或加权**。该工具的名称/版本/调用参数与完整依赖清单见 `logs/ADMET-AI_Chemprop_环境与版本记录.md` |
| `results/结构文件/*.pdb` | 小分子–靶点复合物 |
| `logs/*.log` | gnina / dock_prep 原始日志（含完整命令行） |
| **`logs/ADMET-AI_Chemprop_环境与版本记录.md`** | **第三方工具 ADMET-AI 的运行环境记录**（附件5 第一章/第三章要求）：系统 / GPU / CUDA / Python 环境 / 组件版本（chemprop 2.3.1、admet_ai 2.0.1、torch 2.11.0+cu128）+ 命令行入口自检 + **102 个包的完整 `pip freeze`**。由 `collect_env_md.py` **实跑采集、非人工转录** |
| **`logs/collect_env_md.py`** | **上表的采集脚本，随包提交** —— 可在配置环境中重跑核对（**不硬编码绝对路径**，输出路径由命令行参数或当前目录决定） |

### 推荐区域是怎么定下来的（**三条与对接无关的几何事实**）

本项目**不找“口袋”**。推荐区域的正当性由**三条不依赖任何对接结果的纯几何事实**判定 ——
它们**不看打分、不看姿态、不看筛选结果**，因此**不构成“因为分数好所以说区域对”的循环论证**。

| # | 几何事实 | 实测值 | 怎么算 |
| --- | --- | --- | --- |
| ① | **盒心到最近受体原子的距离** | **4.04 Å**（受体侧为 **ASP380 OD2**） | 遍历受体全部重原子取最小距离 |
| ② | **最近 6 个残基是否连续落在同一结构域内** | **ASP380 · ASN404 · ARG382 · SER402 · ILE379 · LYS378**（依次 4.04 / 4.09 / 4.42 / 4.51 / 4.53 / 4.72 Å），**全部落在 LRR13–16（353–446）内** | 按距离排序取前 6 |
| ③ | **与已知功能界面的距离与重合度** | 距 **RSPO1** 界面 **27.57 Å**；**13 个接触残基**与 RSPO1 / RNF43 界面的**重合均为 0** | 界面集取自 `results/U3_RNF43界面关系.tsv`（**LGR5(A)–RSPO1(M/P) 共 36 残基**，判据：任一重原子对距离 **< 5.0 Å**）；距离 = **区域中心到该残基集的最小重原子距离** |

**另有一条独立的结构学旁证**：全库 **13 个接触残基**中，**12 个属 LRR13–16（353–446）**，
第 13 个（**LYS449**）位于同一凹槽的 **C 端功能面一侧** ——
即该区域是**同一 LRR 凹槽的延续段**，而非一个孤立的新位点。

> 数据与判定过程见随包报告 **`results/报告_U3结合位点论证.md`**（文件名沿用历史键名）。

---

## 五、结果说明：排序逻辑与不确定性（附件5 第三章）

### 5.1 排序逻辑

判定为**四道过滤全部满足**才算通过（阈值见 `config.yaml` 的 `filters`）：

1. 最近受体接触 **≥ 2.5 Å**（排除原子穿透）
2. **affinity < 0** kcal/mol
3. **CNNscore ≥ 0.60**
4. **距位点中心 ≤ 6.0 Å**

通过者在前、未通过者在后，各自按 CNNscore 降序 —— 因此存在**两个不同的"名次"口径**，
必须并列阅读：

- **CNN 名次**：全部 109 个分子按 CNNscore 降序
- **候选表名次**：上述两段式排序内的名次

### 5.2 稳定性与不确定性分析

#### ① 两种配体起跑口径 —— **引用 001 的数字前必读**

001（刺芒柄花素）的数字**有两套，都对**，差别只在**配体从哪里开始跑**：

| | **基准起跑** —— 本结果文件用这套，**对外引用以此为准** ← **先读这一列** | 无冲突起跑（L1）—— 另一种起跑方式的对照；⚠️ **该次运行未随包** |
| --- | --- | --- |
| CNNscore | **0.8350** | 0.8473 |
| 结合能 | **−7.76 kcal/mol** | −7.75 kcal/mol |
| 全库 CNN 名次 | **第 7 / 109** | 第 6 / 109 |
| 起跑方式 | 109 个配体的质心**全部叠在区域同一点**起跑 | 配体被放到离区域中心 **8.31 Å（中位）之外**起跑 |

**两个口径下，001 都在其分子量箱（MW 250–300，n=16）内排第 1。**

**⇒ 对外引用一律以【基准起跑】（即本结果文件）为准：CNNscore 0.8350 / 结合能 −7.76 kcal/mol / 全库第 7。**
无冲突起跑（L1）那一轮起跑更严格、数字更好，但**其运行记录未随包提交**，因此不单独对外引用。

#### ② 通过闸门的计数不作性能指标

同一条件**只换随机种子**或**只换配体起跑位置**，通过数就会变：

| 配置 | 4 个随机种子的通过数 | 范围 |
| --- | --- | --- |
| 基准（质心叠盒心） | 19 / 13 / 13 / 19 | **13–19** |
| 无冲突起跑 | 10 / 19 / 14 / 21 | **10–21** |

**改问「哪些分子在所有条件下都通过」** —— 8 次运行（2 口径 × 4 种子）中：

| 配置 | 四种子共同通过者 |
| --- | --- |
| 基准 | 001 · 002 · 013 · 071 |
| 无冲突起跑 | 001 · 010 · 013 · 044 |
| **两口径交集（8 次全通过）** | **001 刺芒柄花素 · 013 大豆苷元** |

#### ③ 尺寸偏倚（n = 109）

| 相关性 | Spearman |
| --- | --- |
| 分子量 ↔ CNNscore | **−0.826**（CNN 偏小分子；源头 `results/screen_unified.txt:85`） |
| 分子量 ↔ 结合能 | **−0.712**（结合能偏大分子；U3 区域，源头 `results/screen_unified.txt:85`） |

库的分子量跨度 **126–1550** ⇒ 跨分子比较须限定分子量档位，本表同时列两个指标。

#### ④ 无黄金对照

本靶点**不存在**实验解析的小分子配体，无法做重对接 RMSD 验证，只能以几何事实与四道闸门替代。

> 完整论述见随包报告：`results/报告_U3结合位点论证.md`（文件名沿用历史）。

### 5.3 推荐结果（本项目结论）

- **推荐区域：LRR13–16 凹面**（残基 353–446；膜近端 7.9 Å）：13 个接触残基、
  **13 个接触残基与 RSPO1 / RNF43 界面的重合均为 0**；距 RSPO1 界面 **27.57 Å**。
（口径说明：**「无重合」限定在接触残基层面**。RNF43 界面的 19 个残基中有 10 个本就位于 LRR13–16（353–446）区间内，因此不宜表述为「整段区域与界面无重合」。）
- **推荐分子 001 刺芒柄花素（Formononetin，CAS 485-72-3）**：在推荐区域（LRR13–16 凹面）通过四道闸门，候选表名次 3（该轮通过 19 个；**计数不作性能指标**），
  在 **MW 200–450 同尺寸子集（n=43）** 内排第 1。

---

## 六、已知局限（附件5 第三章 Model Card 同步标注）

1. **受体为刚性**，且**仅胞外域**（残基 30–543），**不含 7TM 与胞内区**，**484–536 共 53 残基未建模**。
2. **本靶点无任何实验解析的小分子配体**（RCSB 全库 LGR5 结构 7 个，全部是胞外域，
   HETATM 仅 NAG/离子/水）→ **无法做重对接的黄金对照**，只能以几何门禁替代。
3. **目标区域为浅表面**：埋藏比约 **0.30**（常规深口袋在 0.6 以上）。
   > 注：`src/pocket_scan.py` 给出的「类药口袋」判定标准过严，其结论不用于对外表述。
4. **打分受分子尺寸影响**（见 5.2），大分子皂苷/三萜的负 affinity 主要是接触面积效应。
5. **CNNscore 不可跨轮或跨受体网格比较**，本表只做同一轮内排序。
6. **立体化学**：原始分子库中有 9 个分子的手性在 PubChem 源头记录里未指定，
   3D 化会任选构型并静默固化。本项目已用修正版构型重跑 11 个、剔除 1 个无法确定构型者
   （详见 `results/` 内相关说明）。
7. **结果为计算假说**，不代表已验证的结合位点或活性分子；未做 MD / MM-GBSA 精修，
   未做突变或竞争实验。

---

## 七、许可（附件5 第三条）

### 7.1 本项目自身的许可

**本项目代码以 Apache-2.0 许可发布** —— 许可全文见仓库根目录 **[`LICENSE`](LICENSE)**（Apache License 2.0，含完整条款与免责声明）。

> 关于 `LICENSE` 的措辞：根据 Apache-2.0 附录的实践，版权行按「项目名义」填写，**版权人署名为项目本身而非个人或单位**。

### 7.2 第三方软件与数据

本项目**不修改也不重新分发**下列第三方代码，仅按各自许可调用；完整版本与来源见下表。

| 软件 / 数据 | 版本 | 用途 | 许可 |
| --- | --- | --- | --- |
| **gnina** | 源码编译版（构建 2026-09-11，sm_120；**未内嵌版本号**） | 分子对接与 CNN 打分 | Apache-2.0 |
| **RDKit** | 2026.03.6 | 配体 3D 化、立体化学判定 | BSD-3-Clause |
| **NumPy** | 1.21.5 | 数值计算 | BSD-3-Clause |
| **openpyxl** | ≥3.0 | 生成 `results.xlsx` | MIT |
| **RCSB PDB 4KNG** | X-ray 2.50 Å | 受体结构 | PDB 数据 **CC0** |
| **PubChem** | 下载日期**未单独随包记录** | 配体结构与物性 | 公有领域 |
| **ADMET-AI** | **2.0.1** | ADME（成药性/安全性）预测，产出 `results/ADME预测_109个分子.xlsx` | **MIT**（取自环境内包元数据 `License-Expression`；来源 `https://github.com/swansonk14/admet_ai`） |
| **Chemprop** | **2.3.1** | ADMET-AI 的底层消息传递神经网络框架 | **MIT**（同上，环境内包元数据） |

调用参数：gnina 每次运行的完整命令行见 `logs/*.log` 首部（**日志内不含时间戳**）。
**ADMET-AI 为独立于对接流程的一次性调用**（产出时间见 `results/ADME预测_109个分子.xlsx` 文件时间戳）；
其**调用方式、模型版本、GPU/CUDA、Python 环境与 102 个依赖包的完整清单**见
**`logs/ADMET-AI_Chemprop_环境与版本记录.md`**（由 `collect_env_md.py` 实跑采集）。
**该预测不在对接流程内，不参与任何排序、筛选或判定，也不可与 gnina 分数比较或加权。**

---

## 八、复现性核对清单（供评测方）

```bash
python3 src/config.py          # 打印解析出的路径与参数，确认环境正确
python3 predict.py --demo      # 一键演示（约 1 分钟）
python3 predict.py --report-only   # 用已附带的结果重算，核对与 results.xlsx 是否一致
```

- 全部随机种子固定为 **42**（`config.yaml` 的 `docking.seed` 与 `ligand_prep.seed`）
- 路径与参数统一在 `config.yaml`
- `results/` 内的对接结果由 `predict.py` 同一协议产出，可用 `--report-only` 复算核对
- ⚠️ **前置依赖（重跑前必读）**：重跑 `--report-only` 前**必须装 `openpyxl`**（`pip install openpyxl`）。
  缺它时：① 生不出 `results.xlsx`；② **读不了 `results/ADME预测_109个分子.xlsx`，结果表会少 5 个 ADME 列（23 → 18 列）**。
  **为避免把不完整的表误当交付件，`make_results.py` 在缺 openpyxl 时会改用
  `results_无ADME_非交付件.csv` 输出，不会覆盖 23 列的交付件。**
- **换行符的平台差异（不影响内容）**：包内 `results/*.tsv`、`results/*.txt` 在 WSL/Linux 下生成，
  换行符为 LF。若在 **Windows** 上重跑 `--report-only`，Python 会按平台把换行写成 CRLF，
  于是这些文件与附带版本的**字节**不同（文件体积会大出「行数」个字节）。
  **文本内容逐行一致，数值不受影响**；出现该差异时以包内附带版本为准，不视为复现失败。
  在 WSL/Linux 下重跑则与附带版本逐字节相同。

### 7.3 引用方式与版本记录

- **引用信息**：见仓库根目录 **[`CITATION.cff`](CITATION.cff)**（Citation File Format 1.2.0）。
  在 GitHub 上会显示为仓库页的「Cite this repository」按钮。
- **版本记录**：见 **[`CHANGELOG.md`](CHANGELOG.md)**（按 Keep a Changelog 组织，
  版本号沿用交付日期标签；首个正式版本为 `v1.0`）。
- **版本管理**：本仓库使用 git 管理，`v1.0` 为区赛交付标签。
- **在线仓库**：<https://github.com/echopharmsu/lgr5-tcm-virtual-screening>（公开，含完整 git 历史与 `v1.0` 标签）。
- **许可**：Apache-2.0（见 [`LICENSE`](LICENSE)）；第三方软件与数据许可见 §7.2。
