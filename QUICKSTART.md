# 5 分钟跑通（QUICKSTART）

> 目标：**用最短路径亲眼看到这个流程能跑、结果能核对**。
> 完整说明见 [`README.md`](README.md)；本文只给「照抄就能跑」的最小集。

---

## 一、环境要求

| 项 | 要求 | 备注 |
| --- | --- | --- |
| 操作系统 | **Ubuntu 22.04（WSL2 亦可）** | 本项目实际构建环境 |
| Python | **3.10+** | 本机实测 3.10.12 |
| **GPU** | **仅「完整重跑」需要** | 演示与核对模式 **CPU 也能跑** |
| **gnina** | **必须单独装**（不在 PyPI） | 见 [README §2.3](README.md) |

**Python 依赖**（`requirements.txt`）：
```bash
pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple
```
> ⚠️ **务必装上 `openpyxl`** —— 缺它时读不了 `results/ADME预测_109个分子.xlsx`，
> 结果表会从 **23 列掉到 18 列**。为避免把不完整的表误当交付件，缺它时脚本会改用
> `results_无ADME_非交付件.csv` 输出，**不会覆盖交付件**。

---

## 二、四条命令（按需选）

```bash
cd project

# ① 环境自检：打印解析出的路径与参数（零依赖，先跑这个）
python3 src/config.py

# ② 演示：5 个分子 × LRR13–16 凹面，约 1 分钟，输出到 results/demo/
python3 predict.py --demo

# ③ 核对：复用包内已附带的对接结果，秒级重算候选清单与结果文件
python3 predict.py --report-only

# ④ 完整重跑：全库 109 个分子（需 GPU，本机实测约 66 分钟）
python3 predict.py
```

**⇒ 只想确认「能跑、能复现」的话，跑 ② 和 ③ 就够了** —— **不需要 GPU，总共 2 分钟。**

---

## 三、怎么判断「复现成功」

**跑完 ③ 之后**，对比 `results/results.xlsx`：

| 检查项 | 期望 |
| --- | --- |
| 工作表「候选清单」行数 | **109 行** |
| 列数 | **23 列** |
| **刺芒柄花素（001）** | **CNNscore 0.8350 · affinity −7.76 · 全库第 7 / 109** |
| 最终候选数 | **19** |

> **注意**：包内 `results/*.tsv`、`results/*.txt` 在 Linux 下生成为 **LF** 换行。
> 若在 **Windows** 上重跑，Python 会按平台写成 **CRLF**，于是这些文件**字节不同**
> （文本内容逐行一致，数值不受影响）。**出现该差异不视为复现失败**，以包内附带版本为准。

---

## 四、目录导航（想深入看哪里）

| 想知道 | 去哪 |
| --- | --- |
| 整个流程怎么设计 | [`docs/流程说明.md`](docs/流程说明.md) |
| 用了什么模型、有什么局限 | [`models/README.md`](models/README.md)（**Model Card**） |
| 数据从哪来、什么许可 | [`data/README.md`](data/README.md) |
| 结果每个字段什么意思 | [`results/README.md`](results/README.md) |
| 推荐区域是怎么定的 | [README §四](README.md)「**三条与对接无关的几何事实**」 |
| 排序逻辑与不确定性 | [README §五](README.md) |
| 已知局限 | [README §六](README.md) |
| 许可与第三方软件 | [README §七](README.md) · [`LICENSE`](LICENSE) |
| 复现性核对清单 | [README §八](README.md) |

---

## 五、常见问题

**Q：没有 GPU 能跑吗？**
能。`--demo` 和 `--report-only` 都可在 CPU 上跑；只有「完整重跑 109 个分子」建议用 GPU。

**Q：`predict.py` 报找不到 gnina？**
gnina 不在 PyPI，需单独编译安装。见 [README §2.3](README.md)。
**只想核对结果的话不必装 gnina** —— 直接跑 `--report-only` 即可。

**Q：为什么 `design.py` / `train.py` 是空的？**
**不是空的，是「说明型入口」**：本项目做的是**虚拟筛选**，不含生成式分子设计、也未训练模型。
两个文件各自说明了原因并指向真正的入口 `predict.py`。详见 [`models/README.md`](models/README.md)。

**Q：`--sites U3` 里的 `U3` 是什么？**
**历史内部键名**，现统一表述为 **LRR13–16 凹面**（受体残基 353–446）。
键名仅保留在**文件名与参数**中，见 [`config.yaml`](config.yaml) 的 `sites`。

**Q：随机种子是多少？**
**全部固定为 `42`**（`config.yaml` 的 `docking.seed` 与 `ligand_prep.seed`）。

---

## 六、引用

见 [`CITATION.cff`](CITATION.cff)（GitHub 仓库页有 **"Cite this repository"** 按钮）。
许可：**Apache-2.0**（[`LICENSE`](LICENSE)）。
