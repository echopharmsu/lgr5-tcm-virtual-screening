# 变更记录

本项目遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/) 的组织方式；
版本号沿用**交付日期标签**（`YYYYMMDD`），首个正式版本记为 `v1.0`。

## [v1.0] - 2026-10-12

### 新增
- `CITATION.cff` —— 引用信息（CFF 1.2.0）
- `CHANGELOG.md` —— 本文件
- `.gitignore`
- **建立 git 版本管理**：初始化仓库、首次提交、打标签 `v1.0`

### 变更
- `README.md` 增补「引用方式与版本记录」小节（§7.3）
  - 顶部与 §7.3 补入**在线仓库地址** <https://github.com/echopharmsu/lgr5-tcm-virtual-screening>

## [2026-10-11] - 2026-10-11

### 新增
- `results/复现_8次运行/` —— 8 次运行的原始对接结果与说明（2 种起跑口径 × 4 个随机种子）
- `results/blind_apo_e128.sdf.gz`、`results/focused_S1..S8.sdf.gz` —— 探索轮与聚焦轮的原始姿态

### 变更
- `README.md` 增补**姿态口径**说明（第 1 条姿态 = CNNscore 最高，非 affinity 最低）
- `README.md` 增补**敏感性披露**：2 vs 9 全通过者的差异
- `README.md` 统一「距 RSPO1 界面」口径为 **27.57 Å**（原 28.06 Å 为单链口径，已作废）
- 推荐区域表述统一为 **LRR13–16 凹面**（原内部键名 `U3` 仅保留在文件名与参数中）

## [2026-09-30] - 2026-09-30

### 新增
- 首次打包交付：109 个中药小分子库、配体准备与对接脚本、`results.xlsx` / `results.csv`
- 四道统一筛选判据（最近接触 ≥2.5 Å · 结合能 <0 · CNNscore ≥0.60 · 距区域中心 ≤6 Å）
- `predict.py` 主运行入口（含 `--demo` 演示模式）
- Apache-2.0 许可与第三方软件/数据许可表

## [公开] - 2026-10-12

### 新增
- **首次公开发布**：仓库 <https://github.com/echopharmsu/lgr5-tcm-virtual-screening>（Public · Apache-2.0 · 含 `v1.0` 标签）
