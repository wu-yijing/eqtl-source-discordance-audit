# figure_scripts_officialZ_20260917

**本目录是 2026-09-17 P0 修复后唯一的出图脚本来源。**
2026-09-20 起所有路径收敛到 `paths_config.py`（不再有写死的绝对路径），并补入了两个随包输入/产物 json。

## 为什么需要它

原 `figure_scripts/` 脚本集读取 `data/processed/*.csv`，
而该目录产生于**实现缺陷修正之前**（缺 σᵢ 表达方差因子 + PLINK 2-bit 解码错误），会系统性放大密模型基因的 |Z|：

| 基因 | 陈旧层值 | 官方 MetaXcan v0.8.1 值 |
|---|---|---|
| TUBB（eQTLGen, DR） | 48.52 | 11.89 |
| RNH1（eQTLGen, DR） | 13.32 | 2.31 |
| RNH1（GTEx Nerve_Tibial, DR） | 13.82 | 2.67 |
| CKAP4（eQTLGen, DR） | 4.75 | 0.99 |

后果：Fig. 3 / 4 / 5 / 6 / 7 / 8 的**图形**与稿件**正文数值**不一致——正文早已改用官方
MetaXcan v0.8.1 的 Z，图却没有跟着重出。本目录的脚本修复了这一点。
**旧脚本仍能运行，但会复现被取代的数值，请勿用于复现本文。**

## 路径入口（2026-09-20 起）

所有脚本只从 `paths_config.py` 取路径，**克隆仓库后可直接运行**；默认按"脚本目录的上一级 = 仓库根"推导：

| 环境变量 | 含义 | 默认 |
|---|---|---|
| `TWAS_REPO` | 仓库根（含 `data/`、`figures/`） | 脚本目录的上一级 |
| `TWAS_DATA_Z` | 官方 MetaXcan Z 数据层 | `<TWAS_REPO>/data/processed_officialZ` |
| `FIG_OUT_MAIN` | 正文图输出目录（Fig1–Fig8） | `<TWAS_REPO>/figures` |
| `FIG_OUT_SUPP` | 补充图输出目录（FigS1–FigS6） | 同 `FIG_OUT_MAIN` |
| `AF1_DOCX` | Additional file 1 的 `.docx` | `<TWAS_REPO>/additional_file_1/Additional file 1.docx` |
| `FIG_RESULTS` | 随包 json（输入/产物）目录 | 本脚本目录 |

> ⚠️ `FIG_OUT_MAIN` / `FIG_OUT_SUPP` 是**运行时输出目录**：脚本会把重出的图写进
> `<TWAS_REPO>/figures/`。该目录**不随本仓库分发**（已列入 `.gitignore`）——
> 图件属手稿与 *Additional file 1*，不在此仓库内。重跑后请勿 `git add`。

跑之前可以先看一遍解析结果：

```bash
python paths_config.py
#   REPO      = .../eqtl-source-discordance-audit
#   DATA_Z    = .../data/processed_officialZ
#   OUT_MAIN  = .../figures
#   OUT_SUPP  = .../figures
#   AF1       = .../additional_file_1/Additional file 1.docx
#   RES       = .../figure_scripts_officialZ_20260917
```

⚠️ **Additional file 1 不从本仓库分发**（仓库不含 docx，遵守"正文与补充材料正文不随代码仓库发布"的规则）。
从期刊补充材料下载后指定：

```bash
AF1_DOCX="/path/to/Additional file 1.docx" python 00_build_officialZ_data_layer.py
```

## 执行顺序

1. `00_build_officialZ_data_layer.py` — 从 Additional file 1 的官方表（S1/S2/S4/S12/S14/S17）导出 `data/processed_officialZ/`
2. `01_redraw_Fig5_Fig7.py` — 候选基因 dumbbell + RNH1 跨队列（读 `data/processed_officialZ/`）
3. `02_redraw_Fig3.py` — 头条散点图（96 对）
4. `04_redraw_Fig8.py` — 跨性状三基因集
5. `06_redraw_Fig4.py` — 轴分解 + Δρ 配对 bootstrap（B = 5,000, seed 20260915；结果写回 `fig4_bootstrap_officialZ.json`）
6. `08_redraw_Fig6_labels_20260920.py` — **Fig. 6 的权威脚本**（含 SI Table S18 逐格回归断言）
7. `10_redraw_FigS6_20260920.py` — Fig. S6（读 `m15_positive_control.json`；含 SI Table S19 逐点回归断言）
8. `05_recompute_arms.py` — 核验 Table 3b 与 45 基因共同宇宙（只读、只打印）
9. `07_verify_SCZ_denominators.py` — 核验 SCZ 分母与经验零假设（只读、只打印）
10. `11_figure_precheck.py` — **投稿格式预检（只读，可重复运行）**：遍历两个图集目录，逐张报出页宽(mm) / 字体类型 / 最小字号 / PNG 尺寸·色彩模式·dpi，并按 BMC 判据（页宽 ≤170 mm、最小字号 ≥6 pt、字体非 Type3、PNG 为 RGB 且 ≥300 dpi）给出 `OK / ⚠`。**每次重出图件后必跑。**
    ```bash
    python 11_figure_precheck.py            # 默认扫 OUT_MAIN + OUT_SUPP
    python 11_figure_precheck.py --dir <某目录>
    ```

### 已弃用（保留仅作审计）

**护栏的粒度 = 被取代的粒度**（两种情形不要混用同一种做法）：

| 脚本 | 被取代的范围 | 护栏 | 为什么 |
|---|---|---|---|
| `03_redraw_Fig6.py` | **整个脚本**（其唯一输出 Fig. 6 已由 `08` 重新出具） | **硬弃用**：默认 `SystemExit`，退出码 1，需 `FIG6_FROM_03=1` 强跑。提示里直接给出应运行的命令 | 运行它没有任何意义（不产出任何文件）；"跑完什么都不做还返回 0"会让批量运行把弃用误判为成功 |
| `02_redraw_Fig3.py` 内的 Fig6 段 | **仅一路输出**（脚本的主职责 Fig. 3 仍然有效） | **只跳过该输出**：打印 `SKIP Fig6 …` 后继续，退出码 0，需 `FIG6_FROM_02=1` 强跑 | 脚本还有有效职责，不应因一路输出被取代而整体失败 |

> ⚠️ 因此**官方执行顺序里 `03` 的退出码 1 是设计如此，不是故障**：Fig. 6 请由 `08` 生成。`08` 已在该顺序中（第 6 步），照顺序跑不会缺图。
> **不要**去掉这两个护栏；`02` 的 Fig6 段若被强跑会覆盖正确版。

其余一次性脚本：`_recon_Fig6.py`（对账，`_` 前缀，不随包发布）。

## 随包的 json（`FIG_RESULTS`，默认本目录）

| 文件 | 角色 | 说明 |
|---|---|---|
| `m15_positive_control.json` | **输入** | 端点校准与阳性对照的原始结果（BH 检测边界、零校准、单基因/组间 power）。`10_redraw_FigS6_20260920.py` 全部绘图值取自此处，并与 SI Table S19 逐点断言。 |
| `fig4_bootstrap_officialZ.json` | **产物** | Fig. 4(b) 的配对基因簇 bootstrap 结果（Δρ 点估计/95% CI/P）。由 `06` 确定性重算（B = 5,000，seed 20260915）并覆盖写回。 |

> 这两个文件原先只存在于会话工作目录（`2026-09-1x`），2026-09-20 起随脚本目录归档，
> 使 Fig. 4 与 Fig. S6 在**只有本仓库/本包**的情况下也能复现。

## 数据纪律

所有脚本**只读** `data/processed_officialZ/` 与 Additional file 1 的表格，
**不读** `data/processed/`。该陈旧目录已加弃用声明（`_DEPRECATED_勿用_修正前数据_20260917.md`），
保留仅作历史对照。

## 复现结果（2026-09-17 核验，与稿件逐值一致）

| 产物 | 复现值 | 状态 |
|---|---|---|
| Fig. 3 | 一致率 68.8%，ρ = 0.3896，max absZ = 2.83 | 一致 |
| Fig. 4 | Table 3b 三臂 ρ = 0.490 / 0.414 / 0.442；Δρ(dual−panel) = −0.020 [−0.127, +0.090]，P = 0.76；Δρ(dual−tissue) = +0.034 [−0.214, +0.280]，P = 0.78 | 一致 |
| Fig. 5 | TUBB eQTLGen +11.89（原图 +48.5）；CKAP4 +0.99（原图 +4.75） | 一致 |
| Fig. 6 | GTEx housekeeping 逐表型 4/29 + 0/29 + 1/29 = 5/87（= Table S7） | 一致 |
| Fig. 7 | RNH1 GTEx Nerve_Tibial +2.67 / +1.70 / +1.63 | 一致 |
| Fig. 8 | 0.414 (n = 138) / 0.636 (n = 72) / 0.418 (n = 8,890) | 一致 |
| SCZ 分母 | 基因池 10,357 / panel 9,048 / tissue 8,890 / 四臂 complete-case 8,315 | 一致 |

## 2026-09-20：路径重构与随包输入（本轮）

| 项 | 内容 |
|---|---|
| 路径 | 12 个脚本的 **10 处硬编码绝对路径 + 3 处会话临时目录依赖**全部收敛到 `paths_config.py`；`07` 的输入改指仓库内 `data/processed_officialZ/scz_z_4arm_official.csv`（与旧临时文件**逐字节相同**） |
| 随包输入 | `m15_positive_control.json`、`fig4_bootstrap_officialZ.json` 移入本目录并随包/随仓库发布 |
| AF1 指向 | 由已失效的 `Additional file 1.docx` / 被取代的 `…_20260917.docx` 改为可配置的 `AF1_DOCX`（本次核验用 `…_20260920.docx`；两版之间仅 S15/S16 变动，脚本所读表未变） |
| 验证 | 用 `AF1_DOCX=…0920.docx` + 定稿图集/补充图作为输出目录**实跑整条管线**：7 张图（Fig3/4/5/6/7/8/S6）重出后与重构前**逐像素 0 差异**（PNG mean\|Δ\| = 0.0000），PDF 文本层相同 |
| 连带修正 | `00` 重建数据层时发现仓库内 `gtex_official_Z.csv`（q 列 217/222 行）与 `crosscohort_TableS4_official.csv`（pooled Z、Q、I²、行标签）**与 AF1 Table S2/S4 不符**（陈旧世代）。重建后**逐格等于 AF1**（0 处不符）；图件不受影响。 |

## 2026-09-20 格式批次（只改格式参数，未改数据）

`figsize` 宽度 / 最小字号 / `pdf.fonttype` / PNG 色彩模式，一次处理 7 张不合规图（Fig3/4/5/6/7/S4/S6）+ Fig1 的 P2-4 限定语。
细则与前后对照见 `BMC Genomics投稿资料\定稿图集_Fig1-8_20260914\图件重命名映射与核对清单.md` 第 0.5 节。

## 未覆盖项

- **Fig. 2**（协变量平衡 SMD 图）：绘制的是基因长度 / GC 含量 / eQTL SNP 数，不含 Z 值，
  不受本次问题影响；其生成脚本不在本目录（`定稿资料` 侧的 `gen_figs4.py` 系列）。
- Fig. 1（流程图）与 Fig. S1–S5 各有自己的生成脚本，不在本目录；本目录覆盖 Fig. 3–8 与 Fig. S6。
