# ⚠️ 本目录（data/processed/）已作废 —— 请勿用于出图或统计

**作废日期**：2026-09-17
**原因**：本目录下的 S-PrediXcan 输出产生于 **实现缺陷修正之前**。该版本的自研 S-PrediXcan 例程
（1）在多 SNP 模型的分母中遗漏了 σᵢ 表达方差因子，（2）PLINK 2-bit 基因型解码错误。
两个缺陷都会**系统性放大密模型基因的 |Z|**。

**证据（同基因同一统计量对比）**

| 基因（DR） | 本目录的值 | 官方 MetaXcan v0.8.1 的值 | 倍数 |
|---|---|---|---|
| TUBB（eQTLGen） | 48.516 | 11.8932 | ×4.1 |
| RNH1（eQTLGen） | 13.318 | 2.3064 | ×5.8 |
| RNH1（GTEx Nerve_Tibial） | 13.8246 | 2.6675 | ×5.2 |
| CKAP4（eQTLGen） | 4.7549 | 0.9874 | ×4.8 |
| CKAP4（GTEx Nerve_Tibial） | −3.2266 | −0.5622 | ×5.7 |
| DDX5（eQTLGen） | −3.8681 | −0.0049 | ×789 |

**受影响的下游产物（均已确认或需重绘）**

| 产物 | 读取的本目录文件 | 状态 |
|---|---|---|
| Figure 5（候选基因 dumbbell） | `eqtlgen_spredixcan_harmonized_results.csv`, `gtex_Nerve_Tibial_DR.csv` | ✅ 已重绘（2026-09-17） |
| Figure 7（RNH1 跨队列） | `eqtlgen_spredixcan_harmonized_results.csv`, `gtex_Nerve_Tibial_{DR,DN,DPN}.csv` | ✅ 已重绘（2026-09-17） |
| Figure 6（FDR 富集三面板） | `eqtlgen_spredixcan_harmonized_results.csv`（eQTLGen 侧） | ⏳ 待按同一流程重绘 |
| Figure 8（跨性状三基因集） | `gtex_Whole_Blood_{t}.csv`, `gtex_Nerve_Tibial_{t}.csv` | ⏳ 待按同一流程重绘 |
| Figure 2 / 4（协变量平衡 / 双轴分解） | 需确认是否读取本目录 | ⏳ 待确认 |
| Figure S2 / S3 / S5 | `_data_testbed_eqtlgen.csv`（同源） | ✅ 已重绘（2026-09-17） |
| 主稿 Table 4 的 RNH1 / CKAP4 行 | 同源 | ✅ 已更正（2026-09-17） |
| Additional file 1 Table S4 的 √N_e 行 | 旧 Z 对导出 | ✅ 已更正（2026-09-17） |

**请改用**：`data/processed_officialZ/`（由官方 MetaXcan v0.8.1 结果导出，与稿件 SI 中的
Table S2 / S12 / S14 / S17 / Table S4 逐值一致）。

**建议**：出图脚本一律从 `data/processed_officialZ/` 读数；本目录保留仅作历史对照，
并在下次归档（Zenodo/GitHub）时于 README 中标注为 superseded。
