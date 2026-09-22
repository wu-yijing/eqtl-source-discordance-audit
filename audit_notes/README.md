# Audit notes

Dated audit notes **shipped with this repository**. They record cross-carrier consistency checks and
deprecation decisions that bear directly on the data archived here.

They are deliberately kept separate from `analysis_reports/`, which is a **local-only working
directory** and is listed in `.gitignore` (internal working notes — not published).

> **Read the notes as of their own date.** The notes below were written on 2026-09-20, before this
> repository was rebuilt (see the addendum at the end). They therefore quote directory names and
> repository paths **as they stood at the time** — e.g. a `figures/` directory, `_DEPRECATED_…`
> script directories, and the former repository name. Those paths are historical references and no
> longer resolve; they are left unedited on purpose so the audit chain stays traceable.

| Note | Date | What it records |
|---|---|---|
| `AF1与主稿图表同步核查报告_20260920.md` | 2026-09-20 | *Additional file 1* vs the manuscript: figure/table numbering, captions, cross-references and values all in sync; one **duplicated figure** was found and removed, with before/after verification at both the docx and PDF layers |
| `补充图目录版本谱系判定_20260920.md` | 2026-09-20 | Which supplementary-figure directory is current, decided by **content hash** rather than timestamps: the current directory was byte-identical to the images embedded in *Additional file 1* (6/6, vs 0/16 for the superseded directory). The superseded directory was given a `_DEPRECATED_勿用_` prefix |
| `S1重建与SCZ改名收口执行记录_20260920.md` | 2026-09-20 | The record behind commits `d7de5a8` and `14d2289`: the in-house SCZ pipeline retired from the distributed tree, the S1 section rebuilt onto `data/processed_officialZ/`, and the manuscript title unified across every carrier (README, `.zenodo.json`, Zenodo record) |

No `.docx` manuscript or supplementary file is distributed with this repository; the notes refer to
them only by file name and by values already public in `README.md`.

> 本目录的文件以中文撰写，与项目内部的审计记录体例一致。图集与数值的权威来源仍是
> `data/processed_officialZ/`；`README.md` 是这些结论的英文摘要。

---

## Addendum — 2026-09-23: repository rebuilt as a code-and-data-only release

This repository was rebuilt from a snapshot of its predecessor's `main` tree. The rebuild was a
**scope change**, not an analysis change: no Z-score, denominator, endpoint or reported value was
recomputed or altered. What changed:

| Change | Detail |
|---|---|
| **Figure and table artefacts removed** | All figure images (PNG/PDF/SVG), figure captions, table files (CSV/XLSX/DOCX) and table legends were dropped from the distributed tree. They belong to the manuscript and to *Additional file 1*, and are not mirrored here. `figures/`, `tables/`, `figs/`, `*_caption.txt` and `*_legend.txt` are now `.gitignore`d, so re-running the figure pipeline cannot re-introduce them. |
| **Superseded script directories removed** | The three `_DEPRECATED_*` directories (the pre-2026-09-17 five-script figure set, the old supplementary-figure scheme, and the in-house SCZ TWAS estimator) were dropped from the distributed tree. Their findings are still recorded in the notes above and in `README.md`. |
| **Two inherited defects repaired** | (i) Four files in `scripts/python/s1_cluster_robustness/` carried a mangled line — `if nxt == p:/n            break` — in which the newline between `p:` and `break` had been eaten by an earlier edit, making them unparseable. The source repository's `HEAD` blob was verified to contain the same defect, so it predates this rebuild; the obvious intent (two lines) was restored and all 61 Python files in this repository now compile. (ii) `s1_scz_cluster.py` pointed at a results file inside a now-removed directory; the hard-coded path was replaced with an env-overridable input (`S1_SCZ_SNAPSHOT`) and the fact that the snapshot is **not distributed** is stated at the point of use. |
| **Repository name in paths and metadata normalised** | Absolute paths inside the early-generation scripts, the `Dockerfile` label, `data/README.md` and `data/processed/gigadb_metadata_form.csv` were updated to the current repository name. This was a pure string substitution — no path semantics changed, and no script behaviour changed. |
| **Historical records intentionally not rewritten** | This file and the dated notes above, and `analyses/logs/*` (run logs), still show the former repository name and the paths valid when they were written. Rewriting them would break the audit chain. |

Verification performed at the rebuild: all 61 Python files compile; zero figure/table/caption/legend
files remain in the tree; no `.docx` is present.
