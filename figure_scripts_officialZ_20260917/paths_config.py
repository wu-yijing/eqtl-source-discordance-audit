# -*- coding: utf-8 -*-
"""officialZ 出图管线的**统一路径入口**（2026-09-20）。

所有脚本只从这里取输入/输出路径，不再各自写死绝对路径。
默认值按"脚本所在目录的上一级 = 仓库根"推导，因此**克隆仓库后可直接运行**；
在其他目录结构下用环境变量覆盖，无需改代码：

| 环境变量 | 含义 | 默认 |
|---|---|---|
| `TWAS_REPO`    | 仓库根（含 `data/`） | 脚本目录的上一级 |
| `TWAS_DATA_Z`  | 官方 MetaXcan Z 数据层 | `<TWAS_REPO>/data/processed_officialZ` |
| `FIG_OUT_MAIN` | 正文图输出目录（Fig1–Fig8，**运行时输出，已 gitignore**） | `<TWAS_REPO>/figures` |
| `FIG_OUT_SUPP` | 补充图输出目录（FigS1–FigS6，同上） | 同 `FIG_OUT_MAIN` |
| `AF1_DOCX`     | Additional file 1 的 .docx | `<TWAS_REPO>/additional_file_1/Additional file 1.docx` |
| `FIG_RESULTS`  | 中间结果与随包输入（json）目录 | 本脚本目录 |

⚠️ **Additional file 1 不从本仓库分发**（仓库不含 docx）。请从期刊补充材料下载后
用 `AF1_DOCX` 指向它，例如：

```bash
AF1_DOCX="/path/to/Additional file 1.docx" python 00_build_officialZ_data_layer.py
```
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('TWAS_REPO') or os.path.dirname(HERE)
DATA_Z = os.environ.get('TWAS_DATA_Z') or os.path.join(REPO, 'data', 'processed_officialZ')
OUT_MAIN = os.environ.get('FIG_OUT_MAIN') or os.path.join(REPO, 'figures')
OUT_SUPP = os.environ.get('FIG_OUT_SUPP') or OUT_MAIN
AF1 = os.environ.get('AF1_DOCX') or os.path.join(REPO, 'additional_file_1', 'Additional file 1.docx')
RES = os.environ.get('FIG_RESULTS') or HERE


def need(path, what='输入'):
    """存在性断言：缺失时给出可操作的提示，而不是抛 traceback。"""
    if not os.path.exists(path):
        sys.exit(
            '[路径缺失] %s 不存在：\n  %s\n'
            '  覆盖方式（任选其一）：\n'
            '    · 设 TWAS_REPO 指向你的仓库克隆（默认取脚本目录的上一级）\n'
            '    · 设 AF1_DOCX 指向从期刊下载的 Additional file 1.docx\n'
            '    · 设 FIG_OUT_MAIN / FIG_OUT_SUPP 指向你的图片输出目录'
            % (what, path))
    return path


def report():
    for k in ('REPO', 'DATA_Z', 'OUT_MAIN', 'OUT_SUPP', 'AF1', 'RES'):
        print('  %-9s = %s' % (k, globals()[k]))


if __name__ == '__main__':
    report()
