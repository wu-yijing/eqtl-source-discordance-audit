# -*- coding: utf-8 -*-
"""BMC 图件格式预检（只读）：页面尺寸 / 字体类型 / 最小字号 / PNG dpi。
用于 figure_scripts_officialZ_20260917 批次重出前后的对照。

判据（BMC / Springer Nature）：
  - 页宽 <= 170 mm（单栏 90 mm 亦可）；页高 <= 225 mm
  - 字体类型应为 Type0(TrueType 子集) / TrueType，**不应为 Type3**
  - 正文最小字号 >= 6 pt（在最终印刷尺寸下）
  - 位图 >= 300 dpi；本项目统一 600 dpi
用法：python 11_figure_precheck.py [--dir <图集目录> ...]
"""
import sys, os, glob

import fitz  # PyMuPDF
from PIL import Image
import paths_config as P  # 统一路径入口（2026-09-20）

MM = 25.4
W_LIMIT, H_LIMIT, PX_LIMIT = 170.0, 225.0, 6.0

DEFAULT_DIRS = [P.OUT_MAIN, P.OUT_SUPP]


def pdf_report(p):
    d = fitz.open(p)
    pg = d[0]
    w, h = pg.rect.width / 72 * MM, pg.rect.height / 72 * MM
    fonts = sorted({f[3].split('+')[-1] for f in pg.get_fonts(full=True)} or [])
    types = sorted({f[2] for f in pg.get_fonts(full=True) or []})
    spans = []
    for b in pg.get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                if (s['text'] or '').strip():
                    spans.append((s['size'], s['text'].strip()[:26]))
    spans.sort()
    mn = spans[0] if spans else (None, '')
    # 第二小（用于判断是否只有孤立的极小标注）
    mn2 = spans[1] if len(spans) > 1 else (None, '')
    n_small = sum(1 for s, _ in spans if s < PX_LIMIT)
    info = dict(pages=d.page_count, w=w, h=h, types=types, fonts=fonts,
                min_size=mn[0], min_text=mn[1], min2=mn2[0], min2_text=mn2[1],
                n_small=n_small, n_spans=len(spans))
    d.close()
    return info


def png_report(p):
    im = Image.open(p)
    dpi = im.info.get('dpi', (None, None))
    return dict(size=im.size, mode=im.mode, dpi=dpi)


def main(dirs):
    print('=' * 78)
    print('%-10s %-8s %-22s %-9s %-9s %s' % ('图', '载体', '页尺寸(mm) / 像素', '字体类型', '最小字号', '判定'))
    print('=' * 78)
    problems = []
    for dd in dirs:
        for stem in sorted({os.path.splitext(os.path.basename(f))[0]
                            for f in glob.glob(os.path.join(dd, '*.pdf'))
                            + glob.glob(os.path.join(dd, '*.png'))}):
            p = os.path.join(dd, stem + '.pdf')
            if not os.path.exists(p):
                continue
            r = pdf_report(p)
            flags = []
            if r['w'] > W_LIMIT + 0.05:
                flags.append('宽%.1f>170' % r['w'])
            if r['h'] > H_LIMIT + 0.05:
                flags.append('高%.1f>225' % r['h'])
            if any('Type3' in t for t in r['types']):
                flags.append('Type3')
            if r['min_size'] is not None and r['min_size'] < PX_LIMIT - 1e-6:
                flags.append('最小%.2fpt' % r['min_size'])
            verdict = 'OK' if not flags else '⚠ ' + ' / '.join(flags)
            if flags:
                problems.append((stem, flags))
            print('%-10s %-8s %-22s %-9s %-9s %s' % (
                stem, 'pdf',
                '%.1f x %.1f (%.2f:1)' % (r['w'], r['h'], r['w'] / r['h']),
                ','.join(r['types']) or '-',
                '%.2f [%s]' % (r['min_size'], r['min_text']) if r['min_size'] else '-',
                verdict))
            png = os.path.join(dd, stem + '.png')
            if os.path.exists(png):
                g = png_report(png)
                print('%-10s %-8s %-22s %-9s %-9s %s' % (
                    '', 'png', '%d x %d' % g['size'], g.get('mode', ''),
                    'dpi=%s' % (g['dpi'][0],),
                    '' if (g['dpi'][0] or 0) >= 300 else '⚠ dpi<300'))
            if r['n_small']:
                print('%22s  <6pt 的 span 数 = %d（次小 %.2f [%s]）'
                      % ('', r['n_small'], r['min2'], r['min2_text']))
    print('=' * 78)
    print('有问题的图：', problems if problems else 'NONE')
    return problems


if __name__ == '__main__':
    args = sys.argv[1:]
    ds = [args[args.index('--dir') + 1]] if '--dir' in args else DEFAULT_DIRS
    main(ds)
