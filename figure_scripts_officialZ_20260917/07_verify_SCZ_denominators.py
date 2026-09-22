# -*- coding: utf-8 -*-
"""验证 SCZ 四臂官方数据（分母与 ρ），并纳入官方数据层"""
import csv, os, shutil
import numpy as np
from scipy import stats
import paths_config as P  # 统一路径入口（2026-09-20）

SRC = P.need(os.path.join(P.DATA_Z, 'scz_z_4arm_official.csv'), 'SCZ 四臂官方数据层')
rows = list(csv.DictReader(open(SRC, encoding='utf-8-sig')))
def arr(k):
    return np.array([float(r[k]) if r[k] not in ('', 'NA') else np.nan for r in rows])
eq, wb, nt, mt = arr('eqZ'), arr('wbZ'), arr('ntZ'), arr('multiZ')
print('genes %d | eq %d  wb %d  nt %d  multi %d' % (len(rows), (~np.isnan(eq)).sum(),
      (~np.isnan(wb)).sum(), (~np.isnan(nt)).sum(), (~np.isnan(mt)).sum()))
cc = ~(np.isnan(eq) | np.isnan(wb) | np.isnan(nt) | np.isnan(mt))
tw = ~(np.isnan(wb) | np.isnan(nt))
pw = ~(np.isnan(eq) | np.isnan(wb))
print('complete-case (4 arms) n = %d   <-- 稿件称 8,315' % cc.sum())
print('tissue-only available n  = %d   <-- 稿件称 8,890' % tw.sum())
print('panel-only available n   = %d   <-- 稿件称 9,048' % pw.sum())

def rep(lab, x, y, m):
    a, b = x[m], y[m]
    r = stats.spearmanr(a, b)[0]
    s = ((a > 0) == (b > 0)).mean() * 100
    print('  %-32s n=%5d  rho=%+.4f  consistency=%.1f%%' % (lab, len(a), r, s))

print('\n=== 稿件 §3.6 声称：panel 0.469 / tissue 0.420 / dual 0.447 ===')
print('complete-case (n = 8,315 宇宙):')
rep('panel-only (eqZ vs wbZ)', eq, wb, cc)
rep('tissue-only (wbZ vs ntZ)', wb, nt, cc)
rep('dual (eqZ vs multiZ)', eq, mt, cc)
print('available-case:')
rep('panel-only', eq, wb, pw)
rep('tissue-only (n = 8,890)', wb, nt, tw)

# 经验零假设（min|Z| < 0.5）复算
mn = np.minimum(np.abs(eq), np.abs(wb))
m05 = (mn < 0.5) & pw
s = ((eq > 0) == (wb > 0))
print('\n经验零假设（panel-only, min|Z| < 0.5）: %d/%d = %.1f%%   <-- 稿件称 54.9%% (1,987/3,617)'
      % (s[m05].sum(), m05.sum(), 100 * s[m05].mean()))

print('\n（本脚本只读归档层 %s；旧版会从会话临时目录复制到此处，2026-09-20 起不再需要）' % SRC)
