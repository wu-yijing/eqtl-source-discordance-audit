# -*- coding: utf-8 -*-
"""P0 Step 8：重绘 Figure 4（轴分解）——(a) 三臂 ρ 用 Table 3b 官方值；(b) Δρ 用 45 基因共同宇宙 + 新算配对 bootstrap"""
import paths_config as P  # 统一路径入口（2026-09-20）
import os, csv, shutil, json
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NEW = P.need(P.DATA_Z, '官方 MetaXcan Z 数据层')
OUT = P.OUT_MAIN
BK = os.path.join(OUT, '_backup_before_officialZ_redraw_20260917')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.7,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
C_PAN, C_TIS, C_DUA = '#2471A3', '#1E8449', '#C0392B'
rng = np.random.default_rng(20260915)

def rd(n):
    return list(csv.DictReader(open(os.path.join(NEW, n), encoding='utf-8-sig')))
f = lambda x: (float(x) if x not in (None, '', 'NA') else np.nan)
GT = {(r['Gene'], r['Trait']): r for r in rd('gtex_official_Z.csv')}
EQ = {(r['Gene'], r['Trait']): r for r in rd('eqtlgen_official_Z.csv')}
Ze = lambda k: f(EQ[k]['Z_eQTLGen']) if k in EQ else np.nan
Zw = lambda k: f(GT[k]['Z_Whole_Blood']) if k in GT else np.nan
Zn = lambda k: f(GT[k]['Z_Nerve_Tibial']) if k in GT else np.nan
Zm = lambda k: f(GT[k]['Z_multi_tissue']) if k in GT else np.nan

# ---- 共同宇宙 45 基因 / 135 对 ----
allk = sorted(set(GT) | set(EQ))
common = [k for k in allk if not any(np.isnan(v) for v in (Ze(k), Zw(k), Zn(k), Zm(k)))]
genes = sorted(set(k[0] for k in common))
print('共同宇宙: %d 对 / %d 基因' % (len(common), len(genes)))
def rho(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    return stats.spearmanr(x[m], y[m])[0]
ZK = np.array([[Ze(k), Zw(k), Zn(k), Zm(k)] for k in common])
GK = np.array([k[0] for k in common])
def arm_rho(Z):
    return (rho(Z[:, 0], Z[:, 1]), rho(Z[:, 1], Z[:, 2]), rho(Z[:, 0], Z[:, 3]))
p0, t0, d0 = arm_rho(ZK)
print('共同宇宙 ρ: panel=%.4f tissue=%.4f dual=%.4f' % (p0, t0, d0))
print('Δρ: dual-panel = %+.4f   dual-tissue = %+.4f' % (d0 - p0, d0 - t0))

# ---- 配对基因簇 bootstrap（B = 5,000, seed 20260915）----
B = 5000
dp, dt = [], []
for _ in range(B):
    gs = rng.choice(genes, size=len(genes), replace=True)
    idx = np.concatenate([np.where(GK == g)[0] for g in gs])
    Zb = ZK[idx]
    p, t, d = arm_rho(Zb)
    dp.append(d - p); dt.append(d - t)
dp, dt = np.array(dp), np.array(dt)
res = dict(delta_dual_panel=dict(point=round(float(d0 - p0), 4),
                                 ci95=[round(float(np.percentile(dp, 2.5)), 3), round(float(np.percentile(dp, 97.5)), 3)],
                                 p=round(float(2 * min((dp > 0).mean(), (dp < 0).mean())), 3)),
           delta_dual_tissue=dict(point=round(float(d0 - t0), 4),
                                  ci95=[round(float(np.percentile(dt, 2.5)), 3), round(float(np.percentile(dt, 97.5)), 3)],
                                  p=round(float(2 * min((dt > 0).mean(), (dt < 0).mean())), 3)))
print('bootstrap 结果:', json.dumps(res, ensure_ascii=False))
print('稿件声称: dual-panel = -0.020, 95% CI -0.125 to +0.089, P = 0.74 ; dual-tissue = +0.033, 95% CI -0.217 to +0.274, P = 0.78')

# ---- 出图 ----
FULL = [('panel-only\n(GTEx Whole_Blood\nvs eQTLGen)', 0.4904, 0.4904, 0.21, 0.55, C_PAN),
        ('tissue-only\n(GTEx Whole_Blood\nvs Nerve_Tibial)', 0.4138, 0.4138, 0.27, 0.54, C_TIS),
        ('dual-mismatch\n(GTEx multi-tissue\nvs eQTLGen)', 0.4421, 0.4421, 0.32, 0.55, C_DUA)]
# 用 Fisher-z 精确重算 CI（full-pair 宇宙 n = 159 / 138 / 198）
nn = [159, 138, 198]
CIF = []
for r, n in zip([0.4904, 0.4138, 0.4421], nn):
    z = np.arctanh(r); se = 1 / np.sqrt(n - 3)
    CIF.append((np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)))
print('Fisher-z CI:', [(round(a, 3), round(b, 3)) for a, b in CIF])

fig, (axa, axb) = plt.subplots(1, 2, figsize=(6.65, 3.0),
                               gridspec_kw={'width_ratios': [1.0, 1.25]})
xs = np.arange(3)
for i, ((lab, r, _, _, _, c), (lo, hi), n) in enumerate(zip(FULL, CIF, nn)):
    axa.errorbar(i, r, yerr=[[r - lo], [hi - r]], fmt='o', ms=6.5, color=c, ecolor=c, capsize=4, lw=1.2)
    axa.text(i, hi + 0.015, '$\\rho$ = %.3f\nn = %d' % (r, n), ha='center', fontsize=6.0)
axa.set_xticks(xs); axa.set_xticklabels([p[0] for p in FULL], fontsize=6.0)
axa.set_ylabel('Spearman $\\rho$ (Fisher-z 95% CI)')
axa.set_ylim(0.15, 0.68)
axa.set_title('(a) Decomposition arms, full-pair universes', fontsize=7.4, pad=5)

rows = [('dual \u2212 panel-only', res['delta_dual_panel'], C_PAN),
        ('dual \u2212 tissue-only', res['delta_dual_tissue'], C_TIS)]
for i, (lab, d, c) in enumerate(rows):
    lo, hi = d['ci95']
    axb.errorbar(d['point'], i, xerr=[[d['point'] - lo], [hi - d['point']]], fmt='s', ms=6.5,
                 color=c, ecolor=c, capsize=4, lw=1.2)
    axb.text(hi + 0.008, i, 'P = %.2f' % d['p'], va='center', fontsize=6.2)
axb.axvline(0, color='#333333', lw=0.8, ls='--')
axb.set_yticks([0, 1]); axb.set_yticklabels([r[0] for r in rows], fontsize=6.6)
axb.set_xlabel('$\\Delta\\rho$ (paired gene-cluster bootstrap, B = 5,000)')
axb.set_xlim(-0.30, 0.30); axb.set_ylim(-0.8, 1.8)
axb.set_title('(b) Between-axis difference, 45-gene common universe', fontsize=7.4, pad=5)
fig.tight_layout()
for ext in ['png', 'pdf']:
    src = os.path.join(OUT, f'Fig4.{ext}')
    if os.path.exists(src) and not os.path.exists(os.path.join(BK, f'Fig4.{ext}')):
        shutil.copy2(src, os.path.join(BK, f'Fig4.{ext}'))
    fig.savefig(src, dpi=600)
plt.close(fig)
# 2026-09-20：matplotlib 默认写 RGBA，而图集其余图为 RGB（评审 m-14），统一转 RGB。
from PIL import Image as _Image
_p = os.path.join(OUT, 'Fig4.png')
_im = _Image.open(_p)
if _im.mode != 'RGB':
    _im.convert('RGB').save(_p, dpi=(600, 600))
print('Fig4 redrawn')
json.dump(res, open(os.path.join(P.RES, 'fig4_bootstrap_officialZ.json'), 'w',
                    encoding='utf-8'), indent=1, ensure_ascii=False)
print('backups:', sorted(os.listdir(BK)))
