# -*- coding: utf-8 -*-
"""P0 收口 Step 2：从官方 Z 重绘 Figure 5 与 Figure 7"""
import os, csv, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import paths_config as P  # 统一路径入口（2026-09-20）

NEW = P.need(P.DATA_Z, '官方 MetaXcan Z 数据层')
OUT = P.OUT_MAIN
BK = os.path.join(OUT, '_backup_before_officialZ_redraw_20260917')
os.makedirs(BK, exist_ok=True)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.linewidth': 0.7,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
C_GTEX, C_EQTL, C_POOL = '#C0392B', '#2471A3', '#7D3C98'

def rd(name, key):
    with open(os.path.join(NEW, name), encoding='utf-8-sig') as f:
        return {tuple(r[k] for k in key): r for r in csv.DictReader(f)}

GT = rd('gtex_official_Z.csv', ['Gene', 'Trait'])
EQ = rd('eqtlgen_official_Z.csv', ['Gene', 'Trait'])
GRP = rd('gene_groups_TableS1_official.csv', ['Gene'])
f = lambda x: (float(x) if x not in (None, '', 'NA') else np.nan)

def save(fig, name):
    for ext in ['png', 'pdf']:
        src = os.path.join(OUT, f'{name}.{ext}')
        if os.path.exists(src) and not os.path.exists(os.path.join(BK, f'{name}.{ext}')):
            shutil.copy2(src, os.path.join(BK, f'{name}.{ext}'))
        fig.savefig(src, dpi=600)
    # 2026-09-20：matplotlib 默认写 RGBA，而图集其余图为 RGB（评审 m-14），统一转 RGB。
    from PIL import Image as _Image
    _p = os.path.join(OUT, name + '.png')
    _im = _Image.open(_p)
    if _im.mode != 'RGB':
        _im.convert('RGB').save(_p, dpi=(600, 600))
    plt.close(fig)

# ================= Figure 5 =================
cands = [k[0] for k in GRP if 'andidate' in GRP[k]['Group'] and 'Non' not in GRP[k]['Group']]
rows = []
for g in cands:
    e = EQ.get((g, 'DR')); t = GT.get((g, 'DR'))
    if not e:
        continue
    ze = f(e['Z_eQTLGen'])
    if np.isnan(ze):
        continue
    zt = f(t['Z_Nerve_Tibial']) if t else np.nan
    rows.append(dict(gene=g, ze=ze, zt=zt,
                     snps=e['Model_SNPs'],
                     q=float(e['BH_q']) if e['BH_q'] not in ('', 'NA') else np.nan))
rows.sort(key=lambda r: -abs(r['ze']))
top = rows[:11]
print('Figure 5 — 官方 Z 下 DR 候选基因 |Ze| 前 11:')
for r in top:
    print('   %-9s eQTLGen %+7.3f  GTEx NT %s  model %s' %
          (r['gene'], r['ze'], ('%+.3f' % r['zt']) if not np.isnan(r['zt']) else 'no model', r['snps']))

fig, ax = plt.subplots(figsize=(5.4, 3.6))
y = np.arange(len(top))[::-1]
h = 0.36
for i, r in enumerate(top):
    yy = y[i]
    if not np.isnan(r['zt']):
        ax.barh(yy + h / 2, r['zt'], h, color=C_GTEX, edgecolor='black', lw=0.5, label='GTEx v8 Nerve_Tibial' if i == 0 else None)
    ax.barh(yy - h / 2, r['ze'], h, color=C_EQTL, edgecolor='black', lw=0.5, label='eQTLGen whole blood' if i == 0 else None)
ax.axvline(0, color='#333333', lw=0.7)
ax.set_yticks(y); ax.set_yticklabels([r['gene'] for r in top], fontsize=7)
ax.set_xlabel('S-PrediXcan Z-score (DR)')
ax.legend(frameon=False, fontsize=6.6, loc='lower right')
ymax = max(max(abs(r['ze']) for r in top), max((abs(r['zt']) for r in top if not np.isnan(r['zt'])), default=0))
ax.set_xlim(-ymax * 1.15, ymax * 1.15)
for i, r in enumerate(top):
    ax.text(r['ze'] + (0.25 if r['ze'] > 0 else -0.25), y[i] - h / 2,
            '%+.2f' % r['ze'], va='center', ha='left' if r['ze'] > 0 else 'right', fontsize=6.2)
    if not np.isnan(r['zt']):
        ax.text(r['zt'] + (0.25 if r['zt'] > 0 else -0.25), y[i] + h / 2,
                '%+.2f' % r['zt'], va='center', ha='left' if r['zt'] > 0 else 'right', fontsize=6.2)
fig.tight_layout(); save(fig, 'Fig5')
print('Fig5 redrawn (official Z)')

# ================= Figure 7 =================
PH = ['DR', 'DN', 'DPN']
gtex_z = [f(GT[('RNH1', t)]['Z_Nerve_Tibial']) for t in PH]
eqtl_z = [f(EQ[('RNH1', t)]['Z_eQTLGen']) for t in PH]
print('\nFigure 7 panel (a) — 官方 Z: GTEx NT', np.round(gtex_z, 3), ' eQTLGen', np.round(eqtl_z, 3))
zF, zU = 2.31, 0.72
# panel (b) 数值：直接取归档发布值，与 crosscohort_TableS4_official.csv 第 1 行（= Table S4a）逐字一致。
# 2026-09-23 修正：原为 m = 1.515、se = 0.795 —— 那是把归档的 2 位小数（1.51 / 0.79）**再补一位**得来的
#   中间值，而 '%.2f' % 0.795 恰好进位成 0.80，于是图内写 SE 0.80 而归档与正文都是 SE 0.79
#   （正文：pooled Z = +1.51, SE 0.79, P = 0.056；1.51/0.79 → P = 0.056 ✓，1.51/0.80 → P = 0.059 ✗）。
#   改用归档原值后不再有二次舍入。
m = 1.51; se = 0.79; pi = (-0.33, 3.36)
fig, (axa, axb) = plt.subplots(1, 2, figsize=(6.65, 2.9))
x = np.arange(3); w = 0.36
axa.bar(x - w / 2, gtex_z, w, color=C_GTEX, edgecolor='black', lw=0.5, label='GTEx v8 Nerve_Tibial')
axa.bar(x + w / 2, eqtl_z, w, color=C_EQTL, edgecolor='black', lw=0.5, label='eQTLGen whole blood')
for i in range(3):
    axa.text(x[i] - w / 2, gtex_z[i] + (0.08 if gtex_z[i] >= 0 else -0.30), '%+.2f' % gtex_z[i], ha='center', fontsize=6.2)
    axa.text(x[i] + w / 2, eqtl_z[i] + (0.08 if eqtl_z[i] >= 0 else -0.30), '%+.2f' % eqtl_z[i], ha='center', fontsize=6.2)
axa.axhline(0, color='#333333', lw=0.7)
axa.set_xticks(x); axa.set_xticklabels(PH)
axa.set_ylabel('RNH1 S-PrediXcan Z')
axa.set_ylim(0, max(max(gtex_z), max(eqtl_z)) * 1.28)
axa.legend(frameon=False, fontsize=6.2, loc='upper right')
axa.set_title('(a) RNH1 by phenotype (official Z)', fontsize=7.6, pad=4)

labels = ['FinnGen R13\n(discovery)', 'UK Biobank\nGCST90043640', 'Random-effects\npooled']
ys = [2, 1, 0]
axb.errorbar(zF, ys[0], xerr=0, fmt='o', ms=7, color=C_GTEX, zorder=4)
axb.errorbar(zU, ys[1], xerr=0, fmt='o', ms=7, color=C_EQTL, zorder=4)
axb.errorbar(m, ys[2], xerr=se, fmt='s', ms=7, color=C_POOL, ecolor=C_POOL, capsize=3, zorder=4)
axb.axvspan(pi[0], pi[1], color=C_POOL, alpha=0.10, zorder=1)
axb.axvline(0, color='#333333', lw=0.7, ls='--')
axb.set_yticks(ys); axb.set_yticklabels(labels, fontsize=6.6)
axb.set_xlabel('RNH1 TWAS Z-score (DR, eQTLGen weights)')
axb.text(m, ys[2] + 0.30, 'Z = %+.2f (SE %.2f)' % (m, se), ha='center', fontsize=6.2, color=C_POOL)
axb.text(pi[1], -0.05, '95%% PI %+.2f to %+.2f' % pi, ha='right', fontsize=6.2, color=C_POOL)
axb.set_ylim(-0.7, 2.6)
axb.set_title('(b) Cross-cohort replication', fontsize=7.6, pad=4)
fig.tight_layout(); save(fig, 'Fig7')
print('Fig7 redrawn (official Z)')
print('\nbackups ->', BK)
print('backup contents:', sorted(os.listdir(BK)))
