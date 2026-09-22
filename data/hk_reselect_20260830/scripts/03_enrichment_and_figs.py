#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
富集率计算 + Table S6 + Figure 9a/9b（看家基因真阴性对照 v2）

要点
 * 多组织整合：ACAT-O（与手稿 Methods 一致，p 截断于 1e-300）
 * FDR：Benjamini-Hochberg，分层 = group × phenotype × eQTL weight source
 * 同时重算 104 panel 三个原始分组，保证 HK 与对照基线同管道、同口径
 * 输出与已发表 Table 2 的逐格对照，暴露不可复现部分
"""
import io, json, math, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

ROOT = r'E:\workbuddy\hk_reselect_20260830'
PROC = r'E:\workbuddy\eqtl-source-discordance-audit\data\processed'
TRAITS = ['DR', 'DN', 'DPN']
TISSUES = ['Nerve_Tibial', 'Whole_Blood']

GROUP_LABEL = {'30 HOTAIR Candidate': 'Candidate',
               '44 Non-Candidate': 'Non-Candidate',
               '30 T2DM Control': 'T2DM Control',
               'HK_v2': 'Housekeeping v2'}


# ---------------------------------------------------------------- ACAT-O
def acat_o(p, trunc=1e-300):
    """ACAT-O: T = Σ tan((0.5-p)π); p = 1/2 - arctan(T)/π

    数值稳定性：
      * p 极小时 tan(π/2 - pπ) 病态，改用恒等近似 tan(π/2-ε) = cot(ε) ≈ 1/ε
      * 末段用恒等式 1/2 - arctan(T)/π = arctan(1/T)/π（T>0 恒成立），
        避免 T 极大时 arctan(T) 饱和到 π/2 而把 p 压成 0
        （已发表 gtex_acat_o_results.csv 正是因此把 p 截断在 2.11e-15）
    """
    p = np.asarray(p, dtype=float)
    p = np.clip(p, trunc, 1.0 - 1e-15)
    terms = np.where(p < 1e-8, 1.0 / (p * np.pi), np.tan((0.5 - p) * np.pi))
    T = float(np.sum(terms))
    if T <= 0:
        return 1.0
    return float(math.atan(1.0 / T) / math.pi)


def bh_q(p):
    """Benjamini-Hochberg q 值"""
    p = np.asarray(p, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    q_sorted = p[order] * m / np.arange(1, m + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(q_sorted, 1.0)
    return q


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    rx, ry = pd.Series(x).rank().values, pd.Series(y).rank().values
    mx, my = rx.mean(), ry.mean()
    num = ((rx - mx) * (ry - my)).sum()
    den = math.sqrt(((rx - mx) ** 2).sum() * ((ry - my) ** 2).sum())
    rho = num / den if den > 0 else float('nan')
    # Fisher z 近似 CI
    if n > 3 and abs(rho) < 1:
        z = 0.5 * math.log((1 + rho) / (1 - rho))
        se = 1 / math.sqrt(n - 3)
        lo, hi = math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se)
    else:
        lo = hi = float('nan')
    t = rho * math.sqrt((n - 2) / (1 - rho ** 2)) if abs(rho) < 1 else float('nan')
    p = 2 * 0.5 * math.erfc(abs(t) / math.sqrt(2)) if not math.isnan(t) else float('nan')
    return rho, p, lo, hi, n


# ---------------------------------------------------------------- 载入
def load_all():
    cov = pd.read_csv(os.path.join(PROC, 'covariate_matrix.csv'))
    cov['G'] = cov['Gene'].str.upper()
    grp = dict(zip(cov['G'], cov['Group']))
    src = dict(zip(cov['G'], cov['Source']))
    unused = dict(zip(cov['G'], cov['PullDown_Unused']))

    panel = pd.read_csv(os.path.join(ROOT, 'validate_panel104.csv'))
    hk = pd.read_csv(os.path.join(ROOT, 'hk_twas_v2_raw.csv'))
    panel['grp'] = panel['gene'].map(grp)
    hk['grp'] = 'HK_v2'
    hk['gene'] = hk['gene'].str.upper()
    df = pd.concat([panel, hk], ignore_index=True)
    return df, src, unused


def build_arms(df):
    """每 gene×trait 生成 NT / WB / ACAT-O 三个口径"""
    rows = []
    for (g, t), sub in df.groupby(['gene', 'trait']):
        rec = {'gene': g, 'trait': t, 'grp': sub['grp'].iloc[0],
               'z_Nerve_Tibial': np.nan, 'z_Whole_Blood': np.nan,
               'n_snps_matched_NT': np.nan, 'n_snps_matched_WB': np.nan}
        for _, r in sub.iterrows():
            rec['z_' + r['tissue']] = r['zscore']
            rec['n_snps_matched_NT' if r['tissue'] == 'Nerve_Tibial' else 'n_snps_matched_WB'] = r['n_snps_matched']
        # ACAT-O
        ps = []
        for ti in TISSUES:
            z = rec['z_' + ti]
            if not np.isnan(z):
                ps.append(math.erfc(abs(z) / math.sqrt(2)))
        if ps:
            p_acat = acat_o(ps) if len(ps) > 1 else ps[0]
            rec['p_ACAT_O'] = p_acat
            rec['n_tissues'] = len(ps)
            # 等效 Z（用于绘图与方向一致性）
            if p_acat > 0:
                rec['z_ACAT_O'] = math.copysign(
                    abs(_z_from_p(p_acat)), rec['z_Nerve_Tibial'] if not np.isnan(rec['z_Nerve_Tibial'])
                    else rec['z_Whole_Blood'])
            else:
                rec['z_ACAT_O'] = math.copysign(40.0, rec['z_Nerve_Tibial'] if not np.isnan(rec['z_Nerve_Tibial'])
                                                else rec['z_Whole_Blood'])
        else:
            rec['p_ACAT_O'] = np.nan
            rec['z_ACAT_O'] = np.nan
            rec['n_tissues'] = 0
        rows.append(rec)
    return pd.DataFrame(rows)


def _z_from_p(p):
    """双尾 p -> |z|（Wichura 近似，够用）"""
    if p <= 0:
        return 40.0
    # 二分
    lo, hi = 0.0, 60.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if math.erfc(mid / math.sqrt(2)) > p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ---------------------------------------------------------------- 富集
def enrichment(arms, pcol):
    out = []
    for (g, t), sub in arms.groupby(['grp', 'trait']):
        s = sub.dropna(subset=[pcol])
        n = len(s)
        if n == 0:
            continue
        q = bh_q(s[pcol].values)
        n_fdr = int((q < 0.05).sum())
        out.append({'Group': GROUP_LABEL.get(g, g), 'Trait': t,
                    'N_tested': n, 'N_FDR': n_fdr, 'FDR_pct': round(100 * n_fdr / n, 1)})
    return pd.DataFrame(out).sort_values(['Group', 'Trait'])


def main():
    df, src, unused = load_all()
    arms = build_arms(df)
    arms.to_csv(os.path.join(ROOT, 'arms_all_groups.csv'), index=False)

    print('=' * 78)
    print('一、各组可测基因数（GTEx v8 MASHR，本次重跑）')
    print('=' * 78)
    piv = arms.pivot_table(index='grp', columns='trait', values='p_ACAT_O', aggfunc='count')
    print(piv.to_string())

    print('\n' + '=' * 78)
    print('二、FDR 富集率（BH q<0.05，分层 = group × trait × source）')
    print('=' * 78)
    tables = {}
    for label, pcol in [('GTEx ACAT-O (multi-tissue)', 'p_ACAT_O'),
                        ('GTEx Nerve_Tibial', 'z_Nerve_Tibial'),
                        ('GTEx Whole_Blood', 'z_Whole_Blood')]:
        if pcol.startswith('z'):
            tmp = arms.copy()
            tmp['_p'] = tmp[pcol].apply(lambda z: math.erfc(abs(z) / math.sqrt(2)) if not np.isnan(z) else np.nan)
            tables[label] = enrichment(tmp, '_p')
        else:
            tables[label] = enrichment(arms, pcol)
        print(f'\n-- {label} --')
        print(tables[label].to_string(index=False))

    # ---------- 与已发表 Table 2 对照 ----------
    print('\n' + '=' * 78)
    print('三、与已发表 Table 2（GTEx ACAT-O 列）逐格对照')
    print('=' * 78)
    pub = pd.read_csv(os.path.join(PROC, 'enrichment_comparison.csv'))
    pub['Group'] = pub['Group_gtex'].replace({
        '30_HOTAIR_Candidate': 'Candidate', '44_NonCandidate_HOTAIR': 'Non-Candidate',
        '30_T2DM_Control': 'T2DM Control'})
    new = tables['GTEx ACAT-O (multi-tissue)']
    print(f'{"Group":<14}{"Trait":<6}{"已发表 N":>9}{"重算 N":>8}{"已发表 FDR%":>13}{"重算 FDR%":>12}   一致?')
    n_same = 0
    for _, r in pub.iterrows():
        m = new[(new.Group == r['Group']) & (new.Trait == r['Trait'])]
        if len(m) == 0:
            continue
        m = m.iloc[0]
        ok = (int(m.N_FDR) == int(r['N_FDR005_gtex'])) and (int(m.N_tested) == int(r['N_Tested_gtex']))
        n_same += ok
        print(f'{r["Group"]:<14}{r["Trait"]:<6}{int(r["N_Tested_gtex"]):>9}{int(m.N_tested):>8}'
              f'{r["GTEx_FDR%"]:>13.1f}{m.FDR_pct:>12.1f}   {"OK" if ok else "x"}')
    print(f'\n复现 {n_same}/9 格')

    # ---------- 组间 Fisher 精确检验 ----------
    print('\n' + '=' * 78)
    print('四、看家基因 vs 各原始分组：Fisher 精确检验（GTEx ACAT-O，三表型合并）')
    print('=' * 78)
    sig = arms.dropna(subset=['p_ACAT_O']).copy()
    sig['q'] = np.nan
    for (g, t), idx in sig.groupby(['grp', 'trait']).groups.items():
        sig.loc[idx, 'q'] = bh_q(sig.loc[idx, 'p_ACAT_O'].values)
    pooled = sig.groupby('grp').apply(
        lambda d: pd.Series({'N': len(d), 'FDR': int((d.q < 0.05).sum()),
                             'pct': round(100 * (d.q < 0.05).mean(), 1)}), include_groups=False)
    print(pooled.to_string())
    print()
    from itertools import product
    try:
        from scipy.stats import fisher_exact
    except ImportError:
        fisher_exact = None
    hk_row = pooled.loc['HK_v2']
    for g in ['30 HOTAIR Candidate', '44 Non-Candidate', '30 T2DM Control']:
        if g not in pooled.index:
            continue
        r = pooled.loc[g]
        tab = [[int(hk_row.FDR), int(hk_row.N - hk_row.FDR)],
               [int(r.FDR), int(r.N - r.FDR)]]
        if fisher_exact:
            orr, p = fisher_exact(tab)
            print(f'  HK_v2 ({int(hk_row.FDR)}/{int(hk_row.N)} = {hk_row.pct}%)  vs  '
                  f'{GROUP_LABEL[g]:<14} ({int(r.FDR)}/{int(r.N)} = {r.pct}%)  '
                  f'OR = {orr:.2f}  Fisher P = {p:.3f}')
        else:
            print(f'  HK_v2 vs {GROUP_LABEL[g]}: {tab}')

    # ---------- Table S6 ----------
    hk_sel = [l.strip() for l in open(os.path.join(ROOT, 'hk_genes_v2.txt'), encoding='utf-8')
              if l.strip() and not l.startswith('#')]
    meta = json.load(open(os.path.join(ROOT, 'hk_selection_meta.json'), encoding='utf-8'))
    t6 = arms[arms.grp == 'HK_v2'].pivot_table(index='gene', columns='trait',
                                               values=['z_Nerve_Tibial', 'z_Whole_Blood', 'p_ACAT_O'])
    t6 = t6.reindex(hk_sel).reset_index()
    t6.to_csv(os.path.join(ROOT, 'TableS6_hk_control_v2.csv'), index=False)
    print(f'\nTable S6（30 个看家基因的逐基因 TWAS 结果）-> TableS6_hk_control_v2.csv')

    # ---------- Figure 9a / 9b ----------
    h = arms[arms.grp == 'HK_v2']
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0))
    specs = [('z_ACAT_O', 'z_Nerve_Tibial', 'GTEx MASHR multi-tissue Z', 'GTEx Nerve_Tibial Z', 'a'),
             ('z_Whole_Blood', 'z_Nerve_Tibial', 'GTEx Whole_Blood Z', 'GTEx Nerve_Tibial Z', 'b')]
    stats = {}
    for ax, (xcol, ycol, xl, yl, lab) in zip(axes, specs):
        s = h.dropna(subset=[xcol, ycol])
        rho, p, lo, hi, n = spearman(s[xcol].values, s[ycol].values)
        stats[lab] = dict(rho=rho, p=p, ci=[lo, hi], n=n)
        ax.scatter(s[xcol], s[ycol], s=42, facecolors='none', edgecolors='#1f77b4', linewidths=1.4)
        lim = max(np.abs(s[xcol]).max(), np.abs(s[ycol]).max()) * 1.15
        ax.plot([-lim, lim], [-lim, lim], ls='--', lw=1, color='#888888')
        ax.axhline(0, lw=0.6, color='#cccccc')
        ax.axvline(0, lw=0.6, color='#cccccc')
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel(xl, fontsize=10)
        ax.set_ylabel(yl, fontsize=10)
        ax.tick_params(labelsize=9)
        txt = f'Spearman rho = {rho:+.2f} (n = {n})'
        ax.text(0.03, 0.96, txt, transform=ax.transAxes, ha='left', va='top', fontsize=10)
        ax.text(0.5, -0.14, f'({lab}) Housekeeping true-negative control (v2)',
                transform=ax.transAxes, ha='center', va='top', fontsize=10, fontweight='bold')
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.20)
    plt.savefig(os.path.join(ROOT, 'Figure9_hk_v2.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(ROOT, 'Figure9_hk_v2.pdf'), bbox_inches='tight')
    print('\nFigure 9 统计：')
    for k, v in stats.items():
        print(f'  ({k}) rho = {v["rho"]:+.3f}  95%CI [{v["ci"][0]:+.2f}, {v["ci"][1]:+.2f}]  '
              f'n = {v["n"]} pairs  P = {v["p"]:.3g}')
    print('-> Figure9_hk_v2.png / .pdf')

    # ---------- 汇总 JSON ----------
    summary = {
        'hk_genes_v2': hk_sel,
        'selection': {k: meta[k] for k in ['seed', 'n_universe_hrt', 'n_eligible_pool', 'selection_rules']},
        'enrichment_acat_o': tables['GTEx ACAT-O (multi-tissue)'].to_dict('records'),
        'enrichment_nerve_tibial': tables['GTEx Nerve_Tibial'].to_dict('records'),
        'enrichment_whole_blood': tables['GTEx Whole_Blood'].to_dict('records'),
        'figure9': stats,
        'reproduced_published_table2_cells': f'{n_same}/9',
    }
    with open(os.path.join(ROOT, 'hk_v2_summary.json'), 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print('\n汇总 -> hk_v2_summary.json')


if __name__ == '__main__':
    main()
