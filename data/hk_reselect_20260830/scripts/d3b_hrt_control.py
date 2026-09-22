#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第二重对照：HRT Atlas 全集 × 仅 Whole_Blood 模型（撤除双组织架构筛选），n=818，全部计算。"""
import csv, io, math, os, random, re, sqlite3, sys, time, json
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
MODEL_DIR = r'E:\workbuddy\BMC Genomics投稿资料\DR，DN，DM芬兰原始数据\mashr_eqtl\eqtl\mashr'
GWAS_DIR = r'E:\workbuddy\BMC Genomics投稿资料\DR，DN，DM芬兰原始数据'
PLINK_PREFIX = r'E:\workbuddy\2026-06-24-05-57-20\tools\ldref\g1000_eur'
REPO = r'E:\workbuddy\eqtl-source-discordance-audit'
OUT = r'E:\workbuddy\2026-09-11-19-30-45\recompute'
HRT_RAW = os.path.join(OUT, 'Human_Mouse_Common_raw.csv')
SEED = 20260911
N_CONTROL = 30
GWAS_MAP = {'DR': 'finngen_R13_DM_RETINOPATHY_EXMORE.gz',
            'DN': 'finngen_R13_DM_NEPHROPATHY.gz',
            'DPN': 'finngen_R13_DM_NEUROPATHY.gz'}
TISSUES = ['Nerve_Tibial', 'Whole_Blood']
_DECODE = np.array([0.0, 1.0, 2.0, np.nan], dtype=np.float64)
def log(*a): print(*a, flush=True)
t0 = time.time()

# ---- pool（与 hrt_verify.py 完全一致） ----
hrt = set()
for line in open(HRT_RAW, encoding='utf-8', errors='replace'):
    line = line.strip()
    if not line or line.lower().startswith('mouse'): continue
    p = line.split(';')
    if len(p) >= 2 and p[1].strip(): hrt.add(p[1].strip().upper())
panel = {r['Gene'].upper() for r in csv.DictReader(open(os.path.join(REPO,'data','processed','covariate_matrix.csv'), encoding='utf-8'))}
EXTRA_FAMILY = re.compile(r'^(MRPS|MRPL|MT-|MTRNR|MTND|MTATP|MTCO|MTCYB)')
def lead(s):
    m = re.match(r'^([A-Za-z]+)', s); return m.group(1).upper() if m else ''
fams = {lead(g) for g in panel if len(lead(g)) >= 3}
DISEASE = set("""ADCY5 ADRA2A ANK1 AP3S2 ARAP1 BCAR1 BCL11A CAMK1D CCND2 CDKAL1 CDKN2A CDKN2B CENTD2 CMIP DGKB DUSP8 FTO GCC1 GCK GCKR GIPR GLIS3 GLP1R GPSM1 GRB14 HHEX HMGA1 HMGA2 HNF1A HNF1B HNF4A IDE IGF1 IGF2BP2 INS INSR IRS1 IRS2 JAZF1 KCNJ11 KCNQ1 KLF14 LEPR MAEA MC4R MNX1 MTNR1B NOTCH2 PAM PDX1 PEPD PIK3R1 PPARG PPARGC1A PRC1 PROX1 PSMD6 RREB1 SLC16A11 SLC2A2 SLC2A4 SLC30A8 ST6GAL1 TCF7L2 THADA TP53INP1 TSPAN8 UBE2E2 WFS1 ZBED3 ZFAND6 ADIPOQ AKT1 AKT2 FOXA2 G6PC2 HK1 MLXIPL NRXN3 SREBF1 TCF7 PPP1R3B TMEM154 SSR1 FITM2 RNF6 ANKH C2CD4A C2CD4B VPS13C CILP2 HNF4G RASGRP1 C5orf67 ZMIZ1 VEGFA EPO AKR1B1 NOS3 ACE AGT TGFB1 SERPINE1 MTHFR APOE ELMO1 ENPP1 UNC13B CPVL CHN2 GREM1 FRMD3 CARS SP3 ITGA2 ITGB3 ADAM10 ICAM1 SELE TNF IL6 CRP AGER RAGE CTGF CCN2 MMP2 MMP9 TIMP1 HIF1A PLGF PGF LEP GCGR PCSK9 LDLR HMGCR SREBF2 FASN ACACA CPT1A PPARA LIPC CETP""".split())
cur = hrt - panel
cur = {g for g in cur if not (any(g.startswith(p) for p in fams) or EXTRA_FAMILY.match(g))}
cur = {g for g in cur if g not in DISEASE}
meta = {}
for t in TISSUES:
    conn = sqlite3.connect(os.path.join(MODEL_DIR, 'mashr_%s.db' % t))
    meta[t] = {gn.upper(): n for _, gn, n in conn.execute('SELECT gene, genename, "n.snps.in.model" FROM extra') if gn}
    conn.close()
pool = {g for g in cur if g in meta['Whole_Blood'] and meta['Whole_Blood'][g] and meta['Whole_Blood'][g] >= 1}
both = {g for g in pool if g in meta['Nerve_Tibial'] and meta['Nerve_Tibial'][g] and meta['Nerve_Tibial'][g] >= 1}
log('HRT x WB-only pool = %d (both-tissue %d, WB-only %d)' % (len(pool), len(both), len(pool) - len(both)))

rnd = random.Random(SEED)
ctrl = rnd.sample(sorted(pool), N_CONTROL)
log('GW-HRT control (seed %d): %s' % (SEED, ', '.join(sorted(ctrl))))

# ---- S-PrediXcan（与已验证 02_spredixcan.py 同算法） ----
class Bed:
    def __init__(self, prefix):
        with open(prefix + '.fam') as f: self.n = sum(1 for _ in f)
        self.bpp = math.ceil(self.n / 4); self.bed_path = prefix + '.bed'; self.idx = {}
    def build_index(self, wanted):
        wanted = set(wanted); idx = {}; self.n_snps = 0
        with open(self.bed_path.replace('.bed', '.bim')) as f:
            for i, line in enumerate(f):
                rsid = line.split('\t', 2)[1]
                if rsid in wanted and rsid not in idx: idx[rsid] = i
                self.n_snps += 1
        self.idx = idx; return idx
    def read(self, rsid):
        i = self.idx.get(rsid, -1)
        if i < 0: return None
        with open(self.bed_path, 'rb') as f:
            f.seek(3 + i * self.bpp); raw = f.read(self.bpp)
        codes = np.frombuffer(raw, dtype=np.uint8)
        shifts = np.arange(4) * 2
        c = ((codes[:, None] >> shifts[None, :]) & 3).ravel()[: self.n]
        d = _DECODE[c]; m = np.nanmean(d)
        if np.isnan(m): m = 0.0
        return np.where(np.isnan(d), m, d)

def load_gwas(pheno, wanted):
    wanted = set(wanted); out = {}
    reader = pd.read_csv(os.path.join(GWAS_DIR, GWAS_MAP[pheno]), sep='\t',
                         usecols=['rsids', 'beta', 'sebeta'],
                         dtype={'rsids': 'string', 'beta': 'float64', 'sebeta': 'float64'},
                         chunksize=2_000_000, low_memory=False)
    for chunk in reader:
        sub = chunk[chunk['rsids'].isin(wanted)]
        for rsid, beta, se in zip(sub['rsids'], sub['beta'], sub['sebeta']):
            if se and se > 0 and rsid not in out: out[rsid] = beta / se
    log('  GWAS %s hits %d' % (pheno, len(out)))
    return out

genes = sorted(pool)
conns = {t: sqlite3.connect(os.path.join(MODEL_DIR, 'mashr_%s.db' % t)) for t in TISSUES}
gid = {t: {gn.upper(): gid for gid, gn in conns[t].execute('SELECT gene, genename FROM extra') if gn} for t in TISSUES}
model_snps = {}
for t in TISSUES:
    for g in genes:
        if g in gid[t] and g in meta[t]:
            rows = conns[t].execute('SELECT rsid, weight FROM weights WHERE gene=?', (gid[t][g],)).fetchall()
            if rows: model_snps[(t, g)] = (meta[t][g], rows)
wanted = {rs for _, w in model_snps.values() for rs, _ in w}
bed = Bed(PLINK_PREFIX)
log('model SNPs dedup = %d ; indexing...' % len(wanted))
bed.build_index(wanted)
log('indexed %d of %d SNPs (%.1f min)' % (len(bed.idx), bed.n_snps, (time.time()-t0)/60))

results = []
for pheno in ['DR', 'DN', 'DPN']:
    gw = load_gwas(pheno, wanted)
    for (t, g), (n_model, mod) in model_snps.items():
        matched = [(rs, w) for rs, w in mod if rs in gw]
        if not matched: continue
        w_a = np.array([w for _, w in matched]); z_a = np.array([gw[rs] for rs, _ in matched])
        if len(matched) == 1:
            d = bed.read(matched[0][0])
            if d is None: continue
            ve = np.var(d, ddof=1)
            if ve <= 0: ve = 2 * 0.05 * 0.95
            pv = w_a[0] ** 2 * ve
            z = w_a[0] * z_a[0] / math.sqrt(pv) if pv > 0 else 0.0
        else:
            dl = [bed.read(rs) for rs, _ in matched]
            if any(x is None for x in dl): continue
            dm = np.vstack(dl).T; cv = np.cov(dm, rowvar=False)
            pv = float(w_a @ cv @ w_a)
            if pv <= 0: continue
            z = float(np.dot(w_a, z_a) / math.sqrt(pv))
        results.append({'gene': g, 'tissue': t, 'trait': pheno, 'zscore': round(z, 4),
                        'pvalue': math.erfc(abs(z)/math.sqrt(2)), 'n_snps_model': n_model,
                        'n_snps_matched': len(matched)})
df = pd.DataFrame(results)
df.to_csv(os.path.join(OUT, 'd3b_hrt_pool_raw.csv'), index=False)
log('TWAS rows=%d (%.1f min)' % (len(df), (time.time()-t0)/60))

# ---- ACAT-O + 分析 ----
def bh(p):
    p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = ~np.isnan(p)
    if ok.sum() == 0: return out
    v = np.clip(p[ok], 1e-300, 1.0); o = np.argsort(v); m = len(v)
    q = v[o]*m/np.arange(1, m+1); q = np.minimum.accumulate(q[::-1])[::-1]
    res = np.empty(m); res[o] = np.minimum(q, 1.0); out[ok] = res
    return out
piv = df.pivot_table(index=['gene','trait'], columns='tissue', values='pvalue')
recs = []
for (g, tr), r in piv.iterrows():
    ps = [r[t] for t in TISSUES if t in r.index and not pd.isna(r[t])]
    if not ps: continue
    if len(ps) == 1: pacat = ps[0]
    else:
        pp = np.clip(np.asarray(ps, float), 1e-300, 1-1e-15)
        T = float(np.sum(np.where(pp < 1e-8, 1.0/(pp*np.pi), np.tan((0.5-pp)*np.pi))))
        pacat = 1.0 if T <= 0 else float(math.atan(1.0/T)/math.pi)
    has_nt = 'Nerve_Tibial' in r.index and not pd.isna(r['Nerve_Tibial'])
    recs.append({'gene': g, 'trait': tr, 'p_acat': pacat, 'n_tissues': len(ps), 'in_both': bool(has_nt)})
A = pd.DataFrame(recs)
A.to_csv(os.path.join(OUT, 'd3b_hrt_acat.csv'), index=False)

L = []
def say(*a): L.append(' '.join(str(x) for x in a))
say('='*88)
say('第二重对照：HRT Atlas 全集 × 仅 Whole_Blood 模型（撤除双组织架构筛选）')
say('='*88)
say('池 = %d 基因（HRT 1129 → R2/R3/R4 后 × WB 模型）；双组织可测 %d，仅 WB %d'
    % (len(pool), len(both), len(pool)-len(both)))
def rate(gs, thr=None):
    s = A[A.gene.isin(gs)]
    if thr is None:
        k = n = 0
        for tr, idx in s.groupby('trait').groups.items():
            q = bh(s.loc[idx, 'p_acat'].values)
            k += int(np.nansum(q < 0.05)); n += len(idx)
        return k, n
    v = s['p_acat'].dropna().values
    return int((v < thr).sum()), len(v)
say('\n-- 30 基因对照（seed %d）--' % SEED)
for lab, thr in [('BH q<0.05', None), ('p<0.05', 0.05), ('p<0.01', 0.01), ('p<0.00385', 0.00385)]:
    k, n = rate(ctrl, thr)
    say('   %-12s %3d/%3d = %5.1f%%' % (lab, k, n, 100*k/n))
say('\n-- 30 基因零分布（从全部 %d 基因重抽, B=10000, seed 同）--' % len(pool))
rng = np.random.default_rng(SEED)
arr = np.array(sorted(pool))
boot = {'bh': [], 'p05': [], 'p00385': []}
for _ in range(10000):
    gs = set(rng.choice(arr, N_CONTROL, replace=False))
    k, n = rate(gs, None); boot['bh'].append(100*k/n)
    k, n = rate(gs, 0.05); boot['p05'].append(100*k/n)
    k, n = rate(gs, 0.00385); boot['p00385'].append(100*k/n)
for lab, key in [('BH q<0.05','bh'), ('p<0.05','p05'), ('p<0.00385','p00385')]:
    b = np.array(boot[key])
    say('   %-12s 中位 %5.1f%%  IQR %5.1f-%5.1f%%  2.5-97.5%% %5.1f-%5.1f%%'
        % (lab, np.median(b), np.percentile(b,25), np.percentile(b,75), np.percentile(b,2.5), np.percentile(b,97.5)))
    for nm, val in [('hk v2 control 54.8%', 54.8), ('candidate 58.0%', 58.0),
                    ('non-candidate 56.8%', 56.8), ('T2DM control 59.6%', 59.6),
                    ('genome-wide ctrl 69.1%', 69.1)]:
        say('        %-24s -> 第 %5.1f 百分位' % (nm, 100*(b < val).mean()))
say('\n-- 池内分层：双组织可测 vs 仅 WB --')
for lab, sub in [('both-tissue', A[A.in_both]), ('WB-only', A[~A.in_both])]:
    if len(sub) == 0: continue
    k, n = rate(set(sub.gene))
    say('   %-12s genes=%3d pairs=%3d  BH %3d/%3d = %5.1f%%' % (lab, sub.gene.nunique(), len(sub), k, n, 100*k/n))
open(os.path.join(OUT, 'd3b_results.txt'), 'w', encoding='utf-8').write('\n'.join(L))
json.dump({'seed': SEED, 'pool_n': len(pool), 'control': sorted(ctrl),
           'null_bh': [float(np.percentile(np.array(boot['bh']), p)) for p in (2.5, 50, 97.5)]},
          open(os.path.join(OUT, 'd3b_summary.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
log('\nDONE (%.1f min)' % ((time.time()-t0)/60))
