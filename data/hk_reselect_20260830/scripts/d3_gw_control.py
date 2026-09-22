#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
D3：基因组范围、非架构筛选的随机对照基因集
--------------------------------------------------------------------------
设计
  POOL_A = mashr_Whole_Blood.db 中所有 n.snps>=1 的基因
           - 104-panel 基因
           - 与 panel 共享 >=3 字符前导字母前缀的家族（含 MRPS*/MRPL*/MT-*）
           - curated T2DM/并发症/代谢黑名单（与 hk v2 选择完全同一份）
           → 只要求 Whole_Blood 有模型；**不**要求双组织模型（即撤除架构筛选）
  从 POOL_A 随机抽 600 个基因（seed 20260911）跑同一 S-PrediXcan 算法
  另外单独抽 30 个作为 "GW control" 点估计

分析
  (a) 30 基因零分布：从 600 中重抽 10,000 次，给出富集率分布；
      看 hk v2 (54.8%) / candidate (58.0%) 落在什么分位
  (b) 同一随机样本内按模型可得性分层：WB-only vs WB+NT → 直接检验"双组织架构筛选"本身
      对富集率的影响
  (c) 固定阈值与 BH 两种口径都报
"""
import csv, gzip, io, math, os, random, re, sqlite3, sys, time, json
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

MODEL_DIR = r'E:\workbuddy\BMC Genomics投稿资料\DR，DN，DM芬兰原始数据\mashr_eqtl\eqtl\mashr'
GWAS_DIR = r'E:\workbuddy\BMC Genomics投稿资料\DR，DN，DM芬兰原始数据'
PLINK_PREFIX = r'E:\workbuddy\2026-06-24-05-57-20\tools\ldref\g1000_eur'
REPO = r'E:\workbuddy\eqtl-source-discordance-audit'
OUT = r'E:\workbuddy\2026-09-11-19-30-45\recompute'
os.makedirs(OUT, exist_ok=True)

GWAS_MAP = {'DR': 'finngen_R13_DM_RETINOPATHY_EXMORE.gz',
            'DN': 'finngen_R13_DM_NEPHROPATHY.gz',
            'DPN': 'finngen_R13_DM_NEUROPATHY.gz'}
TISSUES = ['Nerve_Tibial', 'Whole_Blood']
PHENOS = ['DR', 'DN', 'DPN']
SEED = 20260911
N_SAMPLE = 600      # 计算样本量（用于零分布）
N_CONTROL = 30      # 点估计对照集大小

_DECODE = np.array([0.0, 1.0, 2.0, np.nan], dtype=np.float64)


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- pools
EXTRA_FAMILY = re.compile(r'^(MRPS|MRPL|MT-|MTRNR|MTND|MTATP|MTCO|MTCYB)')
DISEASE = {
    'ADCY5','ADRA2A','ANK1','AP3S2','ARAP1','BCAR1','BCL11A','CAMK1D','CCND2','CDKAL1','CDKN2A',
    'CDKN2B','CENTD2','CMIP','DGKB','DUSP8','FTO','GCC1','GCK','GCKR','GIPR','GLIS3','GLP1R','GPSM1',
    'GRB14','HHEX','HMGA1','HMGA2','HNF1A','HNF1B','HNF4A','IDE','IGF1','IGF2BP2','INS','INSR','IRS1',
    'IRS2','JAZF1','KCNJ11','KCNQ1','KLF14','LEPR','MAEA','MC4R','MNX1','MTNR1B','NOTCH2','PAM','PDX1',
    'PEPD','PIK3R1','PPARG','PPARGC1A','PRC1','PROX1','PSMD6','RREB1','SLC16A11','SLC2A2','SLC2A4',
    'SLC30A8','ST6GAL1','TCF7L2','THADA','TP53INP1','TSPAN8','UBE2E2','WFS1','ZBED3','ZFAND6','ADIPOQ',
    'AKT1','AKT2','FOXA2','G6PC2','HK1','MLXIPL','NRXN3','SREBF1','TCF7','PPP1R3B','TMEM154','SSR1',
    'FITM2','RNF6','ANKH','C2CD4A','C2CD4B','VPS13C','CILP2','HNF4G','RASGRP1','C5orf67','ZMIZ1',
    'VEGFA','EPO','AKR1B1','NOS3','ACE','AGT','TGFB1','SERPINE1','MTHFR','APOE','ELMO1','ENPP1',
    'UNC13B','CPVL','CHN2','GREM1','FRMD3','CARS','SP3','ITGA2','ITGB3','ADAM10','ICAM1','SELE',
    'TNF','IL6','CRP','AGER','RAGE','CTGF','CCN2','MMP2','MMP9','TIMP1','HIF1A','PLGF','PGF',
    'LEP','GCGR','PCSK9','LDLR','HMGCR','SREBF2','FASN','ACACA','CPT1A','PPARA','LIPC','CETP',
}


def leading_alpha(s):
    m = re.match(r'^([A-Za-z]+)', s)
    return m.group(1).upper() if m else ''


def build_pool():
    panel = {r['Gene'].upper() for r in csv.DictReader(open(os.path.join(REPO, 'data', 'processed', 'covariate_matrix.csv'), encoding='utf-8'))}
    fams = {leading_alpha(g) for g in panel if len(leading_alpha(g)) >= 3}
    meta = {}
    for t in TISSUES:
        conn = sqlite3.connect(os.path.join(MODEL_DIR, 'mashr_%s.db' % t))
        meta[t] = {gn.upper(): n for _, gn, n in conn.execute('SELECT gene, genename, "n.snps.in.model" FROM extra') if gn}
        conn.close()
    wb = {g: n for g, n in meta['Whole_Blood'].items() if n and n >= 1}
    flow = [('WB model genes', len(wb))]
    cur = {g for g in wb if g not in panel}
    flow.append(('minus 104-panel', len(cur)))
    cur = {g for g in cur if not (any(g.startswith(p) for p in fams) or EXTRA_FAMILY.match(g))}
    flow.append(('minus panel families', len(cur)))
    cur = {g for g in cur if g not in DISEASE}
    flow.append(('minus disease blacklist = POOL_A', len(cur)))
    both = {g for g in cur if g in meta['Nerve_Tibial'] and meta['Nerve_Tibial'][g] and meta['Nerve_Tibial'][g] >= 1}
    flow.append(('of which have BOTH-tissue models', len(both)))
    return sorted(cur), both, meta, flow


# ---------------------------------------------------------------- S-PrediXcan
class Bed:
    def __init__(self, prefix):
        with open(prefix + '.fam') as f:
            self.n = sum(1 for _ in f)
        self.bpp = math.ceil(self.n / 4)
        self.bed_path = prefix + '.bed'
        self.idx = {}
        self.n_snps = 0

    def build_index(self, wanted):
        wanted = set(wanted)
        idx = {}
        with open(self.bed_path.replace('.bed', '.bim')) as f:
            for i, line in enumerate(f):
                rsid = line.split('\t', 2)[1]
                if rsid in wanted and rsid not in idx:
                    idx[rsid] = i
                self.n_snps += 1
        self.idx = idx
        return idx

    def read(self, rsid):
        i = self.idx.get(rsid, -1)
        if i < 0:
            return None
        with open(self.bed_path, 'rb') as f:
            f.seek(3 + i * self.bpp)
            raw = f.read(self.bpp)
        codes = np.frombuffer(raw, dtype=np.uint8)
        shifts = np.arange(4) * 2
        c = ((codes[:, None] >> shifts[None, :]) & 3).ravel()[: self.n]
        d = _DECODE[c]
        m = np.nanmean(d)
        if np.isnan(m):
            m = 0.0
        return np.where(np.isnan(d), m, d)


def load_gwas(pheno, wanted):
    path = os.path.join(GWAS_DIR, GWAS_MAP[pheno])
    wanted = set(wanted)
    out = {}
    t0 = time.time()
    reader = pd.read_csv(path, sep='\t', usecols=['rsids', 'beta', 'sebeta'],
                         dtype={'rsids': 'string', 'beta': 'float64', 'sebeta': 'float64'},
                         chunksize=2_000_000, low_memory=False)
    for chunk in reader:
        sub = chunk[chunk['rsids'].isin(wanted)]
        for rsid, beta, se in zip(sub['rsids'], sub['beta'], sub['sebeta']):
            if se and se > 0 and rsid not in out:
                out[rsid] = beta / se
    log('    GWAS %s: hits %d (%.1fs)' % (pheno, len(out), time.time() - t0))
    return out


def main():
    t_start = time.time()
    pool, both, meta, flow = build_pool()
    log('=== POOL 构建流程 ===')
    for n, c in flow:
        log('   %-38s %6d' % (n, c))

    rnd = random.Random(SEED)
    sample = rnd.sample(pool, N_SAMPLE)
    ctrl = rnd.sample(pool, N_CONTROL)
    log('\nseed = %d ; 计算样本 %d 基因 ; GW control %d 基因' % (SEED, len(sample), len(ctrl)))
    log('  GW control: %s' % ', '.join(sorted(ctrl)))
    log('  样本中双组织可测比例: %d/%d = %.1f%%' % (sum(1 for g in sample if g in both), len(sample),
                                                 100 * sum(1 for g in sample if g in both) / len(sample)))

    genes = sorted(set(sample) | set(ctrl))
    bed = Bed(PLINK_PREFIX)
    conns = {t: sqlite3.connect(os.path.join(MODEL_DIR, 'mashr_%s.db' % t)) for t in TISSUES}
    gid = {t: {gn.upper(): gid for gid, gn in conns[t].execute('SELECT gene, genename FROM extra') if gn} for t in TISSUES}

    model_snps = {}
    for t in TISSUES:
        for g in genes:
            if g in gid[t] and g in meta[t]:
                rows = conns[t].execute('SELECT rsid, weight FROM weights WHERE gene=?', (gid[t][g],)).fetchall()
                if rows:
                    model_snps[(t, g)] = (meta[t][g], rows)
    wanted = {rs for _, w in model_snps.values() for rs, _ in w}
    log('\n模型 SNP（去重）= %d' % len(wanted))
    t0 = time.time(); bed.build_index(wanted)
    log('PLINK: %d 个体 / %d SNPs，索引到 %d（%.1fs）' % (bed.n, bed.n_snps, len(bed.idx), time.time() - t0))

    results = []
    for pheno in PHENOS:
        gw = load_gwas(pheno, wanted)
        for (t, g), (n_model, mod) in model_snps.items():
            matched = [(rs, w) for rs, w in mod if rs in gw]
            if not matched:
                continue
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
                dm = np.vstack(dl).T
                cv = np.cov(dm, rowvar=False)
                pv = float(w_a @ cv @ w_a)
                if pv <= 0: continue
                z = float(np.dot(w_a, z_a) / math.sqrt(pv))
            p = math.erfc(abs(z) / math.sqrt(2))
            results.append({'gene': g, 'tissue': t, 'trait': pheno, 'zscore': round(z, 4),
                            'pvalue': p, 'n_snps_model': n_model, 'n_snps_matched': len(matched)})
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(OUT, 'd3_gw_sample_raw.csv'), index=False)
    log('\nTWAS 结果 %d 行 -> d3_gw_sample_raw.csv  (%.1f min)' % (len(df), (time.time() - t_start) / 60))

    # ---------------- ACAT-O (stable) ----------------
    piv = df.pivot_table(index=['gene', 'trait'], columns='tissue', values=['zscore', 'pvalue'])
    recs = []
    for (g, tr), r in piv.iterrows():
        ps = [r[('pvalue', t)] for t in TISSUES if ('pvalue', t) in r and not pd.isna(r[('pvalue', t)])]
        if not ps:
            continue
        if len(ps) == 1:
            pacat = ps[0]
        else:
            pp = np.clip(np.asarray(ps, float), 1e-300, 1 - 1e-15)
            T = float(np.sum(np.where(pp < 1e-8, 1.0 / (pp * np.pi), np.tan((0.5 - pp) * np.pi))))
            pacat = 1.0 if T <= 0 else float(math.atan(1.0 / T) / math.pi)
        has_nt = ('pvalue', 'Nerve_Tibial') in r.index and not pd.isna(r[('pvalue', 'Nerve_Tibial')])
        recs.append({'gene': g, 'trait': tr, 'p_acat': pacat, 'n_tissues': len(ps),
                     'in_both': bool(has_nt), 'z_wb': r.get(('zscore', 'Whole_Blood'), np.nan)})
    A = pd.DataFrame(recs)
    A.to_csv(os.path.join(OUT, 'd3_gw_acat.csv'), index=False)
    log('ACAT-O 行数 = %d ；其中双组织 %d' % (len(A), int(A.in_both.sum())))

    # ---------------- analysis ----------------
    L = []
    def say(*a): L.append(' '.join(str(x) for x in a))
    def bh(p):
        p = np.asarray(p, float); out = np.full(len(p), np.nan); ok = ~np.isnan(p)
        if ok.sum() == 0: return out
        v = np.clip(p[ok], 1e-300, 1.0); o = np.argsort(v); m = len(v)
        q = v[o] * m / np.arange(1, m + 1); q = np.minimum.accumulate(q[::-1])[::-1]
        res = np.empty(m); res[o] = np.minimum(q, 1.0); out[ok] = res
        return out

    say('=' * 88)
    say('D3. 基因组范围、非架构筛选随机对照（GTEx v8，cot 形式 ACAT-O）')
    say('=' * 88)
    say('SEED=%d ; POOL_A=%d 基因（只要求 Whole_Blood 有 MASHR 模型）' % (SEED, len(pool)))
    say('计算样本 = %d 基因；样本中双组织可测 %d (%.1f%%)'
        % (len(sample), sum(1 for g in sample if g in both), 100 * sum(1 for g in sample if g in both) / len(sample)))

    def rate_of(genes_sub, thr=None, use_bh=True):
        s = A[A.gene.isin(genes_sub)]
        if thr is None:
            k = 0; n = 0
            for tr, idx in s.groupby('trait').groups.items():
                q = bh(s.loc[idx, 'p_acat'].values)
                k += int(np.nansum(q < 0.05)); n += len(idx)
            return k, n
        v = s['p_acat'].dropna().values
        return int((v < thr).sum()), len(v)

    say('\n-- (a) GW control 30 基因点估计 --')
    for lab, thr in [('BH q<0.05 (per trait)', None), ('p<0.05', 0.05), ('p<0.01', 0.01), ('p<0.00385', 0.00385)]:
        k, n = rate_of(ctrl, thr)
        say('   %-22s %3d/%3d = %5.1f%%' % (lab, k, n, 100 * k / n))

    say('\n-- (b) 30 基因零分布（从 %d 基因样本中重抽, B=10000）--' % len(sample))
    rng = np.random.default_rng(SEED)
    samp_arr = np.array(sample)
    boot_bh = []; boot_f05 = []; boot_f00385 = []
    for _ in range(10000):
        gs = set(rng.choice(samp_arr, N_CONTROL, replace=False))
        k, n = rate_of(gs, None); boot_bh.append(100 * k / n)
        k, n = rate_of(gs, 0.05); boot_f05.append(100 * k / n)
        k, n = rate_of(gs, 0.00385); boot_f00385.append(100 * k / n)
    for lab, boot, obs in [('BH q<0.05', boot_bh, 54.8), ('p<0.05', boot_f05, None), ('p<0.00385', boot_f00385, None)]:
        b = np.array(boot)
        say('   %-12s 零分布: 中位 %5.1f%%  IQR %5.1f-%5.1f%%  2.5-97.5%% %5.1f-%5.1f%%'
            % (lab, np.median(b), np.percentile(b, 25), np.percentile(b, 75), np.percentile(b, 2.5), np.percentile(b, 97.5)))
        for name, val in [('hk v2 control 54.8%', 54.8), ('candidate 58.0%', 58.0), ('non-candidate 56.8%', 56.8), ('T2DM control 59.6%', 59.6)]:
            pct = 100 * (b < val).mean()
            say('        %-22s -> 位于零分布第 %5.1f 百分位' % (name, pct))

    say('\n-- (c) 同一随机样本内：双组织可测 vs 仅 WB --')
    for lab, sub in [('both-tissue testable', A[A.in_both]), ('WB-only', A[~A.in_both])]:
        if len(sub) == 0: continue
        k, n = rate_of(set(sub.gene))
        say('   %-22s genes=%3d pairs=%3d  BH 富集 %3d/%3d = %5.1f%%  |Z|中位=%.2f'
            % (lab, sub.gene.nunique(), len(sub), k, n, 100 * k / n,
               np.nanmedian(np.abs(sub.z_wb.dropna().values)) if sub.z_wb.notna().any() else float('nan')))
    say('\n   注：BOTH 组的 p_acat 由两组织 ACAT-O 得到，WB-only 组即单组织 p，两者阈值口径一致。')

    open(os.path.join(OUT, 'd3_results.txt'), 'w', encoding='utf-8').write('\n'.join(L))
    json.dump({'seed': SEED, 'pool_size': len(pool), 'sample_size': len(sample),
               'control': sorted(ctrl), 'flow': flow,
               'null_bh': [float(np.percentile(np.array(boot_bh), p)) for p in (2.5, 50, 97.5)],
               'null_p05': [float(np.percentile(np.array(boot_f05), p)) for p in (2.5, 50, 97.5)]},
              open(os.path.join(OUT, 'd3_summary.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    log('\nDONE  (%.1f min)' % ((time.time() - t_start) / 60))


if __name__ == '__main__':
    main()
