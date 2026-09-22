#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
S-PrediXcan 复现管道（严格复刻 analyses/run_hk_control.py 的算法）

复刻要点（保证与已发表的 104-panel 结果同口径）：
  * 权重：sqlite 读取 mashr_{tissue}.db 的 weights 表（rsid, weight）
  * GWAS：FinnGen R13，Z = beta / sebeta，仅保留 rsIDs（与原始脚本一致）
  * LD：1000G EUR PLINK bed（SNP-major），就地计算 w' Σ w
  * 统计量：Z_twas = w'z / sqrt(w' Σ w)；单 SNP 时 Σ = var(g)
  * 缺失基因型：以该 SNP 的非缺失均值填补（与原始脚本一致）

相对原始脚本的优化（不改变任何数值结果）：
  * bim 索引由 O(n) 线性扫描改为一次性建字典（原脚本每次查询扫全表，不可扩展）
  * GWAS 读取改为 pandas 分块 + 按需过滤，避免 800MB gz 全量驻留内存
  * bed 解码向量化

用法：
  python 02_spredixcan.py --genes hk_genes_v2.txt --out hk_twas_v2.csv --validate
"""
import argparse, gzip, io, math, os, sqlite3, sys, time
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

MODEL_DIR = r'E:\workbuddy\BMC Genomics投稿资料\DR，DN，DM芬兰原始数据\mashr_eqtl\eqtl\mashr'
GWAS_DIR = r'E:\workbuddy\BMC Genomics投稿资料\DR，DN，DM芬兰原始数据'
PLINK_PREFIX = r'E:\workbuddy\2026-06-24-05-57-20\tools\ldref\g1000_eur'

GWAS_MAP = {'DR': 'finngen_R13_DM_RETINOPATHY_EXMORE.gz',
            'DN': 'finngen_R13_DM_NEPHROPATHY.gz',
            'DPN': 'finngen_R13_DM_NEUROPATHY.gz'}

TISSUES = ['Nerve_Tibial', 'Whole_Blood']
PHENOS = ['DR', 'DN', 'DPN']

# 2-bit -> dosage 的查找表，与 run_hk_control.py 完全一致
_DECODE = np.array([0.0, 1.0, 2.0, np.nan], dtype=np.float64)


# ----------------------------------------------------------------------------
def load_models(tissues):
    """{tissue: {GENE: (gene_id, rows[(rsid, weight), ...], n_model)}}"""
    out = {}
    for t in tissues:
        conn = sqlite3.connect(os.path.join(MODEL_DIR, f'mashr_{t}.db'))
        meta = {}
        for gene, genename, n in conn.execute('SELECT gene, genename, "n.snps.in.model" FROM extra'):
            if genename:
                meta[genename.upper()] = (gene, n)
        out[t] = {'conn': conn, 'meta': meta}
    return out


def gene_weights(conn, gene_id):
    rows = conn.execute(
        'SELECT rsid, weight, ref_allele, eff_allele FROM weights WHERE gene=?', (gene_id,)
    ).fetchall()
    return [(r[0], r[1]) for r in rows]


class Bed:
    """PLINK 1 bed (SNP-major) 读取器，仅索引需要的 rsID"""

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
        # 每个个体占 2 bit，LSB 在前
        shifts = np.arange(4) * 2
        c = ((codes[:, None] >> shifts[None, :]) & 3).ravel()[: self.n]
        d = _DECODE[c]
        m = np.nanmean(d)
        if np.isnan(m):
            m = 0.0
        return np.where(np.isnan(d), m, d)


# ----------------------------------------------------------------------------
def load_gwas(pheno, wanted_rsids):
    """返回 {rsid: Z}，仅保留模型需要的 rsID"""
    path = os.path.join(GWAS_DIR, GWAS_MAP[pheno])
    wanted = set(wanted_rsids)
    out = {}
    t0 = time.time()
    reader = pd.read_csv(path, sep='\t', usecols=['rsids', 'beta', 'sebeta'],
                         dtype={'rsids': 'string', 'beta': 'float64', 'sebeta': 'float64'},
                         chunksize=2_000_000, low_memory=False)
    n = 0
    for chunk in reader:
        sub = chunk[chunk['rsids'].isin(wanted)]
        for rsid, beta, se in zip(sub['rsids'], sub['beta'], sub['sebeta']):
            if se and se > 0 and rsid not in out:
                out[rsid] = beta / se
        n += len(chunk)
    print(f'    GWAS {pheno}: {n:,} 行扫描，命中 {len(out)} 个模型 SNP '
          f'({time.time()-t0:.1f}s)', flush=True)
    return out


# ----------------------------------------------------------------------------
def run(genes, bed, models, phenos=PHENOS, tissues=TISSUES):
    genes = [g.upper() for g in genes]

    # 1) 收集全部模型 SNP
    model_snps = {}
    for t in tissues:
        m = models[t]
        for g in genes:
            if g in m['meta']:
                gene_id, n_model = m['meta'][g]
                model_snps[(t, g)] = (gene_id, n_model, gene_weights(m['conn'], gene_id))
    wanted = {rsid for _, _, w in model_snps.values() for rsid, _ in w}
    print(f'  模型 SNP 总数（去重）: {len(wanted)}', flush=True)

    # 2) 建 bed 索引
    t0 = time.time()
    bed.build_index(wanted)
    print(f'  PLINK: {bed.n} 个体 / {bed.n_snps:,} SNPs，索引到 {len(bed.idx)} 个 '
          f'({time.time()-t0:.1f}s)', flush=True)

    results = []
    for pheno in phenos:
        gw = load_gwas(pheno, wanted)
        for t in tissues:
            m = models[t]
            for g in genes:
                if (t, g) not in model_snps:
                    continue
                gene_id, n_model, mod = model_snps[(t, g)]
                matched = [(rs, w) for rs, w in mod if rs in gw]
                n_matched = len(matched)
                if n_matched == 0:
                    continue
                w_a = np.array([w for _, w in matched], dtype=np.float64)
                z_a = np.array([gw[rs] for rs, _ in matched], dtype=np.float64)

                if n_matched == 1:
                    d = bed.read(matched[0][0])
                    if d is None:
                        continue
                    ve = np.var(d, ddof=1)
                    if ve <= 0:
                        ve = 2 * 0.05 * 0.95
                    pv = w_a[0] ** 2 * ve
                    z = w_a[0] * z_a[0] / math.sqrt(pv) if pv > 0 else 0.0
                else:
                    dl = [bed.read(rs) for rs, _ in matched]
                    if any(x is None for x in dl):
                        continue
                    dm = np.vstack(dl).T
                    cv = np.cov(dm, rowvar=False)
                    pv = float(w_a @ cv @ w_a)
                    if pv <= 0:
                        continue
                    z = float(np.dot(w_a, z_a) / math.sqrt(pv))

                p = 2 * 0.5 * math.erfc(abs(z) / math.sqrt(2))
                results.append({'gene': g, 'tissue': t, 'trait': pheno,
                                'zscore': round(z, 4), 'pvalue': p,
                                'n_snps_model': n_model, 'n_snps_matched': n_matched})
    return pd.DataFrame(results)


# ----------------------------------------------------------------------------
def validate(df_new, ref_dir):
    """与已发布的 gtex_{tissue}_{trait}.csv 逐基因比对"""
    print('\n' + '=' * 72)
    print('验证：与已发布 104-panel GTEx 结果比对')
    print('=' * 72)
    allok = True
    for t in TISSUES:
        for ph in PHENOS:
            ref_path = os.path.join(ref_dir, f'gtex_{t}_{ph}.csv')
            if not os.path.exists(ref_path):
                continue
            ref = pd.read_csv(ref_path)
            new = df_new[(df_new.tissue == t) & (df_new.trait == ph)]
            mrg = ref.merge(new, on='gene', suffixes=('_ref', '_new'))
            if len(mrg) == 0:
                continue
            diff = (mrg['zscore_ref'] - mrg['zscore_new']).abs()
            n_bad = int((diff > 1e-6).sum())
            ok = n_bad == 0
            allok &= ok
            print(f'  {t:<13} {ph:<4} 比对基因 {len(mrg):>3}  '
                  f'最大 |ΔZ| = {diff.max():.3e}  '
                  f'不一致 {n_bad:>2}  {"OK" if ok else "FAIL"}')
            if not ok:
                bad = mrg[diff > 1e-6][['gene', 'zscore_ref', 'zscore_new',
                                        'n_snps_matched_ref', 'n_snps_matched_new']]
                print(bad.head(15).to_string(index=False))
    print('=' * 72)
    print('总体：', '全部一致 PASS' if allok else '存在不一致 FAIL')
    return allok


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genes', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--validate', action='store_true')
    ap.add_argument('--ref-dir', default=r'E:\workbuddy\eqtl-source-discordance-audit\data\processed')
    args = ap.parse_args()

    with open(args.genes, encoding='utf-8') as f:
        genes = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    genes = [g.strip().upper() for g in genes]
    print(f'输入基因: {len(genes)}')

    t0 = time.time()
    bed = Bed(PLINK_PREFIX)
    models = load_models(TISSUES)
    df = run(genes, bed, models)
    for t in TISSUES:
        models[t]['conn'].close()

    df.to_csv(args.out, index=False)
    print(f'\n写出 {len(df)} 行 -> {args.out}   ({time.time()-t0:.1f}s)')

    if args.validate:
        validate(df, args.ref_dir)


if __name__ == '__main__':
    main()
