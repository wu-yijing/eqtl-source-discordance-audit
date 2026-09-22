#!/usr/bin/env python3
"""HK-only S-PrediXcan - inline version, no backtick issues"""
import sqlite3, os, math, gzip, sys, time
import numpy as np
import pandas as pd
from scipy import stats

t0 = time.time()
MODEL_DIR = r'E:\workbuddy\hotair\mashr_eqtl\eqtl\mashr'
PLINK_PREFIX = r'C:\Users\Administrator\.workbuddy\tools\1000g_eur\g1000_eur'
GWAS_DIR = r'E:\workbuddy\hotair'

HK_ONLY = ['B2M','UBC','TBP','HPRT1','GUSB','SDHA','HMBS','YWHAZ','PPIA','IPO8',
           'POLR2A','TFRC','ALDOA','PGK1','LDHA','TPI1','NONO','PUM1','PSMG2',
           'EEF1A1','RPL32','RPS20','CYC1']

for tissue in ['Nerve_Tibial', 'Whole_Blood']:
    for pheno in ['DR', 'DN', 'DPN']:
        print(f'=== {tissue} x {pheno} ===', flush=True)
        
        model_db = os.path.join(MODEL_DIR, f'mashr_{tissue}.db')
        conn = sqlite3.connect(model_db)
        
        extra = {}
        col_n = 'n.snps.in.model'
        q = 'SELECT gene, genename, "' + col_n + '" FROM extra'
        for row in conn.execute(q).fetchall():
            extra[row[1].upper()] = {'id': row[0], 'n_snps': row[2]}
        
        available = {g: extra[g.upper()] for g in HK_ONLY if g.upper() in extra}
        missing = [g for g in HK_ONLY if g.upper() not in extra]
        print(f'  Available: {len(available)}, Missing: {missing}', flush=True)
        
        if len(available) == 0:
            print(f'  No HK genes in model - skipping', flush=True)
            conn.close()
            continue
        
        # PLINK
        class PlinkReader:
            def __init__(self, prefix):
                self.bed = prefix + '.bed'
                self.bim = prefix + '.bim'
                with open(prefix + '.fam') as f: self.n = sum(1 for _ in f)
                self.snps = []
                with open(self.bim) as f:
                    for line in f:
                        p = line.strip().split()
                        self.snps.append({'rsid': p[1]})
                self.n_snps = len(self.snps)
                self.bpp = math.ceil(self.n / 4)
                print(f'  PLINK: {self.n} ind, {self.n_snps} SNPs', flush=True)
            def get_idx(self, rsid):
                for i, s in enumerate(self.snps):
                    if s['rsid'] == rsid: return i
                return -1
            def read(self, idx):
                with open(self.bed, 'rb') as f:
                    f.seek(3 + idx * self.bpp)
                    raw = f.read(self.bpp)
                d = np.zeros(self.n, dtype=np.float64)
                for i in range(self.n):
                    g = (raw[i//4] >> (i%4)*2) & 3
                    d[i] = 0.0 if g==0 else (1.0 if g==1 else (2.0 if g==2 else np.nan))
                m = np.nanmean(d)
                if np.isnan(m): m = 0.0
                return np.nan_to_num(d, nan=m)
        
        plink = PlinkReader(PLINK_PREFIX)
        
        # GWAS
        gw = {}
        gwas_map = {'DR': 'finngen_R13_DM_RETINOPATHY_EXMORE.gz',
                    'DN': 'finngen_R13_DM_NEPHROPATHY.gz',
                    'DPN': 'finngen_R13_DM_NEUROPATHY.gz'}
        gwas_file = os.path.join(GWAS_DIR, gwas_map[pheno])
        print(f'  Loading GWAS...', flush=True)
        with gzip.open(gwas_file, 'rt') as f:
            hdr = f.readline().strip().split('\t')
            cm = {h.lower(): i for i, h in enumerate(hdr)}
            for line in f:
                parts = line.strip().split('\t')
                rsid = parts[cm.get('rsids', 0)]
                if ':' in rsid and not rsid.startswith('rs'):
                    rsid = rsid.split(':')[0]
                try:
                    beta = float(parts[cm['beta']])
                    se = float(parts[cm['sebeta']])
                    if se > 0:
                        gw[rsid] = beta / se
                except:
                    pass
        print(f'  GWAS: {len(gw)} variants', flush=True)
        
        # Compute TWAS
        results = []
        for gname, info in available.items():
            gene = info['id']
            rows = conn.execute('SELECT rsid, weight, ref_allele, eff_allele FROM weights WHERE gene=?', (gene,)).fetchall()
            mod = [{'rsid': r[0], 'weight': r[1]} for r in rows]
            
            w, zs, matched = [], [], []
            for s in mod:
                if s['rsid'] in gw:
                    w.append(s['weight'])
                    zs.append(gw[s['rsid']])
                    matched.append(s)
            
            if len(matched) == 0:
                continue
            
            w_a = np.array(w)
            zs_a = np.array(zs)
            
            if len(matched) == 1:
                idx = plink.get_idx(matched[0]['rsid'])
                if idx < 0: continue
                d = plink.read(idx)
                ve = np.var(d, ddof=1)
                if ve <= 0: ve = 2*0.05*0.95
                pv = w_a[0]**2 * ve
                twas_z = w_a[0] * zs_a[0] / math.sqrt(pv) if pv > 0 else 0
            else:
                dl = []
                for s in matched:
                    idx = plink.get_idx(s['rsid'])
                    if idx < 0: break
                    dl.append(plink.read(idx))
                if len(dl) != len(matched): continue
                dm = np.vstack(dl).T
                cv = np.cov(dm, rowvar=False)
                pv = w_a @ cv @ w_a
                if pv <= 0: continue
                twas_z = np.dot(w_a, zs_a) / math.sqrt(pv)
            
            p = 2 * stats.norm.sf(abs(twas_z))
            results.append({'gene': gname, 'tissue': tissue, 'trait': pheno,
                           'zscore': round(twas_z, 4), 'pvalue': p,
                           'n_snps_matched': len(matched)})
            print(f'  {gname}: Z={twas_z:.2f}', flush=True)
        
        conn.close()
        
        # Save - both possible output paths
        out_paths = [
            rf'data/raw/gtex_raw_output/gtex_{tissue}_{pheno}.csv',
            rf'E:\workbuddy\2026-06-28-17-42-31\gtex_raw_output\gtex_{tissue}_{pheno}.csv'
        ]
        
        df_new = pd.DataFrame(results)
        for out in out_paths:
            if os.path.exists(out):
                df_existing = pd.read_csv(out)
                existing_hk_mask = df_existing['gene'].str.upper().isin([g.upper() for g in HK_ONLY])
                df_remaining = df_existing[~existing_hk_mask]
                df_combined = pd.concat([df_remaining, df_new], ignore_index=True)
                df_combined.to_csv(out, index=False)
                print(f'  Saved to {out}: {len(df_combined)} rows', flush=True)

print(f'\nTotal time: {time.time()-t0:.1f}s', flush=True)
print('ALL DONE!', flush=True)
