import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
def _find_root(start):
    p = start
    for _ in range(8):
        if _os.path.isdir(_os.path.join(p, "data", "processed")):
            return p
        nxt = _os.path.dirname(p)
        if nxt == p:
            break
        p = nxt
    return _os.path.dirname(start)
ROOT = _find_root(HERE)
# -*- coding: utf-8 -*-
"""S1: 重建 primary arm 并做基因层 cluster bootstrap。"""
import sys, numpy as np, pandas as pd
from scipy.stats import spearmanr, binomtest

D = r"E:\workbuddy\eqtl-source-discordance-audit\data\processed"
import os
os.chdir(D)

cov = pd.read_csv("covariate_matrix.csv")
g = cov["Group"].str.lower()
nc = set(cov[g.str.contains("non.candidate", regex=True)]["Gene"])
t2 = set(cov[g.str.contains("t2dm")]["Gene"])
cand = set(cov[(g.str.contains("candidate")) & (~g.str.contains("non"))]["Gene"])
print(f"[groups] NC={len(nc)} T2DM={len(t2)} CAND={len(cand)}", flush=True)

def load_gtex(tis):
    return pd.concat([pd.read_csv(f"gtex_{tis}_{t}.csv").assign(Trait=t)
                      for t in ["DR", "DN", "DPN"]], ignore_index=True)

WB = load_gtex("Whole_Blood")
NT = load_gtex("Nerve_Tibial")
eq = pd.read_csv("eqtlgen_spredixcan_harmonized_results.csv")
eq_ok = eq[eq["grp"] != "Excluded_NonTestbed"]
print(f"[src] WB={len(WB)}/{WB.gene.nunique()}g NT={len(NT)}/{NT.gene.nunique()}g "
      f"eQTLGen(testbed)={len(eq_ok)}/{eq_ok.gene.nunique()}g", flush=True)

def build(pool, gtex_df):
    W = gtex_df[gtex_df.gene.isin(pool)][["gene", "Trait", "zscore"]].rename(
        columns={"gene": "Gene", "zscore": "Z_GTEx"})
    E = eq_ok[eq_ok.gene.isin(pool)][["gene", "trait", "zscore"]].rename(
        columns={"gene": "Gene", "trait": "Trait", "zscore": "Z_eQTLGen"})
    m = W.merge(E, on=["Gene", "Trait"])
    cnt = m.groupby("Gene").size()
    keep = set(cnt[cnt == 3].index)
    m3 = m[m.Gene.isin(keep)].copy()
    m3["Same"] = np.sign(m3.Z_GTEx) == np.sign(m3.Z_eQTLGen)
    return m3

pool = nc | t2
for lbl, gd in (("GTEx Whole_Blood", WB), ("GTEx Nerve_Tibial", NT)):
    m3 = build(pool, gd)
    n = len(m3); k = int(m3.Same.sum())
    if n:
        rho, p = spearmanr(m3.Z_GTEx, m3.Z_eQTLGen)
        print(f"\n[{lbl} vs eQTLGen] genes={m3.Gene.nunique()} pairs={n} consistent={k} "
              f"({100*k/n:.1f}%) rho={rho:.3f} P={p:.4f} binomP={binomtest(k,n,0.5).pvalue:.4f}", flush=True)
        for t in ["DR", "DN", "DPN"]:
            s = m3[m3.Trait == t]
            print(f"    {t}: {int(s.Same.sum())}/{len(s)} ({100*s.Same.mean():.1f}%)", flush=True)
        m3.to_csv(_os.path.join(HERE, f"RECON_primary_{lbl.split()[-1]}.csv"), index=False)

print("\n[target] manuscript: 32 genes / 96 pairs / 61 (63.5%) / rho=0.32 P=0.001 "
      "/ per-trait 21,21,19 of 32", flush=True)
