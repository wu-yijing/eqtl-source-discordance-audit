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
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr, binomtest
D = r"E:\workbuddy\eqtl-source-discordance-audit\data\processed"
os.chdir(D)
cov = pd.read_csv("covariate_matrix.csv"); g = cov["Group"].str.lower()
nc = set(cov[g.str.contains("non.candidate", regex=True)]["Gene"])
t2 = set(cov[g.str.contains("t2dm")]["Gene"])

WB = pd.concat([pd.read_csv(f"gtex_Whole_Blood_{t}.csv").assign(Trait=t) for t in ["DR","DN","DPN"]], ignore_index=True)
eq = pd.read_csv("eqtlgen_spredixcan_harmonized_results.csv")
eq_ok = eq[eq.grp != "Excluded_NonTestbed"]
pool = nc | t2
W = WB[WB.gene.isin(pool)][["gene","Trait","zscore","n_snps_matched","n_snps_model"]].rename(columns={"gene":"Gene","zscore":"Z_GTEx"})
E = eq_ok[eq_ok.gene.isin(pool)][["gene","trait","zscore"]].rename(columns={"gene":"Gene","trait":"Trait","zscore":"Z_eQTLGen"})
m = W.merge(E, on=["Gene","Trait"])
cnt = m.groupby("Gene").size(); keep = set(cnt[cnt==3].index)
m3 = m[m.Gene.isin(keep)].copy()
m3["grp"] = np.where(m3.Gene.isin(t2), "T2DM", "NC")
print("构成:", m3.groupby('grp').Gene.nunique().to_dict(), "| 总基因", m3.Gene.nunique())
print("NaN 检查  Z_GTEx:", m3.Z_GTEx.isna().sum(), " Z_eQTLGen:", m3.Z_eQTLGen.isna().sum())
print("SNP 匹配数 <5 的对数:", (m3.n_snps_matched < 5).sum())
m3["Same"] = np.sign(m3.Z_GTEx) == np.sign(m3.Z_eQTLGen)
# 逐基因一致数
per = m3.groupby(["Gene","grp"]).agg(k=("Same","sum"), n=("Same","size"),
        minSNP=("n_snps_matched","min"), maxSNP=("n_snps_matched","max")).reset_index()
per["all_same"] = per.k == per.n
per["all_diff"] = per.k == 0
print("\n全一致基因:", per.all_same.sum(), "| 全不一致基因:", per.all_diff.sum(), "| 混合:", len(per)-per.all_same.sum()-per.all_diff.sum())
print("\n=== 12 对全部不一致的基因（候选排除对象）===")
print(per[per.all_diff][["Gene","grp","k","n","minSNP","maxSNP"]].to_string(index=False))
print("\n=== 逐基因一致数分布 ===")
print(per.sort_values(["grp","k"]).to_string(index=False))
m3.to_csv(_os.path.join(HERE, "RECON_primary_WB_36g.csv"), index=False)
per.to_csv(_os.path.join(HERE, "RECON_per_gene_36g.csv"), index=False)
