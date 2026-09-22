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
grpof = {}
for _, r in cov.iterrows():
    gl = r["Group"].lower()
    grpof[r["Gene"]] = "T2DM" if "t2dm" in gl else ("NC" if "non-candidate" in gl else "CAND")

WB = pd.concat([pd.read_csv(f"gtex_Whole_Blood_{t}.csv").assign(Trait=t) for t in ["DR","DN","DPN"]], ignore_index=True)
eq = pd.read_csv("eqtlgen_spredixcan_harmonized_results.csv")
eq_ok = eq[eq.grp != "Excluded_NonTestbed"]
pool = nc | t2
W = WB[WB.gene.isin(pool)][["gene","Trait","zscore"]].rename(columns={"gene":"Gene","zscore":"Z_GTEx"})
E = eq_ok[eq_ok.gene.isin(pool)][["gene","trait","zscore"]].rename(columns={"gene":"Gene","trait":"Trait","zscore":"Z_eQTLGen"})
m = W.merge(E, on=["Gene","Trait"])
cnt = m.groupby("Gene").size(); m = m[m.Gene.isin(cnt[cnt==3].index)].copy()

nan_g = sorted(set(m[m.Z_eQTLGen.isna()].Gene))
print("Z_eQTLGen 为 NaN 的基因:", nan_g)
print("其分组:", {x: grpof[x] for x in nan_g})
print("NaN 对数:", m.Z_eQTLGen.isna().sum())
m["Same"] = np.sign(m.Z_GTEx) == np.sign(m.Z_eQTLGen)
print("\n这 4 个基因方向一致情况:")
for x in nan_g:
    s = m[m.Gene == x]
    print(f"  {x:9s} {grpof[x]:5s} pairs={len(s)} nonNaN_Z={s.Z_eQTLGen.notna().sum()}")

# ---- primary arm: 去掉 Z_eQTLGen 为 NaN ----
p = m[m.Z_eQTLGen.notna()].copy()
print(f"\n=== PRIMARY ARM ===  genes={p.Gene.nunique()} pairs={len(p)}")
print("  分组:", p.groupby(p.Gene.map(grpof)).Gene.nunique().to_dict())
n = len(p); k = int(p.Same.sum())
rho, rp = spearmanr(p.Z_GTEx, p.Z_eQTLGen)
print(f"  consistent={k}/{n} ({100*k/n:.2f}%)  rho={rho:.3f} (P={rp:.4f})  binomP={binomtest(k,n,0.5).pvalue:.4f}")
for t in ["DR","DN","DPN"]:
    s = p[p.Trait==t]; print(f"    {t}: {int(s.Same.sum())}/{len(s)} ({100*s.Same.mean():.1f}%)")
print("\n[target] 32 genes / 96 pairs / 61 (63.5%) / rho=0.32 P=0.001 / binomP=0.010")
p.to_csv(_os.path.join(HERE, "PRIMARY_ARM_96.csv"), index=False)
