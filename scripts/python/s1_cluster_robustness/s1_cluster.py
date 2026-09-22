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
"""S1: primary arm 基因层 cluster-robust 重算（numpy 向量化版）。"""
import os, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, binomtest, norm

D = r"E:\workbuddy\eqtl-source-discordance-audit\data\processed"
OUT = r"E:\workbuddy\2026-09-10-18-23-32"
os.chdir(D)
RNG = np.random.default_rng(20260910)
B = 10000
TRAITS = ["DR", "DN", "DPN"]

cov = pd.read_csv("covariate_matrix.csv")
grpof = {}
for _, r in cov.iterrows():
    gl = r["Group"].lower()
    grpof[r["Gene"]] = "T2DM" if "t2dm" in gl else ("NC" if "non-candidate" in gl else "CAND")

WB = pd.concat([pd.read_csv(f"gtex_Whole_Blood_{t}.csv").assign(Trait=t) for t in TRAITS], ignore_index=True)
eq = pd.read_csv("eqtlgen_spredixcan_harmonized_results.csv")
eq_ok = eq[eq.grp != "Excluded_NonTestbed"]

def build(pool):
    W = WB[WB.gene.isin(pool)][["gene","Trait","zscore"]].rename(columns={"gene":"Gene","zscore":"Z_GTEx"})
    E = eq_ok[eq_ok.gene.isin(pool)][["gene","trait","zscore"]].rename(columns={"gene":"Gene","trait":"Trait","zscore":"Z_eQTLGen"})
    m = W.merge(E, on=["Gene","Trait"])
    m = m[m.Z_eQTLGen.notna()].copy()
    c = m.groupby("Gene").size(); m = m[m.Gene.isin(c[c == 3].index)].copy()
    m["Same"] = np.sign(m.Z_GTEx) == np.sign(m.Z_eQTLGen)
    return m.reset_index(drop=True)

nc = {k for k,v in grpof.items() if v == "NC"}
t2 = {k for k,v in grpof.items() if v == "T2DM"}
cand = {k for k,v in grpof.items() if v == "CAND"}
NTB = {"SERPINH1","RPL13","VAT1","RPL17","RPL7A","TPM4"}

ARMS = [("Primary (32 genes / 96 pairs)", build(nc | t2)),
        ("Anchor set (34 genes, incl. 4 non-testbed)", build(nc | t2 | NTB)),
        ("Full testbed (all groups w/ both-source Z)", build(nc | t2 | cand))]

def summarize(df, name):
    per = df.groupby("Gene").agg(k=("Same","sum"), n=("Same","size"))
    genes = per.index.to_numpy(); karr = per.k.to_numpy().astype(float); narr = per.n.to_numpy().astype(float)
    K = len(genes); N = int(narr.sum()); kk = int(karr.sum()); p = kk / N
    naive_p = binomtest(kk, N, 0.5).pvalue
    ci = binomtest(kk, N, 0.5).proportion_ci(confidence_level=0.95, method="exact")
    m_bar = narr.mean()
    msb = (narr * (karr/narr - p)**2).sum() / (K - 1)
    msw = (karr * (1 - karr/narr)**2 + (narr - karr) * (0 - (1 - karr/narr))**2).sum() / (N - K)
    icc = max(0.0, (msb - msw) / (msb + (m_bar - 1) * msw))
    deff = 1 + (m_bar - 1) * icc
    resid = karr - p * narr
    var_cl = (K / (K - 1)) * (resid**2).sum() / (N**2)
    se_cl = np.sqrt(var_cl); z = (p - 0.5) / se_cl; p_cl = 2 * (1 - norm.cdf(abs(z)))

    gi = {g: np.flatnonzero(df.Gene.to_numpy() == g) for g in genes}
    idxs = [gi[g] for g in genes]
    Zt = df.Z_GTEx.to_numpy(); Ze = df.Z_eQTLGen.to_numpy(); Same = df.Same.to_numpy()
    n_k = narr.astype(int)
    props = np.empty(B); rhos = np.empty(B)
    for b in range(B):
        ch = RNG.integers(0, K, K)
        rows = np.concatenate([idxs[i] for i in ch])
        props[b] = Same[rows].mean()
        rhos[b] = spearmanr(Zt[rows], Ze[rows]).statistic
    ci_lo, ci_hi = np.percentile(props, [2.5, 97.5])
    p_boot = 2 * min((props <= 0.5).mean(), (props >= 0.5).mean())

    perm = np.empty(B)
    trait_arr = df.Trait.to_numpy()
    for b in range(B):
        Ze_p = Ze.copy()
        for t in TRAITS:
            m = trait_arr == t
            Ze_p[m] = RNG.permutation(Ze[m])
        perm[b] = (np.sign(Zt) == np.sign(Ze_p)).mean()
    p_perm = (perm >= p).mean()

    rho_o = spearmanr(Zt, Ze).statistic
    rlo, rhi = np.nanpercentile(rhos, [2.5, 97.5])
    return dict(name=name, genes=K, pairs=N, consistent=kk, pct=100*p,
                naive_p=float(naive_p), ci_naive=[100*float(ci.low), 100*float(ci.high)],
                icc=float(icc), deff=float(deff), se_cl=float(se_cl), z=float(z), p_cluster=float(p_cl),
                boot_ci=[float(ci_lo*100), float(ci_hi*100)], p_boot=float(p_boot), p_perm=float(p_perm),
                rho=float(rho_o), rho_ci=[float(rlo), float(rhi)])

res = [summarize(df, nm) for nm, df in ARMS]
print("=" * 104)
for r in res:
    print(f"\n### {r['name']}  [{r['genes']} genes / {r['pairs']} pairs]")
    print(f"  方向一致率              : {r['consistent']}/{r['pairs']} = {r['pct']:.1f}%")
    print(f"  naive 精确二项 P        : {r['naive_p']:.4f}   95%CI(exact) {r['ci_naive'][0]:.1f}-{r['ci_naive'][1]:.1f}%")
    print(f"  ICC / 设计效应 DEFF     : {r['icc']:.3f} / {r['deff']:.3f}")
    print(f"  簇稳健 SE / z / P       : {r['se_cl']:.4f} / {r['z']:.2f} / {r['p_cluster']:.4f}")
    print(f"  基因层 bootstrap 95%CI  : {r['boot_ci'][0]:.1f}-{r['boot_ci'][1]:.1f}%   bootstrap P = {r['p_boot']:.4f}")
    print(f"  基因标签置换 P          : {r['p_perm']:.4f}")
    print(f"  Spearman rho            : {r['rho']:.3f}   rho 簇 bootstrap 95%CI {r['rho_ci'][0]:.3f}-{r['rho_ci'][1]:.3f}")

pass
pass
print("\nsaved: s1_results.json")
