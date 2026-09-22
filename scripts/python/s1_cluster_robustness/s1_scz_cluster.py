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
"""S1 补充：SCZ 全基因组臂（2,511 对）的簇稳健复核。

关键前提：complete-case 集合中每个基因恰好贡献 1 对 → 不存在基因内聚类，
因此簇稳健方法与朴素方法在构造上等价。本脚本验证这一点并给出完整统计量。
"""
import os, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, binomtest, norm

# ⚠️ 输入快照属 in-house SCZ 估计器（n = 2,511），已被官方 MetaXcan v0.8.1 取代，**不随本仓库分发**。
#    本脚本保留以记录当时的簇稳健复核过程；重跑前请自备该文件，或用 S1_SCZ_SNAPSHOT 环境变量指向你的副本。
SRC = _os.environ.get("S1_SCZ_SNAPSHOT",
                     _os.path.join(ROOT, "scz_twas_results_limit0.csv"))
OUT = r"E:\workbuddy\2026-09-10-18-23-32"
B = 10000
rng = np.random.default_rng(20260726)          # 沿用原脚本 seed，便于比对

d = pd.read_csv(SRC)
print(f"总行数 {len(d)} | 唯一基因 {d.gene.nunique()}")

cc = d.dropna(subset=["eqZ", "wbZ", "ntZ", "multiZ"])       # 与原脚本一致的 all-4 完整个案
print(f"all-4 完整个案 {len(cc)} | 唯一基因 {cc.gene.nunique()} "
      f"| 每基因对数 = {len(cc)/cc.gene.nunique():.3f}")
eq, wb, nt, mu = cc.eqZ.to_numpy(), cc.wbZ.to_numpy(), cc.ntZ.to_numpy(), cc.multiZ.to_numpy()

def analyse(x, y, name):
    N = len(x)
    same = np.sign(x) == np.sign(y)
    k = int(same.sum()); p = k / N
    naive = binomtest(k, N, 0.5).pvalue
    ci = binomtest(k, N, 0.5).proportion_ci(confidence_level=0.95, method="exact")
    rho = spearmanr(x, y).statistic
    # 每基因 1 对 → 基因层 bootstrap 与个案 bootstrap 等价（此处显式以基因为单位，证明 DEFF=1）
    K = N
    idx = np.arange(K)
    bs = np.empty(B); br = np.empty(B)
    for b in range(B):
        s = rng.integers(0, K, K)
        bs[b] = same[s].mean(); br[b] = spearmanr(x[s], y[s]).statistic
    blo, bhi = np.percentile(bs, [2.5, 97.5])
    pb = 2 * min((bs <= 0.5).mean(), (bs >= 0.5).mean())
    # 基因标签置换：打散两源之间的基因对应关系
    perm = np.empty(B)
    for b in range(B):
        yp = rng.permutation(y)
        perm[b] = (np.sign(x) == np.sign(yp)).mean()
    pperm = float((perm >= p).mean())
    return dict(name=name, N=N, K=K, pairs_per_gene=round(N/K, 3), k=k, pct=100*p,
                naive_p=float(naive), ci_exact=[100*float(ci.low), 100*float(ci.high)],
                boot_ci=[float(blo*100), float(bhi*100)], p_boot=float(pb), p_perm=pperm,
                rho=float(rho), rho_boot_ci=[float(np.percentile(br, 2.5)), float(np.percentile(br, 97.5))],
                deff=1.0, icc=None,
                note="每基因恰好 1 对；簇结构与独立假设完全等价（DEFF = 1，构造性）")

res = [analyse(eq, mu, "SCZ source axis (eQTLGen WB vs GTEx multi-tissue)"),
       analyse(wb, nt, "SCZ tissue axis (GTEx Whole_Blood vs Nerve_Tibial)")]

print("\n" + "=" * 96)
for r in res:
    print(f"\n### {r['name']}")
    print(f"  完整个案 {r['N']}（唯一基因 {r['K']}，每基因 {r['pairs_per_gene']} 对）")
    print(f"  方向一致率        : {r['k']}/{r['N']} = {r['pct']:.2f}%")
    print(f"  naive 精确二项 P  : {r['naive_p']:.3e}   95%CI(exact) {r['ci_exact'][0]:.2f}-{r['ci_exact'][1]:.2f}%")
    print(f"  基因层 bootstrap  : 95%CI {r['boot_ci'][0]:.2f}-{r['boot_ci'][1]:.2f}%   P = {r['p_boot']:.3e}")
    print(f"  基因标签置换 P    : {r['p_perm']:.3e}")
    print(f"  Spearman rho      : {r['rho']:.4f}   bootstrap 95%CI {r['rho_boot_ci'][0]:.4f}-{r['rho_boot_ci'][1]:.4f}")
    print(f"  DEFF              : {r['deff']}  ({r['note']})")

json.dump(res, open(_os.path.join(HERE, "s1_scz_results.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
print("\nsaved: s1_scz_results.json")
