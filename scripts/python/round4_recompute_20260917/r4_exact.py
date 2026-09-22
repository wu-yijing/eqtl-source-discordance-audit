# -*- coding: utf-8 -*-
"""
R4：严格复现 analysis_round3.py 的「双臂配对基因簇自助法」，但只换官方 Z 输入。
估计量（逐字复刻原脚本）：
  GT[i]  = SI Table S2 的 FDR_q_ACAT_O < 0.05
  EQ[i]  = SI Table S17（testbed eQTLGen）的 FDR 调用
  HKGT   = 对 SI Table S5 的 ACAT-O P 在每个表型内重算 BH
  HKEQ   = SI Table S14（HK eQTLGen）的 BH q < 0.05
  gr(gs)= mean over gene×phenotype ；D1 = gr(HK,GT) − gr(CAND,GT)（housekeeping − candidate）
  paired bootstrap: hkC 与 caC 各自有放回重抽 → 同一次重抽里算 d1、d2
  cov → se_dep = sqrt(v1+v2−2cov) ；r12 = cov/sqrt(v1·v2)
  pooled SE（解析）= sqrt(1/(w1+w2))，S1c/S2c 用 len(gene)*3 作分母
  cluster-only / cluster+cross-arm 的合并 SE
"""
import sys, io, os, csv, json, math, random
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

BASE = r"E:\workbuddy\BMC Genomics投稿资料\定稿资料"
PKG = os.path.join(BASE, "官方Z重建_数据包_20260916")
GRPJ = r"E:\workbuddy\2026-09-16-20-45-42\rewrite\groups.json"
PH = ("DR", "DN", "DPN")

g = json.load(open(GRPJ, encoding="utf-8"))
grp = dict(g["GRP"]); grp.update({x: "HK" for x in g["hk"]})

def fdr_bh(vals):
    p = np.asarray(vals, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank in range(n, 0, -1):
        i = o[rank - 1]; v = min(prev, p[i] * n / rank); q[i] = v; prev = v
    return q

# ---- GT：SI Table S2 的 FDR_q_ACAT_O ----
GT = {}
for r in csv.DictReader(open(os.path.join(PKG, "new_TableS2.csv"), encoding="utf-8")):
    v = (r["FDR_q_ACAT_O"] or "").strip()
    if v == "": continue
    GT.setdefault(r["Gene"], {})[r["Phenotype"]] = 1 if float(v) < 0.05 else 0

# ---- EQ / HKEQ：官方 eQTLGen 逐基因 ----
EQ, HKEQ = {}, {}
for r in csv.DictReader(open(os.path.join(PKG, "new_TableS17_S14.csv"), encoding="utf-8")):
    q = (r["BH_q"] or "").strip()
    if q == "": continue
    hit = 1 if float(q) < 0.05 else 0
    if r["Group"] == "HK":
        HKEQ.setdefault(r["Gene"], {})[r["Phenotype"]] = hit
    else:
        EQ.setdefault(r["Gene"], {})[r["Phenotype"]] = hit

# ---- HKGT：对 SI Table S5 的 ACAT-O P 逐表型重算 BH ----
HKGT = {}
S5 = {}
for r in csv.DictReader(open(os.path.join(PKG, "new_TableS5_GTExPart.csv"), encoding="utf-8")):
    v = (r["P_ACAT_O"] or "").strip()
    if v == "": continue
    S5.setdefault(r["Gene"], {})[r["Phenotype"]] = float(v)
for k in PH:
    names = [x for x in S5 if k in S5[x]]
    qs = fdr_bh([S5[x][k] for x in names])
    for nm, qq in zip(names, qs):
        HKGT.setdefault(nm, {})[k] = 1 if qq < 0.05 else 0

hkG = sorted(x for x in HKGT if len(HKGT[x]) == 3)
hkE = sorted(x for x in HKEQ if len(HKEQ[x]) == 3)
candG = sorted(x for x in GT if grp.get(x) == "Candidate" and len(GT[x]) == 3)
candE = sorted(x for x in EQ if grp.get(x) == "Candidate" and len(EQ[x]) == 3)
print("可分析基因：HK_GTEx=%d HK_eQTLGen=%d CAND_GTEx=%d CAND_eQTLGen=%d"
      % (len(hkG), len(hkE), len(candG), len(candE)))

def gr(gs, tab): return sum(tab[x][t] for x in gs for t in PH) / (len(gs) * 3)
def sel(gs, tab):
    v = [tab[x][t] for x in gs for t in PH]; return sum(v), len(v)

D1 = gr(hkG, HKGT) - gr(candG, GT)
D2 = gr(hkE, HKEQ) - gr(candE, EQ)
s1a, n1a = sel(hkG, HKGT); s1b, n1b = sel(candG, GT)
s2a, n2a = sel(hkE, HKEQ); s2b, n2b = sel(candE, EQ)
S1 = math.sqrt((s1a/n1a)*(1-s1a/n1a)/n1a + (s1b/n1b)*(1-s1b/n1b)/n1b)
S2 = math.sqrt((s2a/n2a)*(1-s2a/n2a)/n2a + (s2b/n2b)*(1-s2b/n2b)/n2b)
w1, w2 = 1/S1**2, 1/S2**2
Dp = (w1*D1 + w2*D2)/(w1+w2); Sp = math.sqrt(1/(w1+w2)); Q = w1*(D1-Dp)**2 + w2*(D2-Dp)**2
print("独立 IV 合并：D1=%+.2fpp  D2=%+.2fpp  pooled=%+.2fpp  Q=%.3f  SE=%.2fpp"
      % (100*D1, 100*D2, 100*Dp, Q, 100*Sp))

hkC = sorted(set(hkG) & set(hkE)); caC = sorted(set(candG) & set(candE))
print("配对重抽共同基因集：HK=%d  CAND=%d" % (len(hkC), len(caC)))
assert hkC and caC

def paired(seed):
    rnd = random.Random(seed)
    B = 10000; bs1 = []; bs2 = []
    for _ in range(B):
        aH = [rnd.choice(hkC) for _ in hkC]; aC = [rnd.choice(caC) for _ in caC]
        bs1.append(gr(aH, HKGT) - gr(aC, GT))
        bs2.append(gr(aH, HKEQ) - gr(aC, EQ))
    m1 = sum(bs1)/B; m2 = sum(bs2)/B
    v1 = sum((x-m1)**2 for x in bs1)/(B-1); v2 = sum((x-m2)**2 for x in bs2)/(B-1)
    cov = sum((bs1[i]-m1)*(bs2[i]-m2) for i in range(B))/(B-1)
    return m1, m2, v1, v2, cov, bs1, bs2

for seed in (20260915, 1):
    m1, m2, v1, v2, cov, bs1, bs2 = paired(seed)
    r12 = cov/math.sqrt(v1*v2)
    se_ind = math.sqrt(v1+v2); se_dep = math.sqrt(v1+v2-2*cov)
    print(f"\n--- seed {seed} ---")
    print(f"  配对重抽均值 d1={100*m1:+.2f}pp  d2={100*m2:+.2f}pp")
    print(f"  r(d1,d2) = {r12:.3f}   SE_indep={100*se_ind:.2f}pp  SE_paired={100*se_dep:.2f}pp  比值={se_dep/se_ind:.3f}")
    d1c = gr(hkC, HKGT) - gr(caC, GT); d2c = gr(hkC, HKEQ) - gr(caC, EQ)
    S1c = math.sqrt(gr(hkC, HKGT)*(1-gr(hkC, HKGT))/(len(hkC)*3) + gr(caC, GT)*(1-gr(caC, GT))/(len(caC)*3))
    S2c = math.sqrt(gr(hkC, HKEQ)*(1-gr(hkC, HKEQ))/(len(hkC)*3) + gr(caC, EQ)*(1-gr(caC, EQ))/(len(caC)*3))
    w1a, w2a = 1/S1c**2, 1/S2c**2
    se_an = 1/math.sqrt(w1a+w2a)
    w1b, w2b = 1/v1, 1/v2
    se_cl = math.sqrt((w1b**2*v1 + w2b**2*v2))/(w1b+w2b)
    se_cc = math.sqrt((w1b**2*v1 + w2b**2*v2 + 2*w1b*w2b*cov))/(w1b+w2b)
    Dpc = (w1b*d1c + w2b*d2c)/(w1b+w2b)
    print(f"  共同基因集点估计：d1c={100*d1c:+.2f}pp  d2c={100*d2c:+.2f}pp")
    print(f"  合并 SE：解析 {100*se_an:.2f}pp | 仅簇 {100*se_cl:.2f}pp | 簇+跨臂 {100*se_cc:.2f}pp")
    print(f"  放大倍数：仅簇 ×{se_cl/se_an:.2f}，簇+跨臂 ×{se_cc/se_an:.2f}")
    lo, hi = 100*(Dpc-1.645*se_cc), 100*(Dpc+1.645*se_cc)
    print(f"  簇+跨臂 90% CI = [{lo:+.2f}, {hi:+.2f}] pp（±15 内：{lo > -15 and hi < 15}）")
    out = dict(seed=seed, hk_common=len(hkC), cand_common=len(caC),
               r_paired=round(r12, 3), cov_pp2=round(cov*1e4, 4),
               d1c_pp=round(100*d1c, 2), d2c_pp=round(100*d2c, 2),
               se_analytic_pp=round(100*se_an, 2), se_cluster_pp=round(100*se_cl, 2),
               se_cluster_crossarm_pp=round(100*se_cc, 2),
               infl_cluster=round(se_cl/se_an, 2), infl_cluster_crossarm=round(se_cc/se_an, 2),
               ci90_cluster_crossarm_pp=[round(lo, 2), round(hi, 2)],
               D1_pp=round(100*D1, 2), D2_pp=round(100*D2, 2), pooled_pp=round(100*Dp, 2),
               se_pooled_pp=round(100*Sp, 2), Q=round(Q, 3))
    json.dump(out, open(os.path.join(PKG, f"r4_dual_bootstrap_seed{seed}.json"), "w",
                        encoding="utf-8"), ensure_ascii=False, indent=1)
print("\nsaved r4_dual_bootstrap_seed*.json")
