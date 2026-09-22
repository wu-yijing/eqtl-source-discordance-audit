# -*- coding: utf-8 -*-
"""
R1 + R4 重算（官方 Z 官方 MetaXcan 化之后）
  R1 : Table S7 表注 —— 4 组 × 2 源 的中位 |Z|(IQR) + 固定阈值下的 Fisher OR/P 全表
  R4 : 双臂配对基因簇自助法 —— r、双臂可测基因数、cluster SE、膨胀因子、90% CI
输入：官方 Z 重建数据包 + groups.json ；不依赖任何旧实现。
"""
import sys, io, os, csv, json, math, collections
import numpy as np
from scipy import stats
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

BASE = r"E:\workbuddy\BMC Genomics投稿资料\定稿资料"
PKG = os.path.join(BASE, "官方Z重建_数据包_20260916")
GRPJ = r"E:\workbuddy\2026-09-16-20-45-42\rewrite\groups.json"
OUT = os.path.join(BASE, "官方Z重建_数据包_20260916")
TRAITS = ["DR", "DN", "DPN"]
GROUP_ORDER = ["Candidate", "NonCandidate", "T2DMcontrol"]
LABEL = {"Candidate": "candidate", "NonCandidate": "non-candidate",
         "T2DMcontrol": "T2DM control", "HK": "housekeeping"}

g = json.load(open(GRPJ, encoding="utf-8"))
GRP = g["GRP"]; HK = g["hk"]
GOF = dict(GRP); GOF.update({x: "HK" for x in HK})

# ---------- 官方 Z ----------
gtex = collections.defaultdict(dict)          # gene -> trait -> P_ACAT_O
for r in csv.DictReader(open(os.path.join(PKG, "new_TableS2.csv"), encoding="utf-8")):
    v = (r["P_ACAT_O"] or "").strip()
    if v:
        gtex[r["Gene"]][r["Phenotype"]] = float(v)
# housekeeping 的 GTEx ACAT-O 在 Table S5 的数据包里
for r in csv.DictReader(open(os.path.join(PKG, "new_TableS5_GTExPart.csv"), encoding="utf-8")):
    v = (r["P_ACAT_O"] or "").strip()
    if v:
        gtex[r["Gene"]][r["Phenotype"]] = float(v)

eqg = collections.defaultdict(dict)           # gene -> trait -> Z
for r in csv.DictReader(open(os.path.join(PKG, "new_TableS17_S14.csv"), encoding="utf-8")):
    eqg[r["Gene"]][r["Phenotype"]] = float(r["Z_eQTLGen"])
    if r["Gene"] not in GOF:
        GOF[r["Gene"]] = r["Group"]

def absz_from_p(p):
    return abs(stats.norm.ppf(max(p, 1e-300) / 2.0))

# =====================================================================
print("=" * 78)
print("R1(a)  中位 |Z| (IQR) —— 4 组 × 2 源（官方）")
print("=" * 78)
RES = {}
for src in ("GTEx", "eQTLGen"):
    print(f"\n--- {src} ---")
    for grp in GROUP_ORDER + ["HK"]:
        genes = [x for x in GOF if GOF[x] == grp]
        vals = []
        for x in genes:
            for t in TRAITS:
                if src == "GTEx":
                    if t in gtex.get(x, {}):
                        vals.append(absz_from_p(gtex[x][t]))
                else:
                    if t in eqg.get(x, {}):
                        vals.append(abs(eqg[x][t]))
        vals = np.array(vals)
        if len(vals) == 0:
            continue
        q1, q3 = np.percentile(vals, [25, 75])
        RES[(src, grp)] = (np.median(vals), q1, q3, len(vals))
        print(f"   {LABEL[grp]:<14} n={len(vals):>3}  median {np.median(vals):.2f}  "
              f"IQR {q1:.2f}\u2013{q3:.2f}")
        print(f"        gene count = {len(genes)}  pairs = {len(vals)}")

# =====================================================================
print()
print("=" * 78)
print("R1(b)  Fisher 精确检验：各组 vs housekeeping（三表型合并）")
print("=" * 78)
THR = [(0.05, "p < 0.05"), (0.01, "p < 0.01"), (0.00385, "p < 0.00385"),
       (0.00179, "p < 0.00179")]
FISH = {}
for src in ("GTEx", "eQTLGen"):
    print(f"\n--- {src} ---")
    pool = {}
    for grp in GROUP_ORDER + ["HK"]:
        genes = [x for x in GOF if GOF[x] == grp]
        for thr, lab in THR:
            k = n = 0
            for x in genes:
                for t in TRAITS:
                    if src == "GTEx":
                        if t not in gtex.get(x, {}):
                            continue
                        pv = gtex[x][t]
                    else:
                        if t not in eqg.get(x, {}):
                            continue
                        pv = 2 * stats.norm.sf(abs(eqg[x][t]))
                    n += 1
                    if pv < thr:
                        k += 1
            pool[(grp, lab)] = (k, n)
    for thr, lab in THR:
        line = []
        for grp in GROUP_ORDER:
            a, na = pool[(grp, lab)]
            c, nc = pool[("HK", lab)]
            b = na - a; d = nc - c
            orr, p = stats.fisher_exact([[a, b], [c, d]])
            FISH[(src, lab, grp)] = (orr, p, a, na, c, nc)
            line.append(f"{LABEL[grp]}: OR={orr:.2f} P={p:.3f}")
        print(f"   {lab:<14} " + " | ".join(line))

# =====================================================================
print()
print("=" * 78)
print("R4  双臂配对基因簇自助法（B = 10,000）")
print("=" * 78)
# 双臂皆可测的基因（GTEx ACAT-O 有值 且 eQTLGen 有值）
def pairset(grp):
    out = []
    for x in GOF:
        if GOF[x] != grp:
            continue
        if not any(t in gtex.get(x, {}) for t in TRAITS):
            continue
        if not any(t in eqg.get(x, {}) for t in TRAITS):
            continue
        out.append(x)
    return sorted(out)

hk_both = pairset("HK")
cand_both = pairset("Candidate")
print(f"  双臂皆可测：housekeeping {len(hk_both)} 个，candidate {len(cand_both)} 个")

def indicators(genes, src):
    """返回 (n_pairs, k_hits) 及每基因的 (k, n) 列表"""
    per = {}
    for x in genes:
        k = n = 0
        for t in TRAITS:
            if src == "GTEx":
                if t not in gtex.get(x, {}):
                    continue
                pv = gtex[x][t]
            else:
                if t not in eqg.get(x, {}):
                    continue
                pv = 2 * stats.norm.sf(abs(eqg[x][t]))
            n += 1
            if pv < 0.05:
                k += 1
        per[x] = (k, n)
    return per

GT_hk, GT_ca = indicators(hk_both, "GTEx"), indicators(cand_both, "GTEx")
EQ_hk, EQ_ca = indicators(hk_both, "eQTLGen"), indicators(cand_both, "eQTLGen")

def rate(per, genes):
    k = sum(per[x][0] for x in genes); n = sum(per[x][1] for x in genes)
    return k, n, (k / n if n else float("nan"))

def analytic(k1, n1, k2, n2):
    p1, p2 = k1 / n1, k2 / n2
    return (p1 - p2) * 100, math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2) * 100

k1, n1, r1 = rate(GT_ca, cand_both); k2, n2, r2 = rate(GT_hk, hk_both)
k3, n3, r3 = rate(EQ_ca, cand_both); k4, n4, r4 = rate(EQ_hk, hk_both)
d1, s1 = analytic(k1, n1, k2, n2)
d2, s2 = analytic(k3, n3, k4, n4)
w1, w2 = 1 / s1**2, 1 / s2**2
D = (w1 * d1 + w2 * d2) / (w1 + w2); S = math.sqrt(1 / (w1 + w2))
print(f"  GTEx 臂  cand {k1}/{n1} vs HK {k2}/{n2} -> {d1:+.2f} pp (SE {s1:.2f})")
print(f"  eQTLGen臂 cand {k3}/{n3} vs HK {k4}/{n4} -> {d2:+.2f} pp (SE {s2:.2f})")
print(f"  合并（解析独立对）D = {D:+.2f} pp  S = {S:.2f} pp")

# bootstrap：每次抽基因（有放回），同一基因同时进入两臂
rng = np.random.default_rng(20260916)
B = 10000
pool_genes = sorted(set(hk_both) | set(cand_both))
bg_gt = {x: GT_ca[x] for x in cand_both}; bg_gt_hk = {x: GT_hk[x] for x in hk_both}
bg_eq = {x: EQ_ca[x] for x in cand_both}; bg_eq_hk = {x: EQ_hk[x] for x in hk_both}
bs1, bs2, bpool = [], [], []
for _ in range(B):
    gs = rng.choice(pool_genes, size=len(pool_genes), replace=True)
    ca = [x for x in gs if x in bg_gt]
    hk = [x for x in gs if x in bg_gt_hk]
    if not ca or not hk:
        continue
    a = sum(bg_gt[x][0] for x in ca); b_ = sum(bg_gt[x][1] for x in ca)
    c = sum(bg_gt_hk[x][0] for x in hk); d_ = sum(bg_gt_hk[x][1] for x in hk)
    if not b_ or not d_:
        continue
    x1, y1 = analytic(a, b_, c, d_)
    a = sum(bg_eq[x][0] for x in ca); b_ = sum(bg_eq[x][1] for x in ca)
    c = sum(bg_eq_hk[x][0] for x in hk); d_ = sum(bg_eq_hk[x][1] for x in hk)
    if not b_ or not d_:
        continue
    x2, y2 = analytic(a, b_, c, d_)
    bs1.append(x1); bs2.append(x2)
    bpool.append((w1 * x1 + w2 * x2) / (w1 + w2))
bs1 = np.array(bs1); bs2 = np.array(bs2); bpool = np.array(bpool)
r = float(np.corrcoef(bs1, bs2)[0, 1])
se_cl = float(np.std(bpool, ddof=1))
lo, hi = np.percentile(bpool, [5, 95])
print(f"  有效重抽次数 = {len(bpool)}")
print(f"  两臂率差的相关系数 r = {r:.3f}")
print(f"  cluster SE = {se_cl:.2f} pp ；解析 SE = {S:.2f} pp ；膨胀 ×{se_cl/S:.2f}")
print(f"  percentile 90% CI = [{lo:+.2f}, {hi:+.2f}] pp")
print(f"  是否落在 ±15 pp 内：{lo > -15 and hi < 15}")

json.dump({
  "median_absZ": {f"{k[0]}|{k[1]}": [float(v[0]), float(v[1]), float(v[2]), int(v[3])]
                  for k, v in RES.items()},
  "fisher": {f"{k[0]}|{k[1]}|{k[2]}": [float(v[0]), float(v[1]), int(v[2]), int(v[3]), int(v[4]), int(v[5])]
             for k, v in FISH.items()},
  "dual_bootstrap": {"B_eff": len(bpool), "r": r, "se_cluster": se_cl,
                     "se_analytic": S, "inflation": se_cl / S,
                     "ci90": [float(lo), float(hi)],
                     "hk_both": len(hk_both), "cand_both": len(cand_both),
                     "arms": [[d1, s1], [d2, s2]], "D": D, "S": S},
}, open(os.path.join(OUT, "r14_recompute.json"), "w", encoding="utf-8"),
   ensure_ascii=False, indent=1)
print("\nSAVED r14_recompute.json")
