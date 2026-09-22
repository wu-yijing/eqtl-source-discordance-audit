# -*- coding: utf-8 -*-
"""Robustness checks for the 2x2 TWAS discordance decomposition.

Pre-empts the 'post-hoc segmentation' criticism:
  The three decomposition arms use DIFFERENT gene x trait coverage
  (1D-panel n=162, 1D-tissue n=150, 2D n=201, anchor n=102). A reviewer could
  argue the apparent 'panel-dominance' is an artifact of comparing different
  gene subsets rather than a real effect of the eQTL-panel / tissue dimension.

Three checks:
  RC-A  Common-universe restriction  -> recompute ALL arms on the SAME
        intersection of gene x trait keys present in every column
        (GTEx WB, GTEx NT, GTEx Multi, eQTLGen). If pattern survives, the
        decomposition is not driven by differential gene coverage.
  RC-B  Trait-stratified             -> recompute arms within DN / DR / DPN
        separately. If panel-dominance holds in each stratum, it is not a
        pooling artifact across heterogeneous traits.
  RC-C  Bootstrap CIs                -> 1000 resamples of the common-universe
        arms; report 95% CI per arm and for (rho_tissue - rho_panel) to show
        the dimension effect is beyond sampling noise.
"""
import csv, os, json, random
from scipy import stats
import numpy as np

BASE   = r"E:/workbuddy/GigaScience投稿/投稿所需资料/additional_files/data"
TABLES = r"E:/workbuddy/GigaScience投稿/2026-07-04_figures_supplement备份/tables"

TRAITS = ["DR", "DN", "DPN"]
random.seed(20260716)
np.random.seed(20260716)

# ---------- load GTEx (WB / NT / Multi) from TableS2 ----------
gtex = {}   # (Gene, Trait) -> dict(WB, NT, Multi)
with open(os.path.join(TABLES, "TableS2_gtex_baseline_results.csv"),
          newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        g = row["Gene"]
        for t in TRAITS:
            def z(col):
                try: return float(row[col])
                except (ValueError, TypeError): return None
            gtex[(g, t)] = {"WB": z(f"Z_Whole_Blood_{t}"),
                            "NT": z(f"Z_Nerve_Tibial_{t}"),
                            "Multi": z(f"Z_Multi_{t}")}

# ---------- load eQTLGen (whole blood) ----------
eq = {}
with open(os.path.join(BASE, "eqtlgen_spredixcan_results.csv"),
          newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        try:
            eq[(row["Gene"], row["Trait"])] = float(row["Z_eQTLGen"])
        except (ValueError, TypeError):
            pass

def spearman_pairs(pairs):
    """pairs: list of (a, b). Returns (rho, p, consistency%, n)."""
    pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
    if len(pairs) < 4:
        return (None, None, None, len(pairs))
    za = [p[0] for p in pairs]; zb = [p[1] for p in pairs]
    rho, p = stats.spearmanr(za, zb)
    cons = sum(1 for a, b in pairs if (a > 0) == (b > 0)) / len(pairs)
    return (round(float(rho), 3), float(p), round(cons*100, 1), len(pairs))

def arms_on(keys):
    """Given a set of (Gene,Trait) keys, return the three arms' (rho,p,cons,n)."""
    panel  = [(gtex[k]["WB"], eq.get(k)) for k in keys]          # same tissue, diff panel
    tissue = [(gtex[k]["WB"], gtex[k]["NT"]) for k in keys]      # same panel, diff tissue
    d2     = [(gtex[k]["Multi"], eq.get(k)) for k in keys]       # diff tissue + diff panel
    return {"1D_panel": spearman_pairs(panel),
            "1D_tissue": spearman_pairs(tissue),
            "2D": spearman_pairs(d2)}

def fmt(res):
    rho, p, cons, n = res
    if rho is None: return f"n={n} (insufficient)"
    return f"n={n}  rho={rho:.3f}  cons={cons:.1f}%"

print("="*74)
print("RC-A  COMMON-UNIVERSE RESTRICTION  (intersection of WB, NT, Multi, eQTLGen)")
print("="*74)
common = set()
for k in gtex:
    if gtex[k]["WB"] is None: continue
    if gtex[k]["NT"] is None: continue
    if gtex[k]["Multi"] is None: continue
    if k not in eq: continue
    common.add(k)
print(f"  common gene x trait keys : {len(common)}  "
      f"(genes={len(set(g for g,t in common))})")
rca = arms_on(common)
for name, res in rca.items():
    print(f"  [{name:<10}] {fmt(res)}")

print("\n" + "="*74)
print("RC-B  TRAIT-STRATIFIED  (per-trait recomputation)")
print("="*74)
rcb = {}
for t in TRAITS:
    keys_t = {k for k in common if k[1] == t}
    rcb[t] = arms_on(keys_t)
    print(f"  -- Trait {t}  (common-universe genes only, n_keys={len(keys_t)}) --")
    for name, res in rcb[t].items():
        print(f"     [{name:<10}] {fmt(res)}")

print("\n" + "="*74)
print("RC-C  BOOTSTRAP 95% CI  (common-universe arms, 1000 resamples)")
print("="*74)
B = 1000
arm_pairs = {
    "1D_panel":  [(gtex[k]["WB"], eq.get(k)) for k in common],
    "1D_tissue": [(gtex[k]["WB"], gtex[k]["NT"]) for k in common],
    "2D":        [(gtex[k]["Multi"], eq.get(k)) for k in common],
}
boot_rho = {a: [] for a in arm_pairs}
diff_tp = []   # rho_tissue - rho_panel
for _ in range(B):
    idx = np.random.randint(0, len(common), size=len(common))
    for a, base in arm_pairs.items():
        samp = [base[i] for i in idx]
        rho, _, _, n = spearman_pairs(samp)
        if rho is not None:
            boot_rho[a].append(rho)
    # paired difference uses same resample index for both arms
    samp_p = [arm_pairs["1D_panel"][i] for i in idx]
    samp_t = [arm_pairs["1D_tissue"][i] for i in idx]
    rp, _, _, _ = spearman_pairs(samp_p)
    rt, _, _, _ = spearman_pairs(samp_t)
    if rp is not None and rt is not None:
        diff_tp.append(rt - rp)

def ci(lst):
    arr = np.array(lst)
    return (round(float(np.percentile(arr, 2.5)), 3),
            round(float(np.percentile(arr, 97.5)), 3),
            round(float(np.mean(arr)), 3))

print("  arm            point-rho   95% CI (bootstrap)")
for a in ["1D_panel", "1D_tissue", "2D"]:
    lo, hi, mean = ci(boot_rho[a])
    print(f"  {a:<12} {rca[a][0]:>8.3f}   [{lo:>6.3f}, {hi:>6.3f}]")
lo, hi, mean = ci(diff_tp)
print(f"\n  rho_tissue - rho_panel  point = {rca['1D_tissue'][0]-rca['1D_panel'][0]:+.3f}"
      f"   95% CI [{lo:+.3f}, {hi:+.3f}]   (excludes 0 => panel effect real)")

# ---------- save ----------
out = {
    "common_universe": {"n_keys": len(common),
                        "n_genes": len(set(g for g, t in common)),
                        "arms": {k: {"rho": v[0], "p": v[1],
                                     "consistency": v[2], "n": v[3]} for k, v in rca.items()}},
    "trait_stratified": {t: {k: {"rho": v[0], "p": v[1],
                                 "consistency": v[2], "n": v[3]} for k, v in rcb[t].items()}
                          for t in TRAITS},
    "bootstrap": {
        "B": B,
        "ci": {a: {"lo": ci(boot_rho[a])[0], "hi": ci(boot_rho[a])[1],
                   "mean": ci(boot_rho[a])[2]} for a in arm_pairs},
        "rho_tissue_minus_rho_panel": {"point": round(rca['1D_tissue'][0]-rca['1D_panel'][0], 3),
                                       "ci_lo": lo, "ci_hi": hi},
    },
}
with open(r"E:/workbuddy/2026-07-15-20-17-22/robustness_results.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
print("\nSaved -> robustness_results.json")
