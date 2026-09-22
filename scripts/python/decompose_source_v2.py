# -*- coding: utf-8 -*-
"""Decomposition analysis aligned to the ORIGINAL main conclusion (rho=0.289, 59.5%).
Uses the complete GTEx v8 MASHR multi-tissue arm (Z_Multi) as the primary 2D comparison,
and splits the discordance into tissue-only (1D) vs panel-only (1D) components.
"""
import csv, os, json
from scipy import stats

# Use repository-relative paths for portability
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
BASE = os.path.join(REPO_ROOT, 'data', 'processed')
ALT  = os.path.join(REPO_ROOT, 'data', 'processed')
TABLES = os.path.join(REPO_ROOT, 'tables')

def load_dict(path, keycols, valcol, cast=float):
    d = {}
    with open(path, newline='', encoding='utf-8-sig') as f:
        r = csv.DictReader(f)
        for row in r:
            key = tuple(row[k] for k in keycols)
            try:
                d[key] = cast(row[valcol])
            except (ValueError, TypeError):
                d[key] = None
    return d

# ---- 2D MAIN ARM: eqtlgen_vs_gtex_comparison.csv (Z_GTEx = GTEx multi-tissue MASHR) ----
cmp_path = os.path.join(ALT, "eqtlgen_vs_gtex_comparison.csv")
rows = []
with open(cmp_path, newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        try:
            ze = float(row["Z_eQTLGen"]); zg = float(row["Z_GTEx"])
        except (ValueError, TypeError):
            continue
        rows.append((row["Gene"], row["Trait"], ze, zg, row["Same_Direction"]))

ze = [r[2] for r in rows]; zg = [r[3] for r in rows]
rho2, p2 = stats.spearmanr(ze, zg)
cons2 = sum(1 for r in rows if (r[2] > 0) == (r[3] > 0)) / len(rows)
n2 = len(rows)

# ---- Load single-tissue + multi-tissue GTEx Z from TableS2 ----
t2 = os.path.join(TABLES, "TableS2_gtex_baseline_results.csv")
gtex = {}  # (Gene, Trait) -> {NT, WB, Multi}
with open(t2, newline='', encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        g = row["Gene"]
        for trait in ["DR", "DN", "DPN"]:
            def z(col):
                try: return float(row[col])
                except: return None
            gtex[(g, trait)] = {
                "NT": z(f"Z_Nerve_Tibial_{trait}"),
                "WB": z(f"Z_Whole_Blood_{trait}"),
                "Multi": z(f"Z_Multi_{trait}"),
            }

# ---- Load eQTLGen Z (whole blood) ----
eq = load_dict(os.path.join(BASE, "eqtlgen_spredixcan_results.csv"),
               ["Gene", "Trait"], "Z_eQTLGen")

def arm(name, getterA, getterB):
    """getterA/B: dict[(Gene,Trait)] -> z or None"""
    pairs = []
    for k in gtex:
        a = getterA(k); b = getterB(k)
        if a is None or b is None: continue
        pairs.append((a, b))
    za = [p[0] for p in pairs]; zb = [p[1] for p in pairs]
    rho, p = stats.spearmanr(za, zb)
    cons = sum(1 for a, b in pairs if (a > 0) == (b > 0)) / len(pairs)
    rev = 1 - cons
    print(f"\n[{name}]")
    print(f"  n pairs      : {len(pairs)}")
    print(f"  Spearman rho : {rho:.3f}  (P={p:.2e})")
    print(f"  consistency  : {cons*100:.1f}%   reversal: {rev*100:.1f}%")
    return dict(name=name, n=len(pairs), rho=round(rho,3), p=p,
                consistency=round(cons*100,1), reversal=round(rev*100,1))

print("="*70)
print("ORIGINAL MAIN CONCLUSION ANCHOR (2D, complete multi-tissue arm)")
print("="*70)
print(f"  n pairs      : {n2}")
print(f"  Spearman rho : {rho2:.3f}  (P={p2:.2e})")
print(f"  consistency  : {cons2*100:.1f}%   reversal: {(1-cons2)*100:.1f}%")
print("  (manuscript reports rho=0.289, 59.5% consistency, 40.5% reversal)")

print("\n" + "="*70)
print("DECOMPOSITION (using complete multi-tissue GTEx arm)")
print("="*70)

# 1D-panel: same tissue (whole blood), different panel/sample size
a_panel = arm("1D-panel : GTEx WB  vs  eQTLGen WB",
              lambda k: gtex[k]["WB"], lambda k: eq.get(k))
# 1D-tissue: same panel (GTEx), different tissue
a_tissue = arm("1D-tissue: GTEx WB  vs  GTEx Nerve_Tibial",
               lambda k: gtex[k]["WB"], lambda k: gtex[k]["NT"])
# 2D: GTEx multi-tissue (complete) vs eQTLGen
a_2d = arm("2D        : GTEx Multi vs  eQTLGen",
           lambda k: gtex[k]["Multi"], lambda k: eq.get(k))

# ---- Summary table ----
summary = {
    "anchor_2D_rho": round(rho2,3), "anchor_2D_consistency": round(cons2*100,1),
    "arms": [a_panel, a_tissue, a_2d]
}
with open(os.path.join(REPO_ROOT, "analyses", "decomposition_results_v2.json"), "w", encoding="utf-8") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)

print("\n" + "="*70)
print("SUMMARY TABLE")
print("="*70)
print(f"{'Arm':<32}{'n':>5}{'rho':>8}{'cons%':>8}{'rev%':>8}")
for a in [a_panel, a_tissue, a_2d]:
    print(f"{a['name']:<32}{a['n']:>5}{a['rho']:>8}{a['consistency']:>8}{a['reversal']:>8}")
print(f"{'ORIGINAL anchor (multi vs eQTLGen)':<32}{n2:>5}{rho2:>8}{cons2*100:>8.1f}{(1-cons2)*100:>8.1f}")
