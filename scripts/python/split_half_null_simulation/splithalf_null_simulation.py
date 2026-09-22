# -*- coding: utf-8 -*-
"""Split-half null simulation for the 96-pair primary arm (per-phenotype available SNP sets)."""
import numpy as np, csv, io, os, sys, math, json, time, traceback

BASE = u"E:\\workbuddy\\2026-09-13-06-39-47"
BED = os.path.join(BASE, u"ld", u"g1000_eur.bed")
BIM_MAP = os.path.join(BASE, u"_bim_map.csv")
CIS = os.path.join(BASE, u"_eqtlgen_cis.csv")
BPS = 126
NS = 503
N_PART = 200
SEED = 20260910
COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def decode_a2(buf, n):
    b = np.frombuffer(buf, dtype=np.uint8)
    codes = np.empty((b.size, 4), dtype=np.uint8)
    codes[:, 0] = b & 3
    codes[:, 1] = (b >> 2) & 3
    codes[:, 2] = (b >> 4) & 3
    codes[:, 3] = (b >> 6) & 3
    codes = codes.reshape(-1)[:n]
    out = np.full(n, -1, dtype=np.int8)
    out[codes == 0] = 0
    out[codes == 2] = 1
    out[codes == 3] = 2
    return out


def sgn(eff, A1, A2):
    if eff == A2:
        return 1
    if eff == A1:
        return -1
    c = COMP.get(eff)
    if c == A2:
        return 1
    if c == A1:
        return -1
    return None


def main():
    t0 = time.time()
    log = []
    out = {}
    try:
        bim = {}
        with io.open(BIM_MAP, encoding='utf-8') as f:
            for r in csv.DictReader(f):
                bim[r['rsid']] = (int(r['snp_index']), r['a1'].upper(), r['a2'].upper())
        gwas = {}
        for ph in ["DR", "DN", "DPN"]:
            d = {}
            with io.open(os.path.join(BASE, u"_gwas_%s.csv" % ph), encoding='utf-8') as f:
                for r in csv.DictReader(f):
                    d[r['rsid']] = (float(r['Z']), r['alt'].upper())
            gwas[ph] = d
        cis = {}
        with io.open(CIS, encoding='utf-8') as f:
            for r in csv.DictReader(f):
                cis.setdefault(r['gene'], []).append(r)
        genes = sorted(cis.keys())

        need = set()
        for g in genes:
            for r in cis[g]:
                if r['snp'] in bim:
                    need.add(r['snp'])
        need = sorted(need)
        idxs = np.array([bim[s][0] for s in need])
        order = np.argsort(idxs)
        GENO = np.zeros((len(need), NS), dtype=np.int8)
        with open(BED, 'rb') as f:
            f.read(3)
            for k in order:
                f.seek(3 + int(idxs[k]) * BPS)
                GENO[k] = decode_a2(f.read(BPS), NS)
        for k in range(GENO.shape[0]):
            m = GENO[k] == -1
            if m.any():
                GENO[k, m] = GENO[k][~m].mean() if (~m).any() else 0
        row = {s: i for i, s in enumerate(need)}
        log.append(u"geno %s  (%.0fs)" % (GENO.shape, time.time() - t0))

        rng = np.random.default_rng(SEED)
        records = []
        for gene in genes:
            keep = []
            for r in cis[gene]:
                rid = r['snp']
                if rid not in row:
                    continue
                A1, A2 = bim[rid][1], bim[rid][2]
                sw = sgn(r['assessed'].upper(), A1, A2)
                if sw is None:
                    continue
                keep.append((rid, float(r['zscore']) * sw, A1, A2))
            if len(keep) < 4:
                log.append(u"  skip %s (%d)" % (gene, len(keep)))
                continue
            gi = np.array([row[k[0]] for k in keep])
            X = GENO[gi].astype(np.float32)
            w = np.array([k[1] for k in keep], dtype=np.float32)
            R = np.corrcoef(X).astype(np.float32)
            R = np.nan_to_num(R)
            np.fill_diagonal(R, 1.0)
            n = len(w)
            for ph in ["DR", "DN", "DPN"]:
                ia, zl = [], []
                for i, (rid, _, A1, A2) in enumerate(keep):
                    g = gwas[ph].get(rid)
                    if g is None:
                        continue
                    zval, alt = g
                    t = sgn(alt, A1, A2)
                    if t is None:
                        continue
                    ia.append(i); zl.append(t * zval)
                na = len(ia)
                if na < 4:
                    continue
                ia = np.array(ia)
                zf = np.zeros(n, dtype=np.float32)
                zf[ia] = zl
                u0 = np.zeros(n, dtype=np.float32)
                u0[ia] = w[ia]
                full_Z = float(u0 @ zf) / math.sqrt(max(float(u0 @ (R @ u0)), 1e-12))
                bools = []
                for rep in range(N_PART):
                    perm = rng.permutation(na)
                    h1 = ia[perm[: na // 2]]
                    h2 = ia[perm[na // 2:]]
                    Zs = []
                    for h in (h1, h2):
                        u = np.zeros(n, dtype=np.float32)
                        u[h] = w[h]
                        Zs.append(float(u @ zf) / math.sqrt(max(float(u @ (R @ u)), 1e-12)))
                    bools.append((Zs[0] > 0) == (Zs[1] > 0))
                records.append({"gene": gene, "ph": ph, "bools": bools,
                                "full_Z": full_Z, "n_cis": n, "n_used": na})
            log.append(u"  %s n_cis=%d recs=%d (%.0fs)" % (gene, n, len(records), time.time() - t0))

        out["records"] = records
        if records:
            allrate = [float(np.mean([r["bools"][k] for r in records])) for k in range(N_PART)]
            out["overall"] = {"mean": float(np.mean(allrate)),
                              "lo": float(np.percentile(allrate, 2.5)),
                              "hi": float(np.percentile(allrate, 97.5)),
                              "median": float(np.median(allrate)),
                              "n_pairs": len(records)}
            byph = {}
            for ph in ["DR", "DN", "DPN"]:
                sub = [r for r in records if r["ph"] == ph]
                if sub:
                    v = [float(np.mean([r["bools"][k] for r in sub])) for k in range(N_PART)]
                    byph[ph] = {"mean": float(np.mean(v)), "lo": float(np.percentile(v, 2.5)),
                                "hi": float(np.percentile(v, 97.5)), "n_pairs": len(sub)}
            out["by_phenotype"] = byph
            # gene-level per-pair concordance
            out["per_pair"] = [{"gene": r["gene"], "ph": r["ph"],
                                "conc": float(np.mean(r["bools"])), "n_used": r["n_used"]}
                               for r in records]
        out["log"] = log
    except Exception:
        log.append(u"EXCEPTION:\n" + traceback.format_exc())
        out = {"log": log}
    with io.open(os.path.join(BASE, u"_splithalf_result.json"), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1)
    summ = [u"\n".join(out.get("log", [])), u"",
            u"overall: %s" % out.get("overall"),
            u"by_phenotype: %s" % out.get("by_phenotype")]
    with io.open(os.path.join(BASE, u"_splithalf_summary.txt"), 'w', encoding='utf-8') as f:
        f.write(u"\n".join(summ))


main()
