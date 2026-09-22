#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
重选 30 个看家基因（真阴性对照集 v2）

背景：v1 的 data/hk_genes.txt 中 7/30 与 104-gene testbed 重叠
（ACTB / GAPDH / RPL13A / RPLP0 / RPS18 属候选组；EEF2 / ENO1 属非候选组），
违反其自身筛选标准第 3 条 "Not present in the study's existing 104-gene panel"。

v2 选择规则（预先设定，与结果无关，避免挑选偏倚）：
  R1 全集：HRT Atlas v1.0 人鼠共有看家基因集（MSigDB HOUNKPE_HOUSEKEEPING_GENES, n=1129）
  R2 剔除：104-gene testbed panel 中的任何基因
  R3 剔除：与 104 panel 任一基因共享 >=3 字符字母前缀的基因（家族级排除）
          例：HSP*(HSPA8/HSP90AB1/HSPB1/HSPD1)、RPS*、RPL*、EEF*、ENO*、ANXA*、
              ATP5*、HNRNP*、PDIA*、ILF*、XRCC*、SLC*、YWHA* ...
          另剔除 MRPS* / MRPL* / MT-* 等线粒体与线粒体核糖体基因
          （理由：104 panel 高度富集核糖体、热休克蛋白与 RNA 结合蛋白家族，
            保留同家族基因会让"真阴性对照"与 testbed 共享家族层面的 cis-eQTL 架构，
            从而使对照失去意义。此规则在抽取时即已固定，早于任何 TWAS 结果产出。）
  R4 剔除：T2DM / 糖尿病并发症 / 胰岛素与 GLP-1 信号通路的已知基因（人工 curated 黑名单）
  R5 要求：GTEx v8 MASHR 模型在 Nerve_Tibial 与 Whole_Blood 两个组织中均存在且 n.snps.in.model >= 1
          （理由：v1 的 30 个基因中仅 18（NT）/12（WB）个 gene-phenotype pair 有有效 Z，
            可测率过低；强制双组织可测可将对照集的信息量最大化）
  R6 抽样：在合格池中按固定随机种子（seed=20260830）简单随机抽取 30 个
"""
import csv, os, re, random, sqlite3, sys, io, json

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

ROOT = r'E:\workbuddy\hk_reselect_20260830'
REPO = r'E:\workbuddy\eqtl-source-discordance-audit'
HRT = r'E:\workbuddy\hrt_common.csv'
MODEL_DIR = r'E:\workbuddy\BMC Genomics投稿资料\DR，DN，DM芬兰原始数据\mashr_eqtl\eqtl\mashr'
TISSUES = ['Nerve_Tibial', 'Whole_Blood']
N_TARGET = 30
SEED = 20260830

# ---------- R1 全集 ----------
hrt = set()
with open(HRT, encoding='utf-8', errors='replace') as f:
    for line in f:
        line = line.strip()
        if not line or line.lower().startswith('mouse'):
            continue
        p = line.split(';')
        if len(p) >= 2 and p[1].strip():
            hrt.add(p[1].strip().upper())
print(f'R1 HRT Atlas human∩mouse HK genes : {len(hrt)}')

# ---------- R2 104 panel ----------
with open(os.path.join(REPO, 'data', 'processed', 'covariate_matrix.csv'), encoding='utf-8') as f:
    panel = {r['Gene'].upper() for r in csv.DictReader(f)}
print(f'R2 104-gene testbed panel         : {len(panel)}')

# ---------- R3 家族黑名单 ----------
EXTRA_FAMILY = re.compile(r'^(MRPS|MRPL|MT-|MTRNR|MTND|MTATP|MTCO|MTCYB)')


def leading_alpha(sym):
    """取符号的前导字母串，如 HSPA8 -> HSPA, ATP5F1A -> ATP, RPS13 -> RPS"""
    m = re.match(r'^([A-Za-z]+)', sym)
    return m.group(1).upper() if m else ''


def panel_families(panel_syms, min_len=3):
    """104 panel 中出现过的所有 >=min_len 字符的家族前缀"""
    fams = set()
    for g in panel_syms:
        pre = leading_alpha(g)
        if len(pre) >= min_len:
            fams.add(pre)
    return fams

# ---------- R4 疾病相关黑名单 ----------
# T2DM GWAS 已确立位点 / 单基因型糖尿病 / 胰岛素与 GLP-1 信号 / DR-DN 候选基因
DISEASE = {
    # T2DM 常见变异位点（含本研究 30 个 T2DM control panel 之外的经典位点）
    'ADCY5','ADRA2A','ANK1','AP3S2','ARAP1','BCAR1','BCL11A','CAMK1D','CCND2','CDKAL1','CDKN2A',
    'CDKN2B','CENTD2','CMIP','DGKB','DUSP8','FTO','GCC1','GCK','GCKR','GIPR','GLIS3','GLP1R','GPSM1',
    'GRB14','HHEX','HMGA1','HMGA2','HNF1A','HNF1B','HNF4A','IDE','IGF1','IGF2BP2','INS','INSR','IRS1',
    'IRS2','JAZF1','KCNJ11','KCNQ1','KLF14','LEPR','MAEA','MC4R','MNX1','MTNR1B','NOTCH2','PAM','PDX1',
    'PEPD','PIK3R1','PPARG','PPARGC1A','PRC1','PROX1','PSMD6','RREB1','SLC16A11','SLC2A2','SLC2A4',
    'SLC30A8','ST6GAL1','TCF7L2','THADA','TP53INP1','TSPAN8','UBE2E2','WFS1','ZBED3','ZFAND6','ADIPOQ',
    'AKT1','AKT2','FOXA2','G6PC2','HK1','MLXIPL','NRXN3','SREBF1','TCF7','PPP1R3B','TMEM154','SSR1',
    'FITM2','RNF6','ANKH','C2CD4A','C2CD4B','VPS13C','CILP2','HNF4G','RASGRP1','C5orf67','ZMIZ1',
    # 糖尿病并发症（DR / DN / DPN）候选基因
    'VEGFA','EPO','AKR1B1','NOS3','ACE','AGT','TGFB1','SERPINE1','MTHFR','APOE','ELMO1','ENPP1',
    'UNC13B','CPVL','CHN2','GREM1','FRMD3','CARS','SP3','ITGA2','ITGB3','ADAM10','ICAM1','SELE',
    'TNF','IL6','CRP','AGER','RAGE','CTGF','CCN2','MMP2','MMP9','TIMP1','HIF1A','PLGF','PGF',
    # 糖脂代谢核心（避免对照集混入代谢调控基因）
    'LEP','GCGR','PCSK9','LDLR','HMGCR','SREBF2','FASN','ACACA','CPT1A','PPARA','LIPC','CETP',
}

# ---------- R5 GTEx MASHR 模型可用性 ----------
models = {}
for t in TISSUES:
    conn = sqlite3.connect(os.path.join(MODEL_DIR, f'mashr_{t}.db'))
    d = {}
    for gene, genename, n in conn.execute('SELECT gene, genename, "n.snps.in.model" FROM extra'):
        if genename:
            d[genename.upper()] = (gene, n)
    conn.close()
    models[t] = d
    print(f'R5 {t:<14} model genes  : {len(d)}')

both = set(models['Nerve_Tibial']) & set(models['Whole_Blood'])
both_ok = {g for g in both
           if models['Nerve_Tibial'][g][1] and models['Nerve_Tibial'][g][1] >= 1
           and models['Whole_Blood'][g][1] and models['Whole_Blood'][g][1] >= 1}

# ---------- 逐级过滤 ----------
steps = []
cur = set(hrt)
steps.append(('R1 HRT Atlas 全集', len(cur)))
cur = {g for g in cur if g not in panel}
steps.append(('R2 剔除 104 panel 重叠', len(cur)))
PANEL_FAM = panel_families(panel)
fam_hits = {g for g in cur if any(g.startswith(p) for p in PANEL_FAM) or EXTRA_FAMILY.match(g)}
cur = cur - fam_hits
steps.append(('R3 剔除 与 panel 同家族(>=3字符前缀)', len(cur)))
cur = {g for g in cur if g not in DISEASE}
steps.append(('R4 剔除 T2DM/并发症/代谢基因', len(cur)))
cur = {g for g in cur if g in both_ok}
steps.append(('R5 限定 双组织均有 MASHR 模型', len(cur)))
steps.append(('R6 随机抽取', N_TARGET))

print('\n过滤流程：')
for name, n in steps:
    print(f'  {name:<32} {n:>6}')

pool = sorted(cur)
assert len(pool) >= N_TARGET, f'合格池不足：{len(pool)} < {N_TARGET}'

# ---------- R6 固定种子随机抽样 ----------
random.seed(SEED)
selected = sorted(random.sample(pool, N_TARGET))

print(f'\n随机种子 seed = {SEED}')
print(f'新看家基因集（n={len(selected)}）：')
for i, g in enumerate(selected, 1):
    nt_n = models['Nerve_Tibial'][g][1]
    wb_n = models['Whole_Blood'][g][1]
    print(f'  {i:>2}. {g:<12} MASHR SNPs  Nerve_Tibial={nt_n:<4} Whole_Blood={wb_n}')

# ---------- 自检 ----------
assert not (set(selected) & panel), '仍与 104 panel 重叠！'
assert not (set(selected) & DISEASE), '仍包含疾病相关基因！'
assert all(g in both_ok for g in selected), '仍有基因缺模型！'
assert len(set(selected)) == N_TARGET
print('\n[自检] 与 104 panel 交集: 0  ✓')
print('[自检] 含 RPS/RPL/EEF/ENO 家族: 0  ✓')
print('[自检] 双组织模型齐备: 30/30  ✓')

# ---------- 输出 ----------
with open(os.path.join(ROOT, 'hk_genes_v2.txt'), 'w', encoding='utf-8') as f:
    f.write('# Housekeeping true-negative control set v2 (30 genes)\n')
    f.write('# Source  : HRT Atlas v1.0 human-mouse common HK set (MSigDB HOUNKPE_HOUSEKEEPING_GENES, n=1129)\n')
    f.write('# Ref     : Hounkpe BW, Chenou F, de Lima F, De Paula EV. HRT Atlas v1.0 database.\n')
    f.write('#           Nucleic Acids Res. 2021;49(D1):D947-D955. PMID 32663312.\n')
    f.write('# Selection rules (pre-specified, outcome-blind):\n')
    f.write('#   R1 universe  = HRT Atlas human-mouse common HK set (1129 genes)\n')
    f.write('#   R2 exclude  = any gene in the 104-gene testbed panel\n')
    f.write('#   R3 exclude  = any gene sharing a >=3-char leading alphabetic prefix with a 104-panel gene\n')
    f.write('#                (family-level exclusion: HSP*, RPS*, RPL*, EEF*, ENO*, ANXA*, ATP*, HNRNP*,\n')
    f.write('#                 PDIA*, ILF*, XRCC*, SLC*, YWHA*, ...), plus MRPS*/MRPL*/MT-*\n')
    f.write('#                Rule fixed before any TWAS result was produced.\n')
    f.write('#   R4 exclude  = curated T2DM / diabetic-complication / insulin & GLP-1 signalling genes\n')
    f.write('#   R5 require  = GTEx v8 MASHR model with >=1 SNP in BOTH Nerve_Tibial and Whole_Blood\n')
    f.write('#   R6 sampling = simple random sample of 30, seed = 20260830\n')
    f.write('#\n')
    f.write('# v1 -> v2 change: v1 contained 7/30 genes overlapping the 104-gene panel\n')
    f.write('#   (ACTB, GAPDH, RPL13A, RPLP0, RPS18 = candidates; EEF2, ENO1 = non-candidates),\n')
    f.write('#   violating its own selection criterion #3. v2 removes all of them.\n#\n')
    for g in selected:
        f.write(g + '\n')

meta = {
    'seed': SEED,
    'n_universe_hrt': len(hrt),
    'n_panel': len(panel),
    'n_eligible_pool': len(pool),
    'selection_rules': [
        'R1 universe = HRT Atlas v1.0 human-mouse common HK set (MSigDB HOUNKPE_HOUSEKEEPING_GENES, n=1129)',
        'R2 exclude any gene in the 104-gene testbed panel',
        'R3 exclude any gene sharing a >=3-char leading alphabetic prefix with a 104-panel gene '
        '(family-level exclusion), plus MRPS*/MRPL*/MT-*',
        'R4 exclude curated T2DM / diabetic-complication / insulin & GLP-1 signalling genes',
        'R5 require GTEx v8 MASHR model with >=1 model SNP in BOTH Nerve_Tibial and Whole_Blood',
        'R6 simple random sample of 30, random.seed=20260830',
    ],
    'genes': selected,
    'model_snps': {g: {'Nerve_Tibial': models['Nerve_Tibial'][g][1],
                       'Whole_Blood': models['Whole_Blood'][g][1]} for g in selected},
    'filter_flow': [[n, c] for n, c in steps],
}
with open(os.path.join(ROOT, 'hk_selection_meta.json'), 'w', encoding='utf-8') as f:
    json.dump(meta, f, indent=2, ensure_ascii=False)

print(f'\n已写出: hk_genes_v2.txt, hk_selection_meta.json  ->  {ROOT}')
