# -*- coding: utf-8 -*-
"""Generate Figure 8: 2x2 decomposition matrix of TWAS discordance.
Values are taken from decomposition_results_v2.json (primary) and the manuscript
v2 text. Produces SVG (vector, submission), PDF (submission) and PNG (600 dpi preview).
"""
import os
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.colors as mcolors
import matplotlib.cm as cm

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- lock numbers from the v2 decomposition (matches manuscript text) ----
rho_panel, cons_panel, n_panel = 0.31, 64.2, 162   # 1D-panel: GTEx WB vs eQTLGen WB
rho_tissue, cons_tissue, n_tissue = 0.45, 68.0, 150  # 1D-tissue: GTEx WB vs GTEx NT
rho_dual, cons_dual, n_dual = 0.26, 58.7, 201        # 2D: GTEx Multi vs eQTLGen
rho_anchor = 0.29                                    # primary GTEx-eQTLGen, n=102
rho_ref = 1.00

CELLS = {
    # (row, col)  row 0 = top (matched tissue), row 1 = bottom (mismatched tissue)
    (0, 0): dict(kind="ref", title="Reference",
                 sub="same source, same tissue",
                 lines=[r"$\rho$ = 1.00 (by construction)", "discordance-free"]),
    (0, 1): dict(kind="cmp", title="Panel-only mismatch",
                 sub="GTEx v8 whole-blood\nvs eQTLGen whole-blood",
                 rho=rho_panel, cons=cons_panel, n=n_panel),
    (1, 0): dict(kind="cmp", title="Tissue-only mismatch",
                 sub="GTEx v8 whole-blood\nvs GTEx v8 Nerve_Tibial",
                 rho=rho_tissue, cons=cons_tissue, n=n_tissue),
    (1, 1): dict(kind="cmp", title="Dual mismatch",
                 sub="GTEx v8 multi-tissue ACAT-O\nvs eQTLGen",
                 rho=rho_dual, cons=cons_dual, n=n_dual),
}

# ---- figure geometry (normalized 0..1 on a full-canvas axes) ----
CW, CH = 0.265, 0.255          # cell width / height
GAP_X, GAP_Y = 0.055, 0.045    # gap between columns / rows
X0 = 0.345                     # left edge of column 0
Y0_BOT = 0.175                 # bottom edge of bottom row
x_left = X0
x_right = X0 + CW + GAP_X
y_bot = Y0_BOT
y_top = Y0_BOT + CH + GAP_Y
cell_pos = {
    (0, 0): (x_left, y_top),
    (0, 1): (x_right, y_top),
    (1, 0): (x_left, y_bot),
    (1, 1): (x_right, y_bot),
}

cmap = matplotlib.colormaps["RdYlGn"]
norm = mcolors.Normalize(vmin=0.20, vmax=0.50)

fig = plt.figure(figsize=(7.4, 6.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

# ---------- title ----------
ax.text(0.5, 0.975, "2×2 decomposition of TWAS discordance",
        ha="center", va="top", fontsize=15, fontweight="bold", color="#1a1a1a")

# ---------- column axis (eQTL-source / panel) ----------
ax.text(0.5, 0.905, "eQTL-source / panel axis", ha="center", va="bottom",
        fontsize=11.5, fontweight="bold", color="#1a1a1a")
ax.text((x_left + CW/2), 0.875, "GTEx v8\n(same source)", ha="center", va="top",
        fontsize=9.5, color="#333333")
ax.text((x_right + CW/2), 0.875, "eQTLGen\n(different source)", ha="center", va="top",
        fontsize=9.5, color="#333333")
# bracket under column label
for xc in (x_left, x_right + CW):
    ax.plot([xc, xc], [0.845, 0.862], color="#555555", lw=1.1)
ax.plot([x_left, x_left + CW], [0.862, 0.862], color="#555555", lw=1.1)
ax.plot([x_right, x_right + CW], [0.862, 0.862], color="#555555", lw=1.1)

# ---------- row axis (tissue-context) ----------
ax.text(0.115, 0.49, "Tissue-context axis", rotation=90, ha="center", va="center",
        fontsize=11.5, fontweight="bold", color="#1a1a1a")
ax.text(0.215, (y_top + CH/2), "Matched\n(whole blood)", rotation=90,
        ha="center", va="center", fontsize=9.5, color="#333333")
ax.text(0.215, (y_bot + CH/2), "Mismatched\n(Nerve_Tibial / multi-tissue)", rotation=90,
        ha="center", va="center", fontsize=9.5, color="#333333")
for yc in (y_top, y_bot + CH):
    ax.plot([0.245, 0.262], [yc, yc], color="#555555", lw=1.1)
ax.plot([0.245, 0.245], [y_bot + CH, y_top], color="#555555", lw=1.1)

# ---------- cells ----------
def draw_cell(x0, y0, info):
    if info["kind"] == "ref":
        fc = "#e6e6e6"
    else:
        fc = cmap(norm(info["rho"]))
    box = FancyBboxPatch((x0, y0), CW, CH,
                         boxstyle="round,pad=0.004,rounding_size=0.012",
                         linewidth=1.3, edgecolor="#3a3a3a", facecolor=fc, zorder=2)
    ax.add_patch(box)
    cx = x0 + CW/2
    # title
    ax.text(cx, y0 + CH - 0.018, info["title"], ha="center", va="top",
            fontsize=10.5, fontweight="bold", color="#1a1a1a", zorder=3)
    # sub
    ax.text(cx, y0 + CH - 0.048, info["sub"], ha="center", va="top",
            fontsize=7.0, color="#444444", zorder=3)
    if info["kind"] == "cmp":
        ax.text(cx, y0 + CH*0.50, r"$\rho$ = %.2f" % info["rho"],
                ha="center", va="center", fontsize=15, fontweight="bold",
                color="#111111", zorder=3)
        ax.text(cx, y0 + CH*0.30, "%.1f%% direction consistency" % info["cons"],
                ha="center", va="center", fontsize=9.0, color="#222222", zorder=3)
        ax.text(cx, y0 + CH*0.155, "n = %d gene–phenotype pairs" % info["n"],
                ha="center", va="center", fontsize=8.4, color="#222222", zorder=3)
    else:
        ax.text(cx, y0 + CH*0.46, info["lines"][0], ha="center", va="center",
                fontsize=11, fontweight="bold", color="#333333", zorder=3)
        ax.text(cx, y0 + CH*0.27, info["lines"][1], ha="center", va="center",
                fontsize=9.5, color="#333333", zorder=3)

for (r, c), pos in cell_pos.items():
    draw_cell(pos[0], pos[1], CELLS[(r, c)])

# ---------- arrow: panel-only bounds eQTL-source axis discordance ----------
ax.annotate("", xy=(x_right + CW + 0.01, y_top + CH*0.5),
            xytext=(x_right + CW + 0.10, y_top + CH*0.5),
            arrowprops=dict(arrowstyle="-|>", color="#777777", lw=1.1))
ax.text(x_right + CW + 0.115, y_top + CH*0.5,
        "bounds discordance\nattributable to the\neQTL-source axis\n(primary ρ=%.2f, n=102)" % rho_anchor,
        ha="left", va="center", fontsize=7.6, color="#666666")

# ---------- colorbar ----------
import matplotlib.colorbar as cb
sm = cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cax = fig.add_axes([X0, 0.085, (CW*2 + GAP_X), 0.022])
cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
cbar.set_ticks([0.20, 0.30, 0.40, 0.50])
cbar.set_ticklabels(["0.20", "0.30", "0.40", "0.50"])
cbar.ax.tick_params(labelsize=8)
cbar.outline.set_visible(False)
ax.text(X0, 0.115, "Spearman ρ (higher = more concordant / less discordant)",
        ha="left", va="center", fontsize=8.4, color="#333333")

# ---------- footnote ----------
ax.text(0.5, 0.05,
        ("Off-diagonal cells report Spearman ρ and direction-consistency of the corresponding pairwise TWAS comparison.\n"
         "Reference (same source, same tissue) cell is discordance-free by construction.  "
         "Primary GTEx v8–eQTLGen comparison (n=102): ρ=%.2f.\n"
         "Robustness: ordering preserved in the common-gene universe (n=144) and all three phenotype strata (DR/DN/DPN)."
         % rho_anchor),
        ha="center", va="top", fontsize=7.3, color="#555555")

# ---------- save ----------
base = os.path.join(HERE, "Figure8_2x2_decomposition")
fig.savefig(base + ".svg", format="svg", dpi=300, bbox_inches="tight")
fig.savefig(base + ".pdf", format="pdf", dpi=300, bbox_inches="tight")
fig.savefig(base + ".png", format="png", dpi=600, bbox_inches="tight")
print("saved:", base + ".svg / .pdf / .png")

# dump a small json manifest of the plotted numbers for traceability
manifest = {
    "panel_only": {"rho": rho_panel, "consistency": cons_panel, "n": n_panel},
    "tissue_only": {"rho": rho_tissue, "consistency": cons_tissue, "n": n_tissue},
    "dual_mismatch": {"rho": rho_dual, "consistency": cons_dual, "n": n_dual},
    "anchor_primary": {"rho": rho_anchor, "n": 102},
    "reference": {"rho": rho_ref},
}
with open(base + "_manifest.json", "w", encoding="utf-8") as f:
    json.dump(manifest, f, indent=2, ensure_ascii=False)
print("manifest saved.")
