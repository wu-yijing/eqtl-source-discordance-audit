"""
Generate Figure 7: TWAS eQTL Source Selection Decision Framework
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
ax.set_xlim(0, 14)
ax.set_ylim(0, 10)
ax.axis('off')

# Colors
C_HEADER = '#2C3E50'
C_START = '#2980B9'
C_YES = '#27AE60'
C_YES_LIGHT = '#D5F5E3'
C_NO = '#E74C3C'
C_NO_LIGHT = '#FADBD8'
C_BODY = '#F5F5F5'
C_BORDER = '#BDC3C7'

def draw_box(ax, x, y, w, h, text, fc='white', ec=C_BORDER, fontsize=15, ha='center', va='center', bold=False, color='black', linewidth=1.5):
    """Draw a rounded rectangle with centered text."""
    box = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15", 
                                    fc=fc, ec=ec, linewidth=linewidth)
    ax.add_patch(box)
    # Handle multi-line text
    lines = text.split('\n')
    total_lines = len(lines)
    line_height = fontsize * 0.4
    start_y = y + h/2 + (total_lines-1)*line_height/2
    for i, line in enumerate(lines):
        ax.text(x + w/2, start_y - i*line_height, line, ha=ha, va='center',
                fontsize=fontsize, fontweight='bold' if bold else 'normal', color=color)
    return box

def draw_arrow(ax, x1, y1, x2, y2, label='', color='#7F8C8D'):
    """Draw an arrow between two points with optional label."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=2), 
                xycoords='data', textcoords='data')
    if label:
        mx, my = (x1+x2)/2, (y1+y2)/2
        ax.text(mx + 0.15, my, label, fontsize=15, color=color, fontweight='bold')

def draw_label(ax, x, y, text, color='#7F8C8D', fontsize=15):
    ax.text(x, y, text, fontsize=fontsize, color=color, ha='center', va='center')

# ============= TITLE =============
ax.text(7, 9.6, 'TWAS eQTL Source Selection Decision Framework', 
        fontsize=16, fontweight='bold', ha='center', color=C_HEADER)

# ============= START =============
draw_box(ax, 5.5, 8.6, 3, 0.7, 'TWAS Study', 
         fc=C_START, ec='#1A5276', color='white', fontsize=16, bold=True)
draw_arrow(ax, 7, 8.6, 7, 7.8)

# ============= DECISION 1: Multiple eQTL sources? =============
draw_box(ax, 4, 7.0, 6, 0.7, 'Do you have access to\nmultiple independent eQTL sources?', 
         fc='#F9E79F', ec='#F1C40F', fontsize=15, bold=True)

# Yes branch (right)
draw_arrow(ax, 8.5, 7.0, 10.5, 6.2, 'Yes', C_YES)
# No branch (left)
draw_arrow(ax, 5.5, 7.0, 3.5, 6.2, 'No', C_NO)

# ============= YES BRANCH (right side) =============
draw_box(ax, 8.5, 5.4, 4.0, 0.7, 'Dual-source sensitivity analysis', 
         fc=C_YES_LIGHT, ec=C_YES, color=C_YES, fontsize=15, bold=True)

# Direction consistency check
draw_box(ax, 8.5, 4.2, 4.0, 0.7, 'Compute direction consistency\nand Spearman \u03c1 between sources', 
         fc='white', ec=C_BORDER, fontsize=15)

draw_arrow(ax, 10.5, 4.2, 10.5, 3.4)

# Three threshold branches
# > 80%
draw_box(ax, 6.5, 2.6, 3.0, 0.6, '\u03c1 > 80%\nRobust', 
         fc=C_YES_LIGHT, ec=C_YES, fontsize=15, bold=True, color=C_YES)
draw_arrow(ax, 9.5, 3.4, 8.0, 3.2, '>80%', C_YES)

# 60-80%
draw_box(ax, 10.0, 1.8, 3.0, 0.6, '\u03c1 = 60\u201380%\nCaveat needed', 
         fc='#FEF9E7', ec='#F39C12', fontsize=15, bold=True, color='#D35400')
draw_arrow(ax, 10.5, 3.4, 11.5, 2.4, '60-80%', '#F39C12')

# < 60%
draw_box(ax, 6.5, 1.0, 3.0, 0.6, '\u03c1 < 60%\nSource-dependent', 
         fc=C_NO_LIGHT, ec=C_NO, fontsize=15, bold=True, color=C_NO)
draw_arrow(ax, 9.5, 3.4, 8.0, 1.6, '<60%', C_NO)

# ============= NO BRANCH (left side) =============
draw_box(ax, 0.5, 5.4, 4.5, 0.7, 'Single-source TWAS only\n\u2014 mandatory caveat', 
         fc=C_NO_LIGHT, ec=C_NO, color=C_NO, fontsize=15, bold=True)

# Source-specific limitations
draw_box(ax, 0.5, 3.8, 4.5, 0.7, 'If GTEx: Note small sample\n(N\u2248600) and tissue specificity', 
         fc='white', ec='#7FB3D8', fontsize=15)
draw_arrow(ax, 2.75, 5.4, 2.75, 4.5)

draw_box(ax, 0.5, 2.6, 4.5, 0.7, 'If eQTLGen: Note whole-\nblood-only limitation', 
         fc='white', ec='#7FB3D8', fontsize=15)
draw_arrow(ax, 2.75, 3.8, 2.75, 3.3)

# ============= BOTTOM - Post-hoc validation =============
# Connect both branches to bottom
draw_label(ax, 7, 0.5, '\u2193', '#7F8C8D', 16)

draw_box(ax, 3.0, -0.5, 8.0, 0.8, 'Post-hoc: Mahalanobis matching \u2502\nCross-population replication \u2502\nFunctional annotation stratification', 
         fc='#EBF5FB', ec='#2E86C1', fontsize=15, bold=False, color='#1A5276')

# ============= FOOTER =============
draw_label(ax, 7, -1.0, 
    'Recommended reporting:\n(1) Report \u03c1 for each source comparison\n(2) Flag source-dependent findings in Results\n(3) Use source-aware language in Discussion', 
    '#7F8C8D', 15)

# Save
plt.tight_layout()
plt.savefig('figs/Fig7_Decision_Framework.png', 
            dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('figs/Fig7_Decision_Framework.pdf', 
            bbox_inches='tight', facecolor='white')
plt.close()

print("Figure 7 saved: Fig7_Decision_Framework.png and Fig7_Decision_Framework.pdf")
