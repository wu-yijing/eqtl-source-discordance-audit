#!/usr/bin/env python3
"""Regenerate Figure 2 (scatter), Figure 4 (Love plot), Figure 5 (RNH1 attenuation)
after the hardcoded RNH1 GTEx value fix. Prints SMD for verification."""
import importlib.util, os, sys

spec = importlib.util.spec_from_file_location(
    "genfigs",
    os.path.join(os.path.dirname(__file__), "scripts", "python", "04_generate_all_figures.py"))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

data = m.load_data()
print("Loaded data. RNH1 GTEx dict:", data['rnh1_gtex'])

# Figure 2 (manuscript) = fig2_scatter -> Fig1_Scatter_Plot
m.fig2_scatter(data)
# Figure 4 (manuscript) = fig5_love_plot -> Fig4_Love_Plot
m.fig5_love_plot(data)
# Figure 5 (manuscript) = fig6_rnh1_attenuation -> Fig5_RNH1_Attenuation
m.fig6_rnh1_attenuation(data)

print("\nGenerated files in", m.FIG_DIR)
for f in sorted(os.listdir(m.FIG_DIR)):
    if f.startswith(('Fig1_Scatter', 'Fig4_Love', 'Fig5_RNH1')):
        print("  ", f)
