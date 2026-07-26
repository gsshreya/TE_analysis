#!/usr/bin/env python3
"""Canonical figure style for the TE manuscript 
"""
import matplotlib as mpl
import matplotlib.font_manager as fm

# ── Cell Press column widths  ──
MM = 1 / 25.4
COL1  = 85  * MM     # 1 column   (3.35 in)
COL15 = 114 * MM     # 1.5 column (4.49 in)
COL2  = 174 * MM     # 2 columns / full width (6.85 in)
MAXH  = 240 * MM     # max page height incl. legend

# ── Font sizes (pt, at final size) ──
FS_TICK, FS_LABEL, FS_LEGEND, FS_TITLE, FS_PANEL = 8, 9.5, 9, 11, 12  # larger, legible, uniform legends

# ── Canonical categorical palette: 8 linguistic groups ──
# Colour-blind-safe, print-legible (Okabe-Ito based; yellow darkened to a
# readable gold, distinct purple added). This REPLACES Figure 1's palette and
# is the one all figures should now use.
GROUP_PALETTE = {
    "AA_T":  "#E69F00",   # orange
    "CAO":   "#56B4E9",   # sky blue
    "DR_NT": "#009E73",   # green
    "DR_T":  "#C9A227",   # gold (legible replacement for pure yellow)
    "IE_NT": "#0072B2",   # blue
    "IE_T":  "#D55E00",   # vermillion
    "TB_NT": "#CC79A7",   # pink
    "TB_T":  "#7B3FA0",   # purple
}
GROUP_ORDER = ["AA_T", "CAO", "DR_NT", "DR_T", "IE_NT", "IE_T", "TB_NT", "TB_T"]

# ── TE-type palette (panels A/B) and 1KGP super-population palette (panel C) ──
TE_PALETTE  = {"Alu": "#4C78A8", "LINE1": "#D65F9E"}
POP_PALETTE = {"EUR": "#4C78A8", "AMR": "#F58518", "SAS": "#E45756",
               "EAS": "#72B7B2", "AFR": "#B279A2"}

def set_style():
    """Apply the shared rcParams. Call once at the top of each panel."""
    have_arial = any(f.name == "Arial" for f in fm.fontManager.ttflist)
    mpl.rcParams.update({
        "font.family": "Arial" if have_arial else "DejaVu Sans",
        "font.size": FS_LABEL,
        "axes.labelsize": FS_LABEL, "axes.titlesize": FS_TITLE,
        "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND, "legend.title_fontsize": FS_LEGEND,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 0.8, "patch.linewidth": 0.5,
        "pdf.fonttype": 42, "ps.fonttype": 42,          # editable text in Illustrator
        "savefig.dpi": 300, "savefig.bbox": "tight", "figure.dpi": 150,
    })

def group_label(ling, tribe):
    """Fig-1 grouping: no Tribe/Linguistic value -> CAO; else <Linguistic>_T / _NT."""
    import pandas as pd
    l = "" if (ling is None or (isinstance(ling, float) and pd.isna(ling))) else str(ling).strip()
    t = "" if (tribe is None or (isinstance(tribe, float) and pd.isna(tribe))) else str(tribe).strip()
    if l == "" or l == "CAO" or t == "" or t == "CAO":
        return "CAO"
    return f"{l}_T" if t == "Yes" else f"{l}_NT"

def save(fig, path_noext):
    for ext in ("pdf", "png"):
        fig.savefig(f"{path_noext}.{ext}")
    print(f"  saved {path_noext}.pdf / .png")
