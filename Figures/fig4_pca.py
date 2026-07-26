#!/usr/bin/env python3
"""
Figure 3 PCA panels — CABYR & CYP2J2 carriers on the GI cohort PCA.

Outputs (into GI_OUTDIR, default = current dir):
    fig_pca_CABYR_PC1PC2.pdf/.png
    fig_pca_CYP2J2_PC1PC2.pdf/.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

#  inputs (HPC defaults) 
_PCA  = "/path/to/file"
_FILT = "/path/to/file"
EIGENVEC = os.environ.get("GI_EIGENVEC", _PCA + "/ALU_LINE1.merged.pca.eigenvec")
EIGENVAL = os.environ.get("GI_EIGENVAL", _PCA + "/ALU_LINE1.merged.pca.eigenval")
POPMAP   = os.environ.get("GI_POPMAP",   _FILT + "/sample_population_map.csv")
OUTDIR   = os.environ.get("GI_OUTDIR",   "/gpfs/data/user/vidushi/TE_ana/figs/results")
os.makedirs(OUTDIR, exist_ok=True)

# canonical style 
MM = 1 / 25.4
COL1, COL15, COL2 = 85 * MM, 114 * MM, 174 * MM      # Cell Press 1 / 1.5 / 2-col widths
FS_TICK, FS_LABEL, FS_LEGEND, FS_TITLE = 8, 9.5, 9, 11   # larger, legible, uniform legends

# colour-blind-safe qualitative palette for carrier sub-populations
CARRIER_POP_COLORS = ["#D55E00", "#009E73", "#0072B2", "#E69F00",
                      "#CC79A7", "#56B4E9", "#7B3FA0", "#000000"]

CARRIERS = {
    "CABYR": ["A", "C"],
    "CYP2J2": ["AI", "AU"],
}

def set_style():
    have_arial = any(f.name == "Arial" for f in fm.fontManager.ttflist)
    mpl.rcParams.update({
        "font.family": "Arial" if have_arial else "DejaVu Sans",
        "font.size": FS_LABEL, "axes.labelsize": FS_LABEL, "axes.titlesize": FS_TITLE,
        "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND, "legend.title_fontsize": FS_LEGEND,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300, "savefig.bbox": "tight",
    })

#  data 
def load():
    ev = pd.read_csv(EIGENVEC, sep=r"\s+", usecols=[0, 1, 2, 3, 4, 5])
    ev.columns = [c.lstrip("#") for c in ev.columns]            # FID IID PC1..PC4
    pop = pd.read_csv(POPMAP)                                    # FID,Code,Linguistic_Group,Tribe
    df = ev.merge(pop, left_on="IID", right_on="FID", how="inner").reset_index(drop=True)
    val = pd.read_csv(EIGENVAL, header=None)[0].to_numpy()
    pct = val / val.sum() * 100.0
    print(f"{len(df)} samples in PCA")
    print(f"PC1-4 variance: {pct[0]:.2f}, {pct[1]:.2f}, {pct[2]:.2f}, {pct[3]:.2f} %")
    for te, ids in CARRIERS.items():
        sub = df[df["IID"].isin(ids)]
        print(f"  {te}: {len(sub)}/{len(ids)} carriers found; populations: "
              f"{sorted(sub['Code'].dropna().unique())}")
    return df, pct

# plotting 
def scatter_pc(ax, df, xpc, ypc, pct, te):
    # grey background cloud = all samples (non-carriers dominate)
    ax.scatter(df[xpc], df[ypc], s=6, alpha=0.30, color="0.6",
               linewidths=0, rasterized=True, zorder=1)
    # carriers as stars, coloured by sub-population (Code)
    sub = df[df["IID"].isin(CARRIERS[te])]
    for i, code in enumerate(sorted(sub["Code"].dropna().unique())):
        s2 = sub[sub["Code"] == code]
        ax.scatter(s2[xpc], s2[ypc], s=90, marker="*",
                   color=CARRIER_POP_COLORS[i % len(CARRIER_POP_COLORS)],
                   edgecolors="black", linewidths=0.5, zorder=5,
                   label=f"Carrier — {code}")
    xi, yi = int(xpc[2:]) - 1, int(ypc[2:]) - 1
    ax.set_xlabel(f"{xpc} ({pct[xi]:.2f}%)")
    ax.set_ylabel(f"{ypc} ({pct[yi]:.2f}%)")
    ax.axhline(0, color="0.6", lw=0.4, ls="--", zorder=0)
    ax.axvline(0, color="0.6", lw=0.4, ls="--", zorder=0)
    ax.tick_params(length=2.5)

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUTDIR, f"{name}.{ext}"))
    print(f"  saved {name}.pdf / .png -> {OUTDIR}")

def main():
    set_style()
    df, pct = load()
    for te in CARRIERS:
        fig, ax = plt.subplots(figsize=(COL15, COL15 * 0.78))
        scatter_pc(ax, df, "PC1", "PC2", pct, te)
        ax.set_title(f"{te} — TE carriers by population", fontweight="bold")
        ax.legend(title="Population", loc="upper right", frameon=True,
                  framealpha=0.9, handletextpad=0.3, borderpad=0.4)
        save(fig, f"fig_pca_{te}_PC1PC2"); plt.close(fig)

if __name__ == "__main__":
    main()
