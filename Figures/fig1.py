#!/usr/bin/env python3
"""
Figure 1 
"""
import os, gzip
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde
from PIL import Image, ImageDraw, ImageFont
import fig_style as S

# paths 
_FILT  = "/path/to/file"
_FIGIN = "/path/to/file"
_PCA   = _FIGIN + "/INS_PCA_common_only"
P = {
    # sample-level insertion VCFs (panels A histogram + D ridge)
    "ALU_INS":     _FILT + "/ALU.filtered.final.vcf.gz",
    "LINE1_INS":   _FILT + "/LINE1.filtered.final.vcf.gz",
    "POPMAP":      _FILT + "/sample_population_map.csv",
    # GI (tagged) + 1KGP for panel B overlap
    "ALU_TAGGED":  _FILT + "/ALU.filtered.final.tagged.vcf.gz",
    "L1_TAGGED":   _FILT + "/LINE1.filtered.final.tagged.vcf.gz",
    "KGP_ALU":     _FIGIN + "/3202_30X_1KGP_INS_ALU.vcf.gz",
    "KGP_L1":      _FIGIN + "/3202_30X_1KGP_INS_LINE1.vcf.gz",
    # common-insertion PCA (panel C) — same PCA as manuscript Fig 1D
    "EIGENVEC":    _PCA + "/ALU_LINE1.merged.pca.eigenvec",
    "EIGENVAL":    _PCA + "/ALU_LINE1.merged.pca.eigenval",
}

for k in P:
    P[k] = os.environ.get(f"FIG_{k}", P[k])
OUTDIR = os.environ.get("FIG_OUTDIR", "TE_ana/figs/results")
os.makedirs(OUTDIR, exist_ok=True)

TE_PAL  = S.TE_PALETTE                       # {"Alu","LINE1"}
POP_PAL = S.POP_PALETTE                      # {"EUR","AMR","SAS","EAS","AFR"}
GRP_PAL = dict(S.GROUP_PALETTE, All="#9a9a9a")
POPS    = ["EUR", "AMR", "SAS", "EAS", "AFR"]
INS_SLOP = 20
MISSING  = {"./.", ".|.", "."}


def count_per_sample(vcf):
    """Per-sample non-reference insertion count, keyed by base sample ID."""
    counts, order = {}, []
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                order = line.rstrip("\n").split("\t")[9:]
                counts = {s.split("_")[0]: 0 for s in order}
                continue
            if line.startswith("#"):
                continue
            for raw, g in zip(order, line.rstrip("\n").split("\t")[9:]):
                gt = g.split(":")[0]
                if gt in MISSING:
                    continue
                if "1" in gt:
                    s = raw.split("_")[0]
                    counts[s] += 1
    return counts

def load_sample_map(path):
    m = pd.read_csv(path)
    grp = {}
    for _, r in m.iterrows():
        grp[str(r["FID"]).split("_")[0]] = S.group_label(r.get("Linguistic_Group"), r.get("Tribe"))
    g2s = defaultdict(list)
    for s, g in grp.items():
        g2s[g].append(s)
    return grp, g2s

def _save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUTDIR, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}.pdf / .png")

# histogram 
def panel_A():
    fig, axes = plt.subplots(1, 2, figsize=(S.COL2 / 2, S.COL2 / 2 * 0.68))
    for ax, label, vcf in zip(axes, ("Alu", "LINE1"), (P["ALU_INS"], P["LINE1_INS"])):
        arr = np.array(list(count_per_sample(vcf).values()))
        ax.hist(arr, bins=80, color=TE_PAL[label], edgecolor="black", linewidth=0.2)
        lo, hi = np.percentile(arr, [0.5, 99.5])
        m = (hi - lo) * 0.05
        ax.set_xlim(lo - m, hi + m)
        ax.set_title(label, fontweight="bold")
        ax.set_xlabel("Insertions per individual")
    axes[0].set_ylabel("# of individuals")
    fig.tight_layout()
    _save(fig, "fig1_A")

#  1KGP overlap 
def _load_gi(vcf):
    import re
    out = []
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fl = line.rstrip().split("\t")
            ch = fl[0] if fl[0].startswith("chr") else "chr" + fl[0]
            out.append((ch, int(fl[1])))
    return out

def _load_kgp(vcf):
    import re
    tags = {"EUR": "EUR", "AMR": "AMR", "ASN": "EAS", "SAN": "SAS", "AFR": "AFR"}
    sites = {}
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fl = line.rstrip().split("\t")
            ch = fl[0] if fl[0].startswith("chr") else "chr" + fl[0]
            info = fl[7]
            pops = []
            for tag, lab in tags.items():
                m = re.search(fr"(?<![A-Z]){tag}_AF=([0-9.eE+\-]+)", info)
                if m and float(m.group(1)) > 0:
                    pops.append(lab)
            if pops:
                sites[(ch, int(fl[1]))] = pops
    return sites

def _overlap(kgp, gi):
    import bisect
    by = defaultdict(list)
    for ch, pos in gi:
        by[ch].append(pos)
    for ch in by:
        by[ch].sort()
    tot, ov = defaultdict(int), defaultdict(int)
    for (ch, pos), pops in kgp.items():
        starts = by.get(ch, [])
        hit = bisect.bisect_right(starts, pos + INS_SLOP) > bisect.bisect_left(starts, pos - INS_SLOP)
        for p in pops:
            tot[p] += 1
            if hit:
                ov[p] += 1
    return tot, ov

def panel_B():
    groups = [("Alu", _overlap(_load_kgp(P["KGP_ALU"]), _load_gi(P["ALU_TAGGED"]))),
              ("LINE1", _overlap(_load_kgp(P["KGP_L1"]), _load_gi(P["L1_TAGGED"])))]
    fig, ax = plt.subplots(figsize=(S.COL2 / 2, S.COL2 / 2 * 0.80))
    bw, inner, gap = 0.08, 0.012, 1.0
    offs = [-(len(POPS) * bw + (len(POPS) - 1) * inner) / 2 + i * (bw + inner) + bw / 2 for i in range(len(POPS))]
    for i, (lab, (tot, ov)) in enumerate(groups):
        c = i * gap
        ax.text(c, 103, lab, ha="center", va="bottom", fontweight="bold")
        for pop, off in zip(POPS, offs):
            t, o = tot.get(pop, 0), ov.get(pop, 0)
            op = (o / t * 100) if t else 0
            x = c + off
            ax.bar(x, op, width=bw, color=POP_PAL[pop])
            ax.bar(x, 100 - op, width=bw, bottom=op, color=POP_PAL[pop], alpha=0.30, linewidth=0)
            ax.text(x, -2, pop, ha="center", va="top", rotation=90)
    ax.set_xticks([]); ax.set_ylim(-14, 108)
    ax.set_ylabel("% 1KGP variants overlapping GI")
    for sp in ("top", "right", "bottom"): ax.spines[sp].set_visible(False)
    ax.yaxis.grid(True, linestyle="--", alpha=0.4); ax.set_axisbelow(True)
    fig.legend(handles=[mpatches.Patch(facecolor="grey", label="Overlaps GI callset"),
                        mpatches.Patch(facecolor="grey", alpha=0.30, label="Not in GI callset")],
               loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    _save(fig, "fig1_B")

# PCA 
def panel_C():
    ev = pd.read_csv(P["EIGENVEC"], sep=r"\s+", usecols=[0, 1, 2, 3, 4, 5])
    ev.columns = [c.lstrip("#") for c in ev.columns]
    ev[]
    df = ev.merge(pd.read_csv(P["POPMAP"]), left_on="IID", right_on="FID", how="inner")
    df["group"] = [S.group_label(l, t) for l, t in zip(df["Linguistic_Group"], df["Tribe"])]
    pct = pd.read_csv(P["EIGENVAL"], header=None)[0].to_numpy(); pct = pct / pct.sum() * 100
    fig, axes = plt.subplots(1, 2, figsize=(S.COL2, S.COL2 * 0.46))
    for ax, (x, y) in zip(axes, [("PC1", "PC2"), ("PC3", "PC4")]):
        for g in S.GROUP_ORDER:
            sub = df[df["group"] == g]
            ax.scatter(sub[x], sub[y], s=4, alpha=0.55, lw=0, color=GRP_PAL[g], rasterized=True)
        xi, yi = int(x[2:]) - 1, int(y[2:]) - 1
        ax.set_xlabel(f"{x} ({pct[xi]:.2f}%)"); ax.set_ylabel(f"{y} ({pct[yi]:.2f}%)")
        ax.set_title(f"{x} vs {y}", fontweight="bold")
        ax.axhline(0, color="0.6", lw=0.4, ls="--"); ax.axvline(0, color="0.6", lw=0.4, ls="--")
    handles = [Line2D([0], [0], marker="o", ls="", ms=6, mfc=GRP_PAL[g], mec="none", label=g) for g in S.GROUP_ORDER]
    fig.legend(handles=handles, title="Linguistic group", ncol=8, loc="lower center",
               bbox_to_anchor=(0.5, -0.04), frameon=False)
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    _save(fig, "fig1_C")

# ridge (variants per sample)
def panel_D():
    grp, g2s = load_sample_map(P["POPMAP"])
    all_groups = sorted(g2s)
    row_order = ["All"] + [g for g in all_groups if g != "CAO"] + (["CAO"] if "CAO" in all_groups else [])
    data = {}
    for lab, vcf in (("Alu", P["ALU_INS"]), ("LINE1", P["LINE1_INS"])):
        c = count_per_sample(vcf)
        d = {"All": np.array(list(c.values()))}
        for g in all_groups:
            vals = [c[s] for s in g2s[g] if s in c]
            if len(vals) > 1:
                d[g] = np.array(vals)
        data[lab] = d
    ROW_H, RIDGE = 1.8, 1.8 * 0.85
    fig, axes = plt.subplots(1, 2, figsize=(S.COL2, S.COL2 * 0.42))
    for ax, lab in zip(axes, ("Alu", "LINE1")):
        vt = data[lab]
        allv = np.concatenate([v for v in vt.values() if len(v) > 1])
        pad = (allv.max() - allv.min()) * 0.08
        grid = np.linspace(allv.min() - pad, allv.max() + pad, 500)
        n = len(row_order)
        for ri, g in enumerate(row_order):
            if g not in vt or len(vt[g]) < 2:
                continue
            yb = (n - 1 - ri) * ROW_H
            dens = gaussian_kde(vt[g], bw_method="scott")(grid)
            yv = yb + dens / dens.max() * RIDGE
            ax.fill_between(grid, yb, yv, alpha=0.82, color=GRP_PAL.get(g, "#999"), linewidth=0)
            ax.plot(grid, yv, color="black", linewidth=0.8)
        ax.axvline(np.median(vt["All"]), color="black", lw=0.9, ls="--", alpha=0.6)
        ax.set_yticks([(n - 1 - i) * ROW_H for i in range(n)])
        ax.set_yticklabels(row_order if lab == "Alu" else [""] * n)
        ax.set_ylim(-ROW_H * 0.3, n * ROW_H)
        ax.set_xlabel("Variants per sample"); ax.set_title(lab, fontweight="bold")
        for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    fig.tight_layout()
    _save(fig, "fig1_D")


def main():
    S.set_style()
    for fn in (panel_A, panel_B, panel_C, panel_D):
        print(f"[run] {fn.__name__}")
        try:
            fn()
        except Exception as e:
            print(f"  [SKIP] {fn.__name__}: {type(e).__name__}: {e}")

if __name__ == "__main__":
    main()
