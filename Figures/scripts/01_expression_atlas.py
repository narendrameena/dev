#!/usr/bin/env python
"""
Figure 1 -- Q1 "is the gene expressed in kidney at all?"
            Q2 "which cell types, and is it specific or ubiquitous?"

RNA modalities only (SC_RNA + SN_RNA). SN_ATAC is excluded by design: those rows
are gene-activity scores, not transcript counts, and pooling them onto an
expression axis would be a category error.

Cell-type means are masked below MIN_CELLS so no point averages a handful of
cells. Nature Genetics layout: 183 mm double column, lettered panels, 5-7 pt type.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).parent))
import kidneyviz as kv  # noqa: E402

MIN_CELLS = 25
DETECT_MIN = 0.01


def load():
    z = np.load(kv.DERIVED / "panel_matrix.npz", allow_pickle=True)
    obs = pd.read_parquet(kv.DERIVED / "panel_obs.parquet")
    return z["lognorm"], z["counts"], list(z["genes"]), obs


def main():
    kv.use_style()
    X, C, genes, obs = load()
    rna = obs["tech"].isin(kv.RNA_TECHS).values
    ct = obs["celltype"].values
    print(f"RNA cells: {rna.sum():,} / {len(obs):,}")

    present = [c for c in kv.CT_ORDER if (rna & (ct == c)).sum() >= MIN_CELLS]

    # ---- per-cell-type mean matrix (genes x celltypes) ---------------------
    M = np.full((len(genes), len(present)), np.nan)
    N = np.zeros(len(present), int)
    for j, c in enumerate(present):
        m = rna & (ct == c)
        N[j] = m.sum()
        M[:, j] = X[m].mean(axis=0)
    pd.DataFrame(M, index=genes, columns=present).to_csv(
        kv.SOURCE / "Q2_celltype_mean_expression.csv")
    pd.DataFrame({"celltype": present, "n_RNA_cells": N}).to_csv(
        kv.SOURCE / "Q2_celltype_n.csv", index=False)

    detect = np.array([(C[rna, i] > 0).mean() for i in range(len(genes))])
    taus = np.array([kv.tau(M[i]) for i in range(len(genes))])
    top_ct = [present[int(np.nanargmax(M[i]))] for i in range(len(genes))]
    peak = np.nanmax(M, axis=1)

    q1 = pd.DataFrame({"gene": genes, "detect_RNA": detect, "tau": taus,
                       "peak_mean": peak, "top_celltype": top_ct})
    q1.to_csv(kv.SOURCE / "Q1_detection_per_gene.csv", index=False)

    # ======================================================= FIGURE 1 =======
    fig = plt.figure(figsize=(kv.W2, 10.2))
    gs = GridSpec(1, 2, width_ratios=[1.0, 1.05], wspace=0.26,
                  left=0.075, right=0.965, top=0.965, bottom=0.605)

    # ---- a  detection ------------------------------------------------------
    axa = fig.add_subplot(gs[0, 0])
    d = q1.sort_values("detect_RNA")
    y = np.arange(len(d))
    cols = [kv.STATUS["critical"] if v < DETECT_MIN else kv.SERIES[0]
            for v in d["detect_RNA"]]
    axa.barh(y, d["detect_RNA"] * 100, color=cols, height=0.78, lw=0)
    axa.set_yticks(y); axa.set_yticklabels(d["gene"].values, fontsize=4.2)
    axa.set_xlabel("cells with ≥1 count (%)")
    axa.axvline(DETECT_MIN * 100, color=kv.INK, lw=0.5, ls=(0, (3, 2)), zorder=4)
    axa.set_ylim(-1, len(d)); axa.set_xscale("symlog", linthresh=1)
    axa.set_xticks([0, 1, 10, 50]); axa.set_xticklabels(["0", "1", "10", "50"])
    kv.despine(axa); kv.hairline_grid(axa, "x"); kv.panel(axa, "a", dx=-0.12)
    axa.set_title("Detection across 282,610 scRNA + snRNA cells", loc="left",
                  pad=4)
    axa.text(1.15, len(d) * 0.02, "1% floor", fontsize=kv.FS_NOTE, color=kv.INK)

    # ---- b  abundance vs specificity --------------------------------------
    # Creative panel: separates "ubiquitous and abundant" from "rare but
    # cell-type restricted" -- the two ways a gene can look unimportant in a
    # bulk detection ranking while being biologically decisive.
    axb = fig.add_subplot(gs[0, 1])
    fam = {c: k for k, (n, mem) in enumerate(kv.LINEAGE.items()) for c in mem}
    cvec = [kv.SERIES[fam.get(t, 7) % len(kv.SERIES)] for t in q1["top_celltype"]]
    axb.scatter(q1["detect_RNA"] * 100, q1["tau"], s=17, c=cvec,
                linewidths=0.35, edgecolors=kv.SURFACE, zorder=3)
    axb.set_xscale("symlog", linthresh=1)
    axb.set_xticks([0, 1, 10, 50]); axb.set_xticklabels(["0", "1", "10", "50"])
    axb.set_xlim(-0.06, 95)
    axb.set_xlabel("cells with ≥1 count (%)")
    axb.set_ylabel("specificity index τ")
    axb.axvline(DETECT_MIN * 100, color=kv.INK, lw=0.5, ls=(0, (3, 2)), zorder=1)
    axb.axhline(0.85, color=kv.MUTED, lw=0.5, ls=":", zorder=1)
    # Hand-placed offsets: the automatic ones collided at the top-left cluster
    # (CRB2/IL1RL1) and ran off the right edge (PAX8/FOS).
    lab = {"CRB2": (4, 3), "IL1RL1": (5, -7), "HAVCR1": (5, 2), "UMOD": (5, 2),
           "CUBN": (5, -1), "HNF4A": (5, 2), "VCAM1": (-20, 3),
           "ACTA2": (5, 2), "COL1A1": (-24, -1), "FOS": (-16, -7),
           "VIM": (5, -6), "PAX8": (-20, 4)}
    for _, r in q1[q1.gene.isin(lab)].iterrows():
        axb.annotate(r.gene, (max(r.detect_RNA * 100, 0.05), r.tau),
                     textcoords="offset points", xytext=lab[r.gene],
                     fontsize=kv.FS_NOTE, color=kv.INK)
    kv.despine(axb); kv.hairline_grid(axb, "y"); kv.panel(axb, "b", dx=-0.115)
    axb.set_title("Rare ≠ unimportant: specificity vs abundance", loc="left",
                  pad=4)
    hs = [plt.Line2D([], [], marker="o", ls="", ms=3.4,
                     mfc=kv.SERIES[k % len(kv.SERIES)], mec="none", label=n)
          for k, n in enumerate(kv.LINEAGE)]
    axb.legend(handles=hs, loc="lower left", ncol=2, handletextpad=0.3,
               columnspacing=0.7, borderpad=0.2, labelspacing=0.25,
               title="peak cell-type lineage", title_fontsize=kv.FS_NOTE)

    # ---- c  the nephron atlas ---------------------------------------------
    # Dot plot, the single-cell convention (and the one the source paper uses):
    #   dot SIZE  = % of cells in that cell type with >=1 count  (detection)
    #   dot COLOUR= z-score of mean log expression across cell types (magnitude)
    # Two channels rather than one. Row-max scaling was rejected: it collapses a
    # gene averaging 0.05 and one averaging 5.0 onto the same 0-1 range, and it
    # cannot separate "a few cells express strongly" from "all express weakly".
    P = np.full((len(genes), len(present)), np.nan)     # percent expressing
    for j, c in enumerate(present):
        m = rna & (ct == c)
        P[:, j] = (C[m] > 0).mean(axis=0) * 100
    pd.DataFrame(P, index=genes, columns=present).to_csv(
        kv.SOURCE / "Q2_celltype_percent_expressing.csv")

    order = np.argsort(-taus)
    Mo, Po, go = M[order], P[order], [genes[i] for i in order]
    # z across cell types, per gene
    mu = np.nanmean(Mo, axis=1, keepdims=True)
    sd = np.nanstd(Mo, axis=1, keepdims=True)
    Z = np.divide(Mo - mu, sd, out=np.zeros_like(Mo), where=sd > 0)

    sub = GridSpec(2, 2, height_ratios=[0.032, 1], width_ratios=[1, 0.075],
                   hspace=0.02, wspace=0.012,
                   left=0.075, right=0.895, top=0.545, bottom=0.075)
    # Lay the columns out with a gutter between lineage blocks, so the nephron
    # grouping is carried by whitespace rather than by the colour bar alone.
    GAP = 1.3
    xpos, blocks, cur = {}, [], 0.0
    for name, mem in kv.LINEAGE.items():
        mem = [c for c in mem if c in present]
        if not mem:
            continue
        for c in mem:
            xpos[c] = cur
            cur += 1.0
        blocks.append((name, xpos[mem[0]], xpos[mem[-1]], len(mem)))
        cur += GAP
    XS = np.array([xpos[c] for c in present])
    XMAX = cur - GAP

    axc = fig.add_subplot(sub[1, 0])
    gx, gy = np.meshgrid(XS, np.arange(len(go)))
    zlim = float(np.nanpercentile(np.abs(Z), 99))
    im = axc.scatter(gx.ravel(), gy.ravel(),
                     s=(Po.ravel() / 100.0) * 11.0 + 0.12,
                     c=Z.ravel(), cmap=kv.cmap_div(), vmin=-zlim, vmax=zlim,
                     linewidths=0.12, edgecolors=kv.INK, zorder=3)
    axc.set_xlim(-0.8, XMAX - 0.2); axc.set_ylim(-0.6, len(go) - 0.4)
    axc.set_xticks(XS)
    axc.set_xticklabels(present, rotation=90, fontsize=4.6)
    axc.set_yticks(np.arange(len(go))); axc.set_yticklabels(go, fontsize=4.3)
    axc.set_yticks(np.arange(-.5, len(go)), minor=True)
    axc.grid(which="minor", axis="y", color=kv.GRID, lw=0.25, zorder=0)
    axc.set_axisbelow(True)
    axc.tick_params(which="minor", length=0)
    axc.tick_params(length=1.2)
    for s in axc.spines.values():
        s.set_visible(False)

    # slim tau panel, so the row ordering is legible rather than merely asserted
    axt = fig.add_subplot(sub[1, 1], sharey=axc)
    axt.barh(np.arange(len(go)), taus[order], color=kv.SERIES[1],
             height=0.78, lw=0)
    axt.set_xlim(0, 1.0); axt.set_xticks([0, 0.5, 1])
    axt.set_xlabel("τ", fontsize=kv.FS_LABEL, labelpad=1)
    axt.tick_params(labelleft=False, labelsize=kv.FS_NOTE, length=1.2)
    axt.axvline(0.85, color=kv.MUTED, lw=0.4, ls=":")
    kv.despine(axt, keep=("bottom",))

    # Lineage band. Labels sit INSIDE each block in reversed ink where the block
    # is wide enough, and are abbreviated otherwise -- placing them above the
    # band ran "Distal / collecting" into "Endothelium" on the narrow blocks.
    SHORT = {"Glomerulus": "Glom", "Proximal tubule": "Prox tubule",
             "Loop of Henle": "LoH", "Distal / collecting": "Distal/CD",
             "Endothelium": "Endo", "Stroma": "Stroma", "Immune": "Immune"}
    axband = fig.add_subplot(sub[0, 0], sharex=axc)
    for k, (name, x_lo, x_hi, n_mem) in enumerate(blocks):
        axband.add_patch(plt.Rectangle((x_lo - 0.5, 0), n_mem, 1, lw=0,
                                       color=kv.SERIES[k % len(kv.SERIES)]))
        if n_mem >= 2:
            axband.text((x_lo + x_hi) / 2, 0.5,
                        name if n_mem >= 5 else SHORT[name], ha="center",
                        va="center", fontsize=kv.FS_NOTE, color="white",
                        fontweight="bold", clip_on=True)
    axband.set_ylim(0, 1); axband.set_xlim(-0.8, XMAX - 0.2)
    axband.axis("off")
    kv.panel(axband, "c", dx=-0.078, dy=1.9)
    # Panel title in FIGURE coords -- in axes-data coords it was drawn far above
    # the band and ran straight through panel a's x-axis label.
    fig.text(0.075, 0.567, "Expression across the nephron, rows ordered by "
             "specificity (τ = 1 restricted to one cell type, τ = 0 ubiquitous)",
             ha="left", fontsize=kv.FS_TITLE, color=kv.INK)

    cax = fig.add_axes([0.945, 0.075, 0.008, 0.13])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("mean expression\n(z-score across cell types)",
                 fontsize=kv.FS_NOTE)
    cb.ax.tick_params(labelsize=kv.FS_NOTE, length=1.5)
    cb.outline.set_visible(False)

    # size legend -- without it the second channel is unreadable
    sax = fig.add_axes([0.935, 0.245, 0.055, 0.085]); sax.axis("off")
    for k, pct in enumerate([5, 25, 50, 100]):
        sax.scatter([0.16], [k * 0.26], s=(pct / 100) * 11.0 + 0.12,
                    c=kv.MUTED, linewidths=0.12, edgecolors=kv.INK)
        sax.text(0.42, k * 0.26, f"{pct}%", fontsize=kv.FS_NOTE,
                 va="center", color=kv.INK)
    sax.text(0, 1.10, "% of cells\nexpressing", fontsize=kv.FS_NOTE,
             va="bottom", color=kv.INK)
    sax.set_xlim(0, 1); sax.set_ylim(-0.2, 0.95)

    kv.save(fig, "Figure1_expression_atlas", caption=(
        "Figure 1 | Expression of the 56-gene panel in the human kidney. "
        "a, Fraction of scRNA+snRNA cells with at least one count per gene; "
        "red marks genes below a 1% detection floor. b, Specificity index tau "
        "against detection rate, coloured by the lineage of each gene's peak "
        "cell type; genes that are rare overall but cell-type restricted sit "
        "top-left. c, Mean log-normalised expression per cell type, each gene "
        "scaled to its own maximum, rows ordered by tau and columns ordered "
        "down the nephron. Cell types with <25 cells were excluded (n = "
        f"{len(present)} of 41 retained)."))

    print(f"  most restricted : {', '.join(np.array(genes)[order][:6])}")
    print(f"  most ubiquitous : {', '.join(np.array(genes)[order][-5:])}")
    print(f"  below 1% detect : {', '.join(q1[q1.detect_RNA < DETECT_MIN].gene)}")


if __name__ == "__main__":
    main()
