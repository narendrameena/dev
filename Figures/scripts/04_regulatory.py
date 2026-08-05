#!/usr/bin/env python
"""
Figure 4 -- Q6 "identify regulatory mechanisms and causal cell types"

The CUBN / HNF4A argument, run for the whole panel:

  1. Is the promoter open, and specifically in which cell type?
     Per-cell-type snATAC coverage (Abedini et al. 2024, susztaklab.com),
     quantified as signal over each track's own genome-wide background so that
     tracks of different depth are comparable.
  2. Does the open region carry the TF's motif?
     FIMO scan with HOCOMOCO v11 HNF4A_HUMAN.H11MO.0.A.
  3. Does the TF and its putative target move together across cell states?

IMPORTANT LIMITATIONS, stated because they bound what may be concluded:
  * These are coverage tracks, not MACS2 peak calls. "Enrichment over track
    background" is a quantitative proxy for accessibility, not a called peak.
  * A motif match is sequence evidence of a possible binding site. It is not
    binding. Confirming occupancy needs ChIP-seq or footprinting.
  * Only the CONTROL tracks are usable. The disease tracks on the portal are
    ~100x smaller and internally inconsistent (see README); no disease-vs-
    control chromatin contrast is attempted here.
"""
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).parent))
import kidneyviz as kv  # noqa: E402

import pyBigWig  # noqa: E402

PROM_UP, PROM_DN = 1000, 1000     # window around the TSS
CTRL = kv.ATAC_DIR / "control"


def tss_table():
    """TSS per panel gene from refGene, strand-aware, canonical chromosomes."""
    want = set(kv.PANEL)
    best = {}
    with gzip.open(kv.REFGENE, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            name2, chrom, strand = f[12], f[2], f[3]
            if name2 not in want or "_" in chrom:
                continue
            start, end = int(f[4]), int(f[5])
            tss = end if strand == "-" else start
            span = end - start
            # keep the longest transcript as the representative
            if name2 not in best or span > best[name2][3]:
                best[name2] = (chrom, strand, tss, span)
    rows = [dict(gene=g, chrom=c, strand=s, tss=t,
                 start=max(t - PROM_UP, 0), end=t + PROM_DN)
            for g, (c, s, t, _) in best.items()]
    return pd.DataFrame(rows).sort_values("gene")


def promoter_matrix(tss):
    """gene x celltype accessibility, normalised to each track's background."""
    tracks = sorted(CTRL.glob("*.bw"))
    cts = [t.stem for t in tracks]
    M = np.full((len(tss), len(cts)), np.nan)
    for j, t in enumerate(tracks):
        bw = pyBigWig.open(str(t))
        h = bw.header()
        bg = h["sumData"] / h["nBasesCovered"] if h["nBasesCovered"] else np.nan
        chroms = bw.chroms()
        for i, r in enumerate(tss.itertuples()):
            if r.chrom not in chroms:
                continue
            e = min(r.end, chroms[r.chrom])
            if r.start >= e:
                continue
            try:
                v = bw.stats(r.chrom, int(r.start), int(e), type="mean")[0]
            except Exception:
                v = None
            M[i, j] = (v or 0.0) / bg if bg and bg > 0 else np.nan
        bw.close()
    return pd.DataFrame(M, index=list(tss.gene), columns=cts)


def fimo_promoters(tss, motif="HNF4A_HUMAN.H11MO.0.A", thresh=1e-4):
    """Scan every panel promoter for the motif; returns hits per gene."""
    import urllib.request
    import json
    fa = []
    for r in tss.itertuples():
        url = ("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;"
               f"chrom={r.chrom};start={r.start};end={r.end}")
        try:
            with urllib.request.urlopen(url, timeout=60) as fh:
                seq = json.load(fh)["dna"].upper()
        except Exception as exc:
            print(f"    ! {r.gene}: {exc}")
            continue
        fa.append(f">{r.gene}\n" + "\n".join(seq[i:i+60]
                                             for i in range(0, len(seq), 60)))
    with tempfile.TemporaryDirectory() as td:
        fpath = Path(td) / "prom.fa"
        fpath.write_text("\n".join(fa) + "\n")
        out = Path(td) / "fimo"
        subprocess.run([kv.FIMO, "--verbosity", "1", "--thresh", str(thresh),
                        "--motif", motif, "--oc", str(out),
                        str(kv.MOTIF_DB), str(fpath)], check=True)
        hits = pd.read_csv(out / "fimo.tsv", sep="\t", comment="#")
    hits = hits.dropna(subset=["sequence_name"])
    return hits


def main():
    kv.use_style()
    tss = tss_table()
    print(f"TSS resolved for {len(tss)}/{len(kv.PANEL)} panel genes")
    missing = sorted(set(kv.PANEL) - set(tss.gene))
    if missing:
        print(f"  not in refGene: {missing}")
    tss.to_csv(kv.SOURCE / "Q6_panel_TSS_hg38.csv", index=False)

    print("quantifying promoter accessibility across 23 cell types...")
    P = promoter_matrix(tss)
    P.to_csv(kv.SOURCE / "Q6_promoter_accessibility.csv")

    print("scanning promoters for the HNF4A motif...")
    hits = fimo_promoters(tss)
    hits.to_csv(kv.SOURCE / "Q6_HNF4A_motif_hits.csv", index=False)
    nhit = hits.groupby("sequence_name").size()
    best = hits.groupby("sequence_name")["p-value"].min()
    print(f"  genes with >=1 HNF4A motif (p<1e-4): {len(nhit)}/{len(tss)}")

    # expression context: HNF4A and CUBN across PT states
    z = np.load(kv.DERIVED / "panel_matrix.npz", allow_pickle=True)
    obs = pd.read_parquet(kv.DERIVED / "panel_obs.parquet")
    X, G = z["lognorm"], list(z["genes"])
    rna = obs["tech"].isin(kv.RNA_TECHS).values
    ct = obs["celltype"].values
    states = ["PT_S1", "PT_S2", "PT_S3", "iPT"]
    expr = {g: [X[rna & (ct == s), G.index(g)].mean() for s in states]
            for g in ("HNF4A", "CUBN")}
    pd.DataFrame(expr, index=states).to_csv(kv.SOURCE / "Q6_PT_state_expression.csv")

    # ======================================================= FIGURE 4 =======
    PT = ["PT_S1", "PT_S2", "PT_S3", "iPT"]
    fig = plt.figure(figsize=(kv.W2, 8.0))
    gs = GridSpec(2, 3, height_ratios=[0.85, 1.3], width_ratios=[1.15, 1, 1],
                  hspace=0.44, wspace=0.34,
                  left=0.075, right=0.935, top=0.93, bottom=0.085)

    # ---- a  CUBN promoter accessibility across cell types ------------------
    axa = fig.add_subplot(gs[0, 0])
    cu = P.loc["CUBN"].sort_values(ascending=False)
    y = np.arange(len(cu))
    cols = [kv.SERIES[1] if c in PT else kv.MUTED for c in cu.index]
    axa.barh(y, cu.values, color=cols, height=0.76, lw=0)
    axa.set_yticks(y); axa.set_yticklabels(cu.index, fontsize=4.6)
    axa.invert_yaxis()
    axa.set_xlabel("CUBN promoter accessibility\n(signal / track background)")
    axa.axvline(1, color=kv.INK, lw=0.5, ls=(0, (3, 2)))
    kv.despine(axa); kv.hairline_grid(axa, "x"); kv.panel(axa, "a", dx=-0.30)
    axa.set_title("CUBN promoter is open in PT", loc="left", pad=4)
    axa.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, fc=kv.SERIES[1], lw=0, label="proximal tubule"),
        plt.Rectangle((0, 0), 1, 1, fc=kv.MUTED, lw=0, label="other cell types")],
        loc="lower right", handlelength=1.0, handleheight=0.7,
        borderpad=0.25, labelspacing=0.3)

    # ---- b  UMOD positive control -----------------------------------------
    axb = fig.add_subplot(gs[0, 1])
    um = P.loc["UMOD"].sort_values(ascending=False)
    TAL = ["C_TAL", "M_TAL"]
    y = np.arange(len(um))
    cols = [kv.SERIES[2] if c in TAL else kv.MUTED for c in um.index]
    axb.barh(y, um.values, color=cols, height=0.76, lw=0)
    axb.set_yticks(y); axb.set_yticklabels(um.index, fontsize=4.6)
    axb.invert_yaxis()
    axb.set_xlabel("UMOD promoter accessibility\n(signal / track background)")
    axb.axvline(1, color=kv.INK, lw=0.5, ls=(0, (3, 2)))
    kv.despine(axb); kv.hairline_grid(axb, "x"); kv.panel(axb, "b", dx=-0.30)
    axb.set_title("Positive control: UMOD in TAL", loc="left", pad=4)
    axb.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, fc=kv.SERIES[2], lw=0,
                      label="thick ascending limb"),
        plt.Rectangle((0, 0), 1, 1, fc=kv.MUTED, lw=0, label="other cell types")],
        loc="lower right", handlelength=1.0, handleheight=0.7,
        borderpad=0.25, labelspacing=0.3)

    # ---- c  HNF4A and CUBN track each other across PT states ---------------
    axc = fig.add_subplot(gs[0, 2])
    xx = np.arange(len(states))
    axc.plot(xx, expr["HNF4A"], "-o", color=kv.SERIES[0], ms=3.2, lw=1.0,
             label="HNF4A (TF)")
    axc.plot(xx, expr["CUBN"], "-s", color=kv.SERIES[1], ms=3.2, lw=1.0,
             label="CUBN (target)")
    axc.set_xticks(xx); axc.set_xticklabels(states, fontsize=5.5)
    axc.set_ylabel("mean log expression")
    axc.legend(loc="center left", bbox_to_anchor=(0.02, 0.55),
               handlelength=1.4, borderpad=0.2)
    kv.despine(axc); kv.hairline_grid(axc, "y"); kv.panel(axc, "c", dx=-0.24)
    _t = axc.set_title("Both fall in injured PT", loc="left", pad=4)
    drop_h = expr["HNF4A"][0] - expr["HNF4A"][-1]
    drop_c = expr["CUBN"][0] - expr["CUBN"][-1]
    axc.set_title(f"Both fall in injured PT\n"
                  f"PT_S1 → iPT:  HNF4A −{drop_h:.2f},  CUBN −{drop_c:.2f}",
                  loc="left", pad=4, fontsize=kv.FS_TITLE)

    # ---- d  panel-wide promoter accessibility ------------------------------
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    subd = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, :],
                                   height_ratios=[0.035, 1], hspace=0.02)
    axd = fig.add_subplot(subd[1])
    cols_order = [c for c in ["Podo", "PEC", "Endo_GC", "PT_S1", "PT_S2",
                              "PT_S3", "iPT", "Des-Thin_Limb", "C_TAL",
                              "M_TAL", "DCT1", "DCT2", "CNT", "PC", "IC_A",
                              "IC_B", "Endo_Peritubular", "Fibroblast_1",
                              "Fibroblast_2", "GS_Stromal", "Mac", "CD4T",
                              "B_Naive"] if c in P.columns]
    Pl = np.log2(P[cols_order].clip(lower=0.05))
    ordr = Pl.mean(axis=1).sort_values(ascending=False).index
    Pl = Pl.loc[ordr]
    # Gutters between lineage blocks. imshow has no column spacing, so insert
    # NaN spacer columns and render them as the chart surface.
    GAPCOLS = 1
    cols_draw, xticks, xlabels, blocks = [], [], [], []
    for name, mem in kv.LINEAGE.items():
        mem = [c for c in mem if c in cols_order]
        if not mem:
            continue
        start = len(cols_draw)
        for c in mem:
            xticks.append(len(cols_draw)); xlabels.append(c)
            cols_draw.append(c)
        blocks.append((name, start, len(cols_draw) - 1, len(mem)))
        cols_draw.extend([None] * GAPCOLS)
    while cols_draw and cols_draw[-1] is None:
        cols_draw.pop()

    D = np.full((len(Pl), len(cols_draw)), np.nan)
    for j, c in enumerate(cols_draw):
        if c is not None:
            D[:, j] = Pl[c].values

    cmap = kv.cmap_div().copy(); cmap.set_bad(kv.SURFACE)
    lim = float(np.nanpercentile(np.abs(Pl.values), 98))
    im = axd.imshow(np.ma.masked_invalid(D), aspect="auto", cmap=cmap,
                    vmin=-lim, vmax=lim, interpolation="nearest")

    SHORT = {"Glomerulus": "Glom", "Proximal tubule": "Prox tubule",
             "Loop of Henle": "LoH", "Distal / collecting": "Distal/CD",
             "Endothelium": "Endo", "Stroma": "Stroma", "Immune": "Immune"}
    axband = fig.add_subplot(subd[0], sharex=axd)
    for k, (name, x_lo, x_hi, n_mem) in enumerate(blocks):
        axband.add_patch(plt.Rectangle((x_lo - 0.5, 0), n_mem, 1, lw=0,
                                       color=kv.SERIES[k % len(kv.SERIES)]))
        if n_mem >= 2:
            axband.text((x_lo + x_hi) / 2, 0.5,
                        name if n_mem >= 5 else SHORT[name], ha="center",
                        va="center", fontsize=kv.FS_NOTE, color="white",
                        fontweight="bold", clip_on=True)
    axband.set_ylim(0, 1); axband.set_xlim(-0.5, len(cols_draw) - 0.5)
    axband.axis("off")
    # mark genes carrying an HNF4A motif in the promoter
    motif_genes = set(hits["sequence_name"])
    xm = len(cols_draw) - 0.1
    for i, g in enumerate(Pl.index):
        if g in motif_genes:
            axd.plot([xm], [i], marker="D", ms=1.5, color=kv.SERIES[0],
                     clip_on=False, zorder=6)
    axd.set_xticks(xticks)
    axd.set_xticklabels(xlabels, rotation=90, fontsize=5)
    axd.set_yticks(np.arange(len(Pl))); axd.set_yticklabels(Pl.index, fontsize=4.0)
    axd.set_yticks(np.arange(-.5, len(Pl)), minor=True)
    axd.set_xticks(np.arange(-.5, len(cols_draw)), minor=True)
    axd.grid(which="minor", color=kv.SURFACE, lw=0.4)
    axd.tick_params(which="minor", length=0); axd.tick_params(length=1.2)
    for s in axd.spines.values():
        s.set_visible(False)
    kv.panel(axband, "d", dx=-0.062, dy=3.0)
    # A motif alone is weak evidence -- HNF4A's site is degenerate and a 2 kb
    # window yields many chance matches. The informative conjunction is
    # PT-restricted accessibility AND a motif, not either alone.
    axband.set_title(
        "Promoter accessibility of every panel gene (snATAC, control donors)\n"
        f"{len(motif_genes)}/{len(Pl)} promoters carry an HNF4A motif "
        "(FIMO p<1e-4), marked ◆ at right.\nA motif in a 2 kb window is weak "
        "alone; the informative signal is PT-restricted accessibility with a "
        "motif.", loc="left", pad=4, fontsize=kv.FS_TITLE)
    # A multi-line title anchors on the FIRST line's baseline, so lines 2+ grow
    # downward -- into the lineage band on an axes this short. va="bottom" makes
    # the whole block sit above the anchor instead.
    axband._left_title.set_va("bottom")

    cax = fig.add_axes([0.978, 0.09, 0.008, 0.15])
    cb = fig.colorbar(im, cax=cax); cb.outline.set_visible(False)
    cb.set_label("log2 accessibility\n(signal / background)", fontsize=kv.FS_NOTE)
    cb.ax.tick_params(labelsize=kv.FS_NOTE, length=1.5)

    kv.save(fig, "Figure4_regulatory", caption=(
        "Figure 4 | Chromatin context of the panel. Per-cell-type snATAC "
        "coverage from Abedini et al. 2024 (control donors), quantified over "
        "TSS +/-1 kb and divided by each track's genome-wide mean so tracks of "
        "differing depth are comparable. a, CUBN promoter, all 23 cell types. "
        "b, UMOD promoter as a positive control. c, Mean expression of HNF4A "
        "and CUBN across proximal-tubule states. d, All panel genes; diamonds "
        "mark promoters containing an HNF4A motif (HOCOMOCO v11 "
        "HNF4A_HUMAN.H11MO.0.A, FIMO p<1e-4). These are coverage tracks, not "
        "called peaks, and a motif match is sequence evidence of a candidate "
        "site, not evidence of binding."))

    print(f"\n  CUBN top: {', '.join(f'{c} {v:.1f}x' for c, v in cu.head(5).items())}")
    print(f"  UMOD top: {', '.join(f'{c} {v:.1f}x' for c, v in um.head(4).items())}")


if __name__ == "__main__":
    main()
