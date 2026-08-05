#!/usr/bin/env python
"""
Figure 2 -- Q3 "is expression enriched in disease-relevant cell types?"
            Q4 "does expression change in CKD vs healthy kidney?"

Two contrasts, kept separate because they answer different questions:

  Q3  injured vs healthy proximal tubule (iPT vs PT_S1/S2/S3), within-tissue.
      A CELL-STATE contrast; it does not use donor group at all.

  Q4  disease vs Control donors, on SAMPLE-LEVEL PSEUDOBULK, with DKD and HKD
      tested SEPARATELY. Donors are the unit of replication; testing per cell
      would pseudoreplicate (280k cells from 36 donors are not 280k independent
      observations). Each contrast is FDR-corrected on its own -- sharing one
      correction across both would make each contrast's threshold depend on
      the other.

THREE CORRECTIONS relative to the first version of this figure, all found by
independent re-derivation from the raw h5ad:

 1. GROUPING uses obs['group'] (Control / DKD / HKD), not obs['Status'].
    obs['Status'] is internally inconsistent: donor HK2770's scRNA cells are
    labelled Disease while its snRNA cells are labelled Control. obs['group']
    is consistent for every donor (asserted below).

 2. The composition test previously de-duplicated on (sample, status), which
    returned 37 rows for 36 donors -- HK2770 appeared in BOTH arms and was
    compared against itself. De-duplication is now on `sample` alone.

 3. The pseudobulk offset is the TRUE library size (row-sum of layers/counts
    over all 34,733 genes, precomputed in 00). The previous offset,
    obs['nCount_RNA'], is not a raw UMI total -- it is integer-valued for only
    5.6% of RNA cells and under-counts the true total by ~2% for scRNA but ~11%
    for snRNA, importing a modality bias into a modality-mixed pseudobulk.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).parent))
import kidneyviz as kv  # noqa: E402

MIN_CELLS_PB = 20     # cells needed for a (donor, cell type) pseudobulk point
MIN_SAMPLES = 4       # donors needed per arm to test a gene
CONTRASTS = [("DKD", "Control"), ("HKD", "Control")]
SHORT = {"Glomerulus": "Glom", "Proximal tubule": "Prox tubule",
         "Loop of Henle": "LoH", "Distal / collecting": "Distal/CD",
         "Endothelium": "Endo", "Stroma": "Stroma", "Immune": "Immune"}
GCOL = {"Control": kv.SERIES[0], "DKD": kv.SERIES[1], "HKD": kv.SERIES[3]}
GMARK = {"Control": "o", "DKD": "s", "HKD": "^"}


def lineage_layout(cts, gap=1.3):
    """Column positions with a gutter between lineage blocks."""
    xpos, blocks, cur = {}, [], 0.0
    for name, mem in kv.LINEAGE.items():
        mem = [c for c in mem if c in cts]
        if not mem:
            continue
        for c in mem:
            xpos[c] = cur
            cur += 1.0
        blocks.append((name, xpos[mem[0]], xpos[mem[-1]], len(mem)))
        cur += gap
    return xpos, blocks, cur - gap


def draw_band(fig, spec, sharex, blocks, xmax):
    ax = fig.add_subplot(spec, sharex=sharex)
    for k, (name, lo, hi, n) in enumerate(blocks):
        ax.add_patch(plt.Rectangle((lo - 0.5, 0), n, 1, lw=0,
                                   color=kv.SERIES[k % len(kv.SERIES)]))
        if n >= 2:
            ax.text((lo + hi) / 2, 0.5, name if n >= 5 else SHORT[name],
                    ha="center", va="center", fontsize=kv.FS_NOTE,
                    color="white", fontweight="bold", clip_on=True)
    ax.set_ylim(0, 1); ax.set_xlim(-0.9, xmax - 0.1); ax.axis("off")
    return ax


def main():
    kv.use_style()
    z = np.load(kv.DERIVED / "panel_matrix.npz", allow_pickle=True)
    obs = pd.read_parquet(kv.DERIVED / "panel_obs.parquet")
    X, C, G = z["lognorm"], z["counts"], np.array(list(z["genes"]))
    rna = obs["tech"].isin(kv.RNA_TECHS).values
    ct, smp = obs["celltype"].values, obs["sample"].values
    lib = obs["libsize_true"].values

    # ---- donor table: one row per donor, group resolved once ---------------
    o = obs[rna]
    chk = o.groupby("sample")["group"].nunique()
    assert (chk == 1).all(), f"group not donor-consistent: {list(chk[chk>1].index)}"
    donors = o[["sample", "group"]].drop_duplicates("sample").reset_index(drop=True)
    assert len(donors) == o["sample"].nunique(), "donor de-duplication is wrong"
    grp_of = dict(zip(donors["sample"], donors["group"]))
    print("donors per group:", donors.groupby("group").size().to_dict(),
          f"(total {len(donors)})")

    # =================================================== Q3  iPT vs PT ======
    healthy_pt = rna & np.isin(ct, ["PT_S1", "PT_S2", "PT_S3"])
    injured_pt = rna & (ct == "iPT")
    print(f"healthy PT n={healthy_pt.sum():,}  iPT n={injured_pt.sum():,}")
    lfc = np.array([X[injured_pt, i].mean() - X[healthy_pt, i].mean()
                    for i in range(len(G))])
    q3 = pd.DataFrame({"gene": G, "delta_log_iPT_vs_PT": lfc}) \
        .sort_values("delta_log_iPT_vs_PT")
    q3.to_csv(kv.SOURCE / "Q3_iPT_vs_healthyPT.csv", index=False)

    # ======================================= Q4  pseudobulk, TRUE offset ====
    rows = []
    for c in kv.CT_ORDER:
        for s in donors["sample"]:
            m = rna & (ct == c) & (smp == s)
            n = int(m.sum())
            if n < MIN_CELLS_PB:
                continue
            rows.append(dict(celltype=c, sample=s, n_cells=n,
                             group=grp_of[s], _lib=float(lib[m].sum()),
                             **{g: v for g, v in zip(G, C[m].sum(axis=0))}))
    pb = pd.DataFrame(rows)
    for g in G:                       # CPM of the donor's own transcriptome
        pb[g] = np.log1p(pb[g] / pb["_lib"].replace(0, np.nan) * 1e6)
    pb.to_csv(kv.SOURCE / "Q4_pseudobulk_by_sample_celltype.csv", index=False)
    print(f"pseudobulk points: {len(pb)} across {pb.celltype.nunique()} cell types")

    res = {}
    for dis, ref in CONTRASTS:
        out = []
        for c in pb["celltype"].unique():
            sub = pb[pb.celltype == c]
            a, b = sub[sub.group == dis], sub[sub.group == ref]
            if len(a) < MIN_SAMPLES or len(b) < MIN_SAMPLES:
                continue
            for g in G:
                x, y = a[g].dropna(), b[g].dropna()
                if len(x) < MIN_SAMPLES or len(y) < MIN_SAMPLES:
                    continue
                out.append(dict(celltype=c, gene=g,
                                delta=float(x.mean() - y.mean()),
                                p=float(stats.mannwhitneyu(
                                    x, y, alternative="two-sided").pvalue),
                                n_dis=len(x), n_ref=len(y)))
        df = pd.DataFrame(out)
        df["fdr"] = multipletests(df["p"], method="fdr_bh")[1]
        res[dis] = df
        df.to_csv(kv.SOURCE / f"Q4_{dis}_vs_Control_pseudobulk.csv", index=False)
        print(f"  {dis} vs Control: {len(df)} tests across "
              f"{df.celltype.nunique()} cell types, "
              f"{(df.fdr < 0.05).sum()} at FDR<0.05, min p={df.p.min():.3g}, "
              f"min FDR={df.fdr.min():.3g}")

    # ---- composition, de-duplicated on donor ONLY --------------------------
    comp = (pd.DataFrame({"celltype": ct[rna], "sample": smp[rna]})
            .groupby(["sample", "celltype"]).size().rename("n").reset_index())
    comp["frac"] = comp["n"] / comp.groupby("sample")["n"].transform("sum")
    comp["group"] = comp["sample"].map(grp_of)
    comp.to_csv(kv.SOURCE / "Q4_composition_by_sample.csv", index=False)

    crows = []
    for c in comp["celltype"].unique():
        sub = comp[comp.celltype == c][["sample", "frac"]]
        full = donors.merge(sub, on="sample", how="left").fillna({"frac": 0.0})
        full["pct"] = full["frac"] * 100
        rec = dict(celltype=c,
                   median_Control=float(full[full.group == "Control"]["pct"].median()))
        for dis, ref in CONTRASTS:
            x = full[full.group == dis]["pct"]
            y = full[full.group == ref]["pct"]
            rec[f"median_{dis}"] = float(x.median())
            rec[f"p_{dis}"] = (float(stats.mannwhitneyu(
                x, y, alternative="two-sided").pvalue)
                if len(x) >= MIN_SAMPLES and len(y) >= MIN_SAMPLES else np.nan)
        crows.append(rec)
    ctest = pd.DataFrame(crows)
    for dis, _ in CONTRASTS:                      # correct each contrast alone
        ok = ctest[f"p_{dis}"].notna()
        ctest.loc[ok, f"fdr_{dis}"] = multipletests(
            ctest.loc[ok, f"p_{dis}"], method="fdr_bh")[1]
    ctest = ctest.sort_values("p_DKD")
    ctest.to_csv(kv.SOURCE / "Q4_composition_test.csv", index=False)
    print("composition, top 6 by DKD P:")
    print(ctest.head(6)[["celltype", "median_Control", "median_DKD",
                         "median_HKD", "p_DKD", "p_HKD"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    # ======================================================= FIGURE 2 =======
    fig = plt.figure(figsize=(kv.W2, 13.4))
    gs = GridSpec(3, 2, height_ratios=[0.72, 1.30, 1.30],
                  width_ratios=[1.25, 1], hspace=0.34, wspace=0.26,
                  left=0.075, right=0.915, top=0.950, bottom=0.050)

    # ---- a  iPT vs healthy PT ---------------------------------------------
    axa = fig.add_subplot(gs[0, 0])
    dd = pd.concat([q3.head(12), q3.tail(12)])
    y = np.arange(len(dd))
    axa.barh(y, dd["delta_log_iPT_vs_PT"], height=0.75, lw=0,
             color=[kv.STATUS["critical"] if v > 0 else kv.SERIES[0]
                    for v in dd["delta_log_iPT_vs_PT"]])
    axa.set_yticks(y); axa.set_yticklabels(dd["gene"], fontsize=5)
    axa.axvline(0, color=kv.INK, lw=0.5)
    axa.set_xlabel("Δ mean log-normalised expression   (iPT − healthy PT)")
    axa.set_ylim(-1, len(dd))
    kv.despine(axa); kv.hairline_grid(axa, "x"); kv.panel(axa, "a", dx=-0.16)
    axa.set_title("Q3  Injured vs healthy proximal tubule\n"
                  f"iPT n={injured_pt.sum():,}, healthy PT n={healthy_pt.sum():,}"
                  " — cell state, independent of donor group",
                  loc="left", pad=4, fontsize=kv.FS_TITLE)

    # ---- b  composition, three groups -------------------------------------
    axb = fig.add_subplot(gs[0, 1])
    key = [c for c in ["PT_S1", "PT_S2", "PT_S3", "iPT", "C_TAL", "M_TAL",
                       "Fibroblast_1", "MyoFib/VSMC", "Mac", "CD4T"]
           if c in set(ctest.celltype)]
    look = ctest.set_index("celltype")
    key = sorted(key, key=lambda c: min(look.loc[c, "p_DKD"],
                                        look.loc[c, "p_HKD"]))
    rng = np.random.default_rng(0)
    FLOOR = 0.02
    for i, c in enumerate(key):
        sub = comp[comp.celltype == c][["sample", "frac"]]
        full = donors.merge(sub, on="sample", how="left").fillna({"frac": 0.0})
        full["pct"] = (full["frac"] * 100).clip(lower=FLOOR)
        for k, g in enumerate(["Control", "DKD", "HKD"]):
            v = full[full.group == g]["pct"].values
            pos = i + (k - 1) * 0.27
            # Box (median + IQR, 1.5*IQR whiskers) behind, individual donors on
            # top. With n = 7-16 per arm the points must stay visible: a box
            # alone would hide how few donors carry each comparison.
            bp = axb.boxplot([v], positions=[pos], widths=0.21, vert=False,
                             patch_artist=True, showfliers=False,
                             manage_ticks=False, zorder=2,
                             medianprops=dict(color=kv.INK, lw=0.7),
                             whiskerprops=dict(lw=0.4, color=kv.INK),
                             capprops=dict(lw=0.4, color=kv.INK),
                             boxprops=dict(lw=0.35, edgecolor=kv.INK))
            for b in bp["boxes"]:
                b.set_facecolor(GCOL[g]); b.set_alpha(0.30)
            yy = pos + rng.normal(0, 0.035, len(v))
            axb.scatter(v, yy, s=3.0, c=GCOL[g], alpha=0.85, lw=0.15,
                        edgecolors=kv.SURFACE, zorder=4)
            axb.scatter([np.median(v)], [pos], s=15, c=GCOL[g],
                        marker=GMARK[g], edgecolors=kv.INK, linewidths=0.4,
                        zorder=5)
    axb.set_xscale("log"); axb.set_xlim(FLOOR * 0.7, 6000)
    axb.set_yticks(np.arange(len(key)))
    axb.set_yticklabels([{"Fibroblast_1": "Fib_1",
                          "MyoFib/VSMC": "MyoFib"}.get(c, c) for c in key],
                        fontsize=5.5)
    axb.invert_yaxis(); axb.set_ylim(len(key) - 0.4, -1.5)
    axb.set_xlabel("% of a donor's RNA cells  (log scale)")
    kv.despine(axb); kv.hairline_grid(axb, "x"); kv.panel(axb, "b", dx=-0.145)
    PX = [110, 1200]
    for j, (dis, _) in enumerate(CONTRASTS):
        axb.text(PX[j], -1.05, f"P {dis}", ha="center", va="center",
                 fontsize=kv.FS_NOTE, color=kv.INK, fontweight="bold")
        for i, c in enumerate(key):
            p, f = look.loc[c, f"p_{dis}"], look.loc[c, f"fdr_{dis}"]
            axb.text(PX[j], i, f"{p:.3f}" if p >= 0.001 else f"{p:.0e}",
                     ha="center", va="center", fontsize=kv.FS_NOTE,
                     color=kv.STATUS["critical"] if f < 0.05 else kv.INK2,
                     fontweight="bold" if f < 0.05 else "normal")
    nC = int((donors.group == "Control").sum())
    nD = int((donors.group == "DKD").sum())
    nH = int((donors.group == "HKD").sum())
    axb.set_title("Q4  Cell-type composition, by disease\n"
                  f"Control n={nC}, DKD n={nD}, HKD n={nH} donors; "
                  "box = median and IQR, whiskers 1.5xIQR; every donor plotted",
                  loc="left", pad=4, fontsize=kv.FS_TITLE)
    hs = [plt.Line2D([], [], marker=GMARK[g], ls="", ms=3.2, mfc=GCOL[g],
                     mec=kv.INK, mew=0.4, label=g)
          for g in ["Control", "DKD", "HKD"]]
    axb.legend(handles=hs, loc="lower left", bbox_to_anchor=(0.0, -0.26),
               ncol=3, handletextpad=0.3, columnspacing=1.0, borderpad=0.2)

    # ---- c, d  the two disease contrasts ----------------------------------
    for row, (dis, ref) in enumerate(CONTRASTS, start=1):
        df = res[dis]
        cts = [c for c in kv.CT_ORDER if c in set(df.celltype)]
        piv = df.pivot(index="gene", columns="celltype",
                       values="delta").reindex(index=list(G), columns=cts)
        sig = df.pivot(index="gene", columns="celltype",
                       values="fdr").reindex(index=list(G), columns=cts)
        pct = np.zeros((len(G), len(cts)))
        for j, c in enumerate(cts):
            m = rna & (ct == c)
            pct[:, j] = (C[m] > 0).mean(axis=0) * 100 if m.sum() else 0.0

        sd = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[row, :],
                                     height_ratios=[0.032, 1], hspace=0.02)
        ax = fig.add_subplot(sd[1])
        xpos, blocks, xmax = lineage_layout(cts)
        XS = np.array([xpos[c] for c in cts])
        lim = float(np.nanpercentile(np.abs(piv.values), 98))
        gx, gy = np.meshgrid(XS, np.arange(len(G)))
        im = ax.scatter(gx.ravel(), gy.ravel(),
                        s=(pct.ravel() / 100.0) * 10.0 + 0.10,
                        c=piv.values.ravel(), cmap=kv.cmap_div(),
                        vmin=-lim, vmax=lim, linewidths=0.12,
                        edgecolors=kv.INK, zorder=3)
        ys, xs = np.where(sig.values < 0.05)
        ax.scatter(XS[xs], ys, s=5.5, facecolors="none", edgecolors=kv.INK,
                   linewidths=0.45, zorder=5)
        ax.set_xlim(-0.9, xmax - 0.1); ax.set_ylim(len(G) - 0.3, -0.7)
        ax.set_xticks(XS); ax.set_xticklabels(cts, rotation=90, fontsize=4.8)
        ax.set_yticks(np.arange(len(G))); ax.set_yticklabels(G, fontsize=4.0)
        ax.set_yticks(np.arange(-.5, len(G)), minor=True)
        ax.grid(which="minor", axis="y", color=kv.GRID, lw=0.25, zorder=0)
        ax.set_axisbelow(True); ax.tick_params(which="minor", length=0)
        ax.tick_params(length=1.2)
        for sp in ax.spines.values():
            sp.set_visible(False)

        band = draw_band(fig, sd[0], ax, blocks, xmax)
        nsig = int((df.fdr < 0.05).sum())
        band.set_title(
            f"Q4  {dis} − Control, per cell type "
            f"(pseudobulk, two-sided Mann–Whitney, FDR within this contrast)\n"
            + (f"{nsig} of {len(df)} gene × cell-type tests reach FDR < 0.05"
               if nsig else
               f"none of {len(df)} gene × cell-type tests reach FDR < 0.05 "
               f"({int(df.n_dis.max())} {dis} vs {int(df.n_ref.max())} Control "
               f"donors; smallest FDR {df.fdr.min():.2f})"),
            loc="left", pad=4, fontsize=kv.FS_TITLE)
        band._left_title.set_va("bottom")
        kv.panel(band, "cd"[row - 1], dx=-0.062, dy=2.4)
        if row == 1:
            cax = fig.add_axes([0.930, 0.30, 0.008, 0.10])
            cb = fig.colorbar(im, cax=cax); cb.outline.set_visible(False)
            cb.set_label("Δ log expression\n(disease − Control)",
                         fontsize=kv.FS_NOTE)
            cb.ax.tick_params(labelsize=kv.FS_NOTE, length=1.5)
            sax = fig.add_axes([0.928, 0.43, 0.045, 0.065]); sax.axis("off")
            for k, q in enumerate([5, 25, 50, 100]):
                sax.scatter([0.16], [k * 0.26], s=(q / 100) * 10.0 + 0.10,
                            c=kv.MUTED, linewidths=0.12, edgecolors=kv.INK)
                sax.text(0.45, k * 0.26, f"{q}%", fontsize=kv.FS_NOTE,
                         va="center", color=kv.INK)
            sax.text(0, 1.12, "% of cells\nexpressing", fontsize=kv.FS_NOTE,
                     va="bottom", color=kv.INK)
            sax.set_xlim(0, 1); sax.set_ylim(-0.2, 0.95)

    kv.save(fig, "Figure2_disease", caption=(
        "Figure 2 | Panel-gene behaviour in kidney disease. Donor group is "
        "taken from obs['group'] (Control/DKD/HKD), which is internally "
        "consistent, rather than obs['Status'], which is not. a, Difference in "
        "mean log-normalised expression between injured (iPT) and healthy "
        "proximal tubule; 12 most decreased and 12 most increased of 56 genes. "
        "b, Per-donor cell-type fractions by group. Boxes show median and "
        "interquartile range with 1.5x IQR whiskers; every donor is plotted as "
        "a point and the large marker is the median. Exact two-sided "
        "Mann-Whitney P per "
        "contrast, each Benjamini-Hochberg corrected within its own contrast. "
        "c, DKD minus Control and d, HKD minus Control differences in log "
        "pseudobulk expression. Counts were summed per (donor, cell type), "
        "divided by the true library size (row-sum of raw counts over all "
        "34,733 genes) and log-transformed; donors are the unit of "
        "replication. Dot size is the percentage of that cell type's cells "
        "expressing the gene, so differences for genes that are off in a cell "
        "type shrink out of visual prominence. Open rings mark FDR < 0.05."))


if __name__ == "__main__":
    main()
