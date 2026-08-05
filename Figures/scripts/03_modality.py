#!/usr/bin/env python
"""
Figure 3 -- Q5 "are expression patterns consistent across scRNA and snRNA?"

Single-cell and single-nucleus assays sample different RNA pools. snRNA-seq
reads the nuclear transcriptome, so it over-represents long/intron-rich and
actively transcribed genes and under-represents mature cytoplasmic mRNA;
scRNA-seq requires tissue dissociation, which lyses fragile cell types and
induces a stress/immediate-early programme (FOS, JUN, immune activation).

Disagreement between the two is therefore expected for specific gene classes
and is informative about the assay, not merely noise. This figure quantifies
it per gene and per cell type instead of assuming concordance.

The comparison is made WITHIN matched cell types, because the two assays also
recover different cell-type proportions -- comparing pooled means would confound
assay chemistry with composition.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import kidneyviz as kv  # noqa: E402

MIN_CELLS = 30


def main():
    kv.use_style()
    z = np.load(kv.DERIVED / "panel_matrix.npz", allow_pickle=True)
    obs = pd.read_parquet(kv.DERIVED / "panel_obs.parquet")
    X, C, G = z["lognorm"], z["counts"], np.array(list(z["genes"]))
    ct = obs["celltype"].values
    tech = obs["tech"].values

    sc = tech == "SC_RNA"
    sn = tech == "SN_RNA"
    print(f"scRNA {sc.sum():,}   snRNA {sn.sum():,}")

    # cell types with enough cells in BOTH assays -- the matched set
    shared = [c for c in kv.CT_ORDER
              if (sc & (ct == c)).sum() >= MIN_CELLS
              and (sn & (ct == c)).sum() >= MIN_CELLS]
    print(f"cell types present in both assays: {len(shared)}")

    rows = []
    for c in shared:
        a, b = sc & (ct == c), sn & (ct == c)
        for i, g in enumerate(G):
            rows.append(dict(celltype=c, gene=g,
                             sc_mean=float(X[a, i].mean()),
                             sn_mean=float(X[b, i].mean()),
                             sc_pct=float((C[a, i] > 0).mean() * 100),
                             sn_pct=float((C[b, i] > 0).mean() * 100),
                             n_sc=int(a.sum()), n_sn=int(b.sum())))
    m = pd.DataFrame(rows)
    m["delta"] = m["sn_mean"] - m["sc_mean"]        # + = higher in snRNA
    m.to_csv(kv.SOURCE / "Q5_sc_vs_sn_by_celltype.csv", index=False)

    # per-gene summary across matched cell types
    per_gene = (m.groupby("gene")
                  .apply(lambda d: pd.Series({
                      "sc_mean": d.sc_mean.mean(),
                      "sn_mean": d.sn_mean.mean(),
                      "delta": d.delta.mean(),
                      "rho": stats.spearmanr(d.sc_mean, d.sn_mean).statistic
                             if d.sc_mean.std() > 0 and d.sn_mean.std() > 0
                             else np.nan,
                      "p": stats.wilcoxon(d.sn_mean, d.sc_mean).pvalue
                           if len(d) > 5 else np.nan,
                  }), include_groups=False)
                  .reset_index())
    from statsmodels.stats.multitest import multipletests
    ok = per_gene["p"].notna()
    per_gene.loc[ok, "fdr"] = multipletests(per_gene.loc[ok, "p"],
                                            method="fdr_bh")[1]
    per_gene = per_gene.sort_values("delta")
    per_gene.to_csv(kv.SOURCE / "Q5_sc_vs_sn_per_gene.csv", index=False)

    rho_all = stats.spearmanr(m.sc_mean, m.sn_mean).statistic
    print(f"overall Spearman rho (gene x celltype) = {rho_all:.3f}")
    print("higher in snRNA:",
          ", ".join(f"{r.gene} {r.delta:+.2f}" for r in per_gene.tail(6).itertuples()))
    print("higher in scRNA:",
          ", ".join(f"{r.gene} {r.delta:+.2f}" for r in per_gene.head(6).itertuples()))

    # ======================================================= FIGURE 3 =======
    fig = plt.figure(figsize=(kv.W2, 9.6))
    gs = GridSpec(2, 3, height_ratios=[0.72, 1.55], width_ratios=[1, 1, 1.15],
                  hspace=0.30, wspace=0.34,
                  left=0.07, right=0.955, top=0.945, bottom=0.075)

    # ---- a  scatter, every gene x cell type -------------------------------
    axa = fig.add_subplot(gs[0, 0])
    axa.scatter(m.sc_mean, m.sn_mean, s=2.2, c=kv.SERIES[0], alpha=0.35, lw=0)
    hi = float(max(m.sc_mean.max(), m.sn_mean.max())) * 1.05
    axa.plot([0, hi], [0, hi], color=kv.INK, lw=0.5, ls=(0, (3, 2)), zorder=4)
    axa.set_xlim(0, hi); axa.set_ylim(0, hi)
    axa.set_xlabel("scRNA  mean log expression")
    axa.set_ylabel("snRNA  mean log expression")
    axa.set_title("Matched cell types agree overall", loc="left", pad=4)
    axa.text(0.04, 0.96, f"ρ = {rho_all:.2f}\n{len(m):,} gene × cell type",
             transform=axa.transAxes, ha="left", va="top",
             fontsize=kv.FS_NOTE, color=kv.INK2)
    kv.despine(axa); kv.panel(axa, "a", dx=-0.20)

    # ---- b  per-gene bias --------------------------------------------------
    axb = fig.add_subplot(gs[0, 1])
    pg = pd.concat([per_gene.head(11), per_gene.tail(11)])
    y = np.arange(len(pg))
    cols = [kv.SERIES[2] if v > 0 else kv.SERIES[1] for v in pg["delta"]]
    axb.barh(y, pg["delta"], color=cols, height=0.75, lw=0)
    axb.set_yticks(y); axb.set_yticklabels(pg["gene"], fontsize=5)
    axb.axvline(0, color=kv.INK, lw=0.5)
    axb.set_xlabel("Δ mean log expression  (snRNA − scRNA)\n"
                   "← higher in scRNA          higher in snRNA →")
    axb.set_ylim(-1, len(pg))
    kv.despine(axb); kv.hairline_grid(axb, "x"); kv.panel(axb, "b", dx=-0.24)
    axb.set_title("Per-gene assay bias", loc="left", pad=4)

    # ---- c  detection rate -------------------------------------------------
    axc = fig.add_subplot(gs[0, 2])
    axc.scatter(m.sc_pct, m.sn_pct, s=2.2, c=kv.SERIES[6], alpha=0.35, lw=0)
    axc.plot([0, 100], [0, 100], color=kv.INK, lw=0.5, ls=(0, (3, 2)), zorder=4)
    axc.set_xlim(0, 100); axc.set_ylim(0, 100)
    axc.set_xlabel("scRNA  % cells detected")
    axc.set_ylabel("snRNA  % cells detected")
    axc.set_title("snRNA detects more per cell type", loc="left", pad=4)
    frac = float((m.sn_pct > m.sc_pct).mean() * 100)
    axc.text(0.04, 0.96, f"{frac:.0f}% of points\nabove the diagonal",
             transform=axc.transAxes, ha="left", va="top",
             fontsize=kv.FS_NOTE, color=kv.INK2)
    kv.despine(axc); kv.panel(axc, "c", dx=-0.20)

    # ---- d  where the disagreement lives -----------------------------------
    axd = fig.add_subplot(gs[1, :])
    piv = m.pivot(index="gene", columns="celltype",
                  values="delta").reindex(index=list(per_gene["gene"]),
                                          columns=[c for c in shared])
    lim = float(np.nanpercentile(np.abs(piv.values), 98))
    im = axd.imshow(piv.values, aspect="auto", cmap=kv.cmap_div(),
                    vmin=-lim, vmax=lim, interpolation="nearest")
    axd.set_xticks(np.arange(len(shared)))
    axd.set_xticklabels(shared, rotation=90, fontsize=4.8)
    axd.set_yticks(np.arange(len(piv)))
    axd.set_yticklabels(piv.index, fontsize=4.4)
    axd.set_yticks(np.arange(-.5, len(piv)), minor=True)
    axd.set_xticks(np.arange(-.5, len(shared)), minor=True)
    axd.grid(which="minor", color=kv.SURFACE, lw=0.4)
    axd.tick_params(which="minor", length=0); axd.tick_params(length=1.2)
    for s in axd.spines.values():
        s.set_visible(False)
    kv.panel(axd, "d", dx=-0.062, dy=1.015)
    axd.set_title("snRNA − scRNA, per gene and cell type "
                  "(rows ordered by mean bias)", loc="left", pad=4)

    cax = fig.add_axes([0.965, 0.11, 0.008, 0.15])
    cb = fig.colorbar(im, cax=cax); cb.outline.set_visible(False)
    cb.set_label("Δ log expression\n(snRNA − scRNA)", fontsize=kv.FS_NOTE)
    cb.ax.tick_params(labelsize=kv.FS_NOTE, length=1.5)

    kv.save(fig, "Figure3_scRNA_vs_snRNA", caption=(
        "Figure 3 | Single-cell versus single-nucleus measurement of the panel. "
        "All comparisons are made within matched cell types (>=30 cells in both "
        f"assays; n = {len(shared)} cell types), because the two assays recover "
        "different cell-type proportions and pooled means would confound assay "
        "with composition. a, Mean log-normalised expression per gene and cell "
        "type; dashed line is identity. b, Mean difference per gene across "
        "matched cell types, 11 most scRNA-biased and 11 most snRNA-biased. "
        "c, Detection rate per gene and cell type. d, Difference per gene and "
        "cell type, rows ordered by mean bias. Source: "
        "Q5_sc_vs_sn_by_celltype.csv, Q5_sc_vs_sn_per_gene.csv."))


if __name__ == "__main__":
    main()
