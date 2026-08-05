#!/usr/bin/env python
"""
00 -- Extract the 56-gene panel out of the 18 GB PreSCVI h5ad.

Streams the CSR matrix once, keeping only the panel's gene columns, and writes a
compact cells x 56 matrix plus the cell metadata every downstream figure needs.

Two matrices are kept and they are NOT interchangeable:
  X       log-normalised values  -> for display / cross-cell comparison
  counts  raw integer counts     -> for pseudobulk aggregation and DE

Modality caveat carried through to every downstream script: for SN_ATAC cells
these values are gene-ACTIVITY (fragments over gene body + promoter), not mRNA.

Run:  <biomni_e1 python> 00_extract_panel.py
"""
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from kidneyviz import DERIVED, H5AD_PRE, PANEL  # noqa: E402

CHUNK = 20_000  # cells per streamed block


def read_cat(h, col):
    """Read a categorical obs column back to string labels."""
    g = h["obs"][col]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([x.decode() if isinstance(x, bytes) else str(x)
                         for x in g["categories"][:]])
        codes = g["codes"][:]
        out = np.where(codes >= 0, cats[np.clip(codes, 0, len(cats) - 1)], "NA")
        return out
    v = g[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in v])


def stream_libsize(h, n_obs):
    """
    TRUE per-cell library size: row-sum of layers/counts over ALL 34,733 genes.

    obs['nCount_RNA'] is NOT this. It is integer-valued for only 5.6% of RNA
    cells and equals the real row-sum for 2 cells in 282,610 -- it is some
    upstream-corrected quantity, and it is biased by modality (true totals run
    ~2.6% above it for scRNA but ~11.1% for snRNA). Using it as a pseudobulk
    offset therefore imports a systematic snRNA-vs-scRNA bias. For SN_ATAC rows
    the same column holds an unrelated quantity entirely.
    """
    grp = h["layers"]["counts"]
    indptr = grp["indptr"][:]
    data = grp["data"]
    out = np.zeros(n_obs, dtype=np.float64)
    t0 = time.time()
    for start in range(0, n_obs, CHUNK):
        stop = min(start + CHUNK, n_obs)
        lo, hi = int(indptr[start]), int(indptr[stop])
        if hi == lo:
            continue
        vals = data[lo:hi]
        counts = np.diff(indptr[start:stop + 1])
        out[start:stop] = np.add.reduceat(
            vals, np.r_[0, np.cumsum(counts)[:-1]]) * (counts > 0)
    print(f"    libsize: done in {time.time()-t0:.0f}s")
    return out


def stream_panel(h, key, col_idx, n_obs):
    """One pass over a CSR matrix, materialising only the panel columns."""
    grp = h[key] if key == "X" else h["layers"][key]
    indptr = grp["indptr"][:]
    data_ds, idx_ds = grp["data"], grp["indices"]

    lookup = np.full(int(idx_ds[:1].max() + 1) if False else 0, -1)  # placeholder
    ncol = int(h["X"].attrs["shape"][1])
    lookup = np.full(ncol, -1, dtype=np.int32)
    lookup[col_idx] = np.arange(len(col_idx), dtype=np.int32)

    out = np.zeros((n_obs, len(col_idx)), dtype=np.float32)
    t0 = time.time()
    for start in range(0, n_obs, CHUNK):
        stop = min(start + CHUNK, n_obs)
        lo, hi = int(indptr[start]), int(indptr[stop])
        if hi == lo:
            continue
        cols = idx_ds[lo:hi]
        keep = lookup[cols] >= 0
        if not keep.any():
            continue
        vals = data_ds[lo:hi][keep]
        tgt = lookup[cols[keep]]
        # row index for each retained nonzero
        rows = np.repeat(np.arange(start, stop),
                         np.diff(indptr[start:stop + 1]))[keep]
        out[rows, tgt] = vals
        if (start // CHUNK) % 5 == 0:
            done = stop / n_obs
            print(f"    {key}: {done:5.1%}  ({time.time()-t0:5.0f}s)", flush=True)
    print(f"    {key}: done in {time.time()-t0:.0f}s")
    return out


def main():
    print(f"opening {H5AD_PRE.name} ({H5AD_PRE.stat().st_size/2**30:.1f} GB)")
    with h5py.File(H5AD_PRE, "r") as h:
        var = np.array([x.decode() for x in h["var"]["_index"][:]])
        n_obs = int(h["X"].attrs["shape"][0])

        pos = {g: i for i, g in enumerate(var)}
        missing = [g for g in PANEL if g not in pos]
        if missing:
            raise SystemExit(f"panel genes absent from object: {missing}")
        col_idx = np.array([pos[g] for g in PANEL], dtype=np.int64)

        print(f"  {n_obs:,} cells x {len(var):,} genes -> panel of {len(PANEL)}")

        obs = pd.DataFrame({
            "cell": [x.decode() for x in h["obs"]["_index"][:]],
            "celltype": read_cat(h, "Cluster_Idents"),
            "tech": read_cat(h, "tech"),
            "status": read_cat(h, "Status"),
            "group": read_cat(h, "group"),
            "sample": read_cat(h, "sample"),
            "sex": read_cat(h, "sex"),
            "orig_ident": read_cat(h, "orig_ident"),
        })
        obs["nCount_RNA"] = pd.to_numeric(read_cat(h, "nCount_RNA"), errors="coerce")

        print("  streaming TRUE library size (row-sum of counts, all genes)...")
        obs["libsize_true"] = stream_libsize(h, n_obs)

        print("  streaming X (log-normalised)...")
        Xn = stream_panel(h, "X", col_idx, n_obs)
        print("  streaming counts (raw)...")
        Xc = stream_panel(h, "counts", col_idx, n_obs)

    # ---- integrity checks, reported rather than assumed --------------------
    print("\n  verification")
    print(f"    X       min={Xn[Xn>0].min():.4f} max={Xn.max():.3f} "
          f"integer={np.allclose(Xn, np.round(Xn))}")
    print(f"    counts  min={Xc[Xc>0].min():.4f} max={Xc.max():.0f} "
          f"integer={np.allclose(Xc, np.round(Xc))}")
    for t in obs["tech"].unique():
        m = (obs["tech"] == t).values
        det = (Xc[m] > 0).mean()
        lt = obs.loc[m, "libsize_true"].values
        nc = obs.loc[m, "nCount_RNA"].values
        with np.errstate(invalid="ignore", divide="ignore"):
            bias = np.nanmedian(lt / nc) - 1
        print(f"    {t:8s} n={m.sum():7,}  panel detection={det:.1%}  "
              f"median true library={np.nanmedian(lt):.0f}  "
              f"nCount_RNA under-counts by {bias:+.1%}")
    li = obs["libsize_true"].values
    print(f"    libsize_true integer-valued: "
          f"{np.isclose(li, np.round(li)).mean():.1%}  (must be 100%)")

    DERIVED.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(DERIVED / "panel_matrix.npz",
                        lognorm=Xn, counts=Xc, genes=np.array(PANEL))
    obs.to_parquet(DERIVED / "panel_obs.parquet")
    print(f"\n  wrote {DERIVED/'panel_matrix.npz'}")
    print(f"  wrote {DERIVED/'panel_obs.parquet'}")


if __name__ == "__main__":
    main()
