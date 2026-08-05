"""
Shared style + paths for the kidney gene-panel figure set.

Palette is the validated reference instance: fixed-order categorical slots,
single-hue sequential (blue), blue<->red diverging with a neutral gray midpoint.
Figures are rendered light-mode only -- they are destined for print / manuscript,
where the dark variant has no consumer.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- paths -----
ROOT = Path("/mnt/beegfs/scratch/miska/nm667/inProgress/dev")
FIGROOT = ROOT / "Figures"
FIGURES = FIGROOT / "figures"
SOURCE = FIGROOT / "source_data"
DERIVED = FIGROOT / "source_data" / "derived"

H5AD_PRE = ROOT / "data/Kidney_SC/GSE211785_Susztak_SC_SN_ATAC_merged_PreSCVI_final.h5ad"
ATAC_DIR = ROOT / "data/Kidney_scATAC"
MOTIF_DB = ROOT / "resources/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
REFGENE = FIGROOT / "source_data" / "hg38.refGene.txt.gz"
FIMO = "/mnt/home1/miska/nm667/meme/bin/fimo"

for _p in (FIGURES, SOURCE, DERIVED):
    _p.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------- the panel ----
# The 56-gene panel this project is built around (from figr_cluster_pipeline).
PANEL = [
    "SLC5A10", "CYFIP2", "RELB", "ACAT1", "ARMC7", "PAX8", "ENPEP", "IL1RL1",
    "DCN", "ALDH8A1", "SLC2A9", "SLC9A3", "UMOD", "PPARA", "AMN", "ACTA2",
    "KCNJ15", "PEPD", "TNF", "EFHD1", "SOST", "SLC17A3", "ALDOB", "MME",
    "LGALS2", "BDH2", "FN1", "BICC1", "ICAM1", "HAVCR1", "TDGF1", "TM7SF3",
    "CD24", "VCAM1", "VIM", "HOXA9", "RGL3", "COL1A1", "JUN", "DPEP1", "TCN2",
    "OCIAD2", "SOX9", "TWIST1", "SNAI1", "BCAM", "CRB2", "CCL2", "CUBN", "FOS",
    "ACSF2", "SHMT1", "AOC1", "HNF4A", "PROM1", "ACY1",
]

# --------------------------------------------------- cell types (verified) --
# The 41 labels below are read off obs['Cluster_Idents'] of the PreSCVI object,
# not assumed. NOTE: the KPMP adaptive/degenerative vocabulary (aPT, aTAL1/2,
# dPT, dC-TAL, aFIB, cycEPI) belongs to the *PostSCVI KPMP-merged* object and is
# NOT present here -- do not reference those labels against this dataset.
CELLTYPES = [
    "Ascending_Thin_LOH", "B_Naive", "B_memory", "Baso/Mast", "CD14_Mono",
    "CD16_Mono", "CD4T", "CD8T", "CNT", "C_TAL", "DCT1", "DCT2",
    "Des-Thin_Limb", "Endo_GC", "Endo_Lymphatic", "Endo_Peritubular",
    "Fibroblast_1", "Fibroblast_2", "GS_Stromal", "IC_A", "IC_B", "M_TAL",
    "Mac", "Macula_Densa", "Mes", "MyoFib/VSMC", "NK", "Neural_Cells",
    "Neutrophil", "PC", "PEC", "PT_S1", "PT_S2", "PT_S3", "Plasma_Cells",
    "Podo", "Prolif_Lym", "RBC", "cDC", "iPT", "pDC",
]

# Anatomical ordering: down the nephron, then vasculature, stroma, immune.
# Ordering the axis by nephron segment (not alphabetically) is what makes a
# segment-restricted expression pattern legible as a band.
LINEAGE = {
    "Glomerulus": ["Podo", "PEC", "Mes", "Endo_GC"],
    "Proximal tubule": ["PT_S1", "PT_S2", "PT_S3", "iPT"],
    "Loop of Henle": ["Des-Thin_Limb", "Ascending_Thin_LOH", "C_TAL", "M_TAL",
                      "Macula_Densa"],
    "Distal / collecting": ["DCT1", "DCT2", "CNT", "PC", "IC_A", "IC_B"],
    "Endothelium": ["Endo_Peritubular", "Endo_Lymphatic"],
    "Stroma": ["Fibroblast_1", "Fibroblast_2", "MyoFib/VSMC", "GS_Stromal",
               "Neural_Cells"],
    "Immune": ["Mac", "CD14_Mono", "CD16_Mono", "cDC", "pDC", "Neutrophil",
               "Baso/Mast", "CD4T", "CD8T", "NK", "Prolif_Lym", "B_Naive",
               "B_memory", "Plasma_Cells", "RBC"],
}
CT_ORDER = [c for group in LINEAGE.values() for c in group]

# Proximal-tubule populations present in THIS object -- the CUBN / HNF4A story.
PT_TYPES = ["PT_S1", "PT_S2", "PT_S3", "iPT"]

# The only injury/disease-associated states annotated in this object.
# iPT  = injured proximal tubule; GS_Stromal = glomerulosclerosis-associated
# stroma. Deliberately short: this object does not carry the broader adaptive/
# degenerative annotation, so claiming more states than are labelled would be
# an assumption, not a finding.
INJURY_TYPES = ["iPT", "GS_Stromal"]

# Modality semantics. SC_RNA / SN_RNA measure transcript abundance. SN_ATAC rows
# in this object are gene-ACTIVITY scores (fragments over gene body + promoter),
# i.e. chromatin accessibility at the locus -- NOT mRNA. Never pool or compare
# them on an "expression" axis; they answer a different question.
RNA_TECHS = ["SC_RNA", "SN_RNA"]
ATAC_TECH = "SN_ATAC"

# ------------------------------------------------------------- palette ------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

# Fixed-order categorical slots. Never cycled: a 9th series folds to "Other".
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#4a3aa7", "#e34948"]

# Single-hue sequential (blue), light -> dark.
SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
       "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

# Diverging blue <-> red, neutral gray midpoint (never a hue at the middle).
DIV_LO, DIV_MID, DIV_HI = "#2a78d6", "#f0efec", "#d03b3b"

STATUS = {"good": "#0ca30c", "warning": "#fab219",
          "serious": "#ec835a", "critical": "#d03b3b"}

# Semantic assignments held constant across the whole figure set, so a colour
# means the same thing in every panel.
C_CONTROL = SERIES[0]     # blue
C_DISEASE = SERIES[1]     # orange
C_DKD = SERIES[1]
C_HKD = SERIES[3]
TECH_COLORS = {"SC_RNA": SERIES[0], "SN_RNA": SERIES[2], "SN_ATAC": SERIES[6]}


def cmap_seq():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("kseq", SEQ)


def cmap_div():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("kdiv", [DIV_LO, DIV_MID, DIV_HI])


# ------------------------------------------------- Nature Genetics style ----
# NG figure widths are fixed by column: 89 mm single, 120 mm 1.5-column,
# 183 mm double. Type is small (5-7 pt), rules are thin (0.5-0.75 pt), panels
# are lettered in bold lowercase at the top-left, and there are no frames.
MM = 1 / 25.4
W1, W15, W2 = 89 * MM, 120 * MM, 183 * MM   # inches

FS_TICK, FS_LABEL, FS_PANEL, FS_TITLE, FS_NOTE = 5.5, 6.5, 8, 7, 6


def use_style():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        # NG house face is Helvetica/Arial; fall back cleanly if absent.
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "Nimbus Sans",
                            "Liberation Sans", "DejaVu Sans"],
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
        "axes.linewidth": 0.5,
        "axes.titlesize": FS_TITLE,
        "axes.titleweight": "normal",
        "axes.labelsize": FS_LABEL,
        "axes.grid": False,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 2.0,
        "ytick.major.size": 2.0,
        "legend.frameon": False,
        "legend.fontsize": FS_TICK,
        "lines.linewidth": 1.0,
        "patch.linewidth": 0.4,
        "figure.dpi": 200,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        # Keep text as text in vector output so the journal can restyle it.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def panel(ax, letter, dx=-0.055, dy=1.045):
    """NG panel letter: bold lowercase, top-left, outside the axes."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=FS_PANEL,
            fontweight="bold", va="bottom", ha="left", color=INK)


def despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)


def hairline_grid(ax, axis="y"):
    """Recessive grid, always behind the marks."""
    ax.set_axisbelow(True)
    ax.grid(True, axis=axis, color=GRID, lw=0.4, zorder=0)


def check_overlap(fig, frac=0.10):
    """
    Report tick-label collisions before the figure ships.

    Renders once, then tests each adjacent pair of visible tick labels.
    Returns a list of colliding pairs; an empty list is the pass condition.

    Collision means the intersection covers more than `frac` of the smaller
    label's area -- a hair of antialiasing touching is not a defect, and for
    90-degree-rotated tick labels the naive test fires on every adjacent pair.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    hits = []

    def box(t):
        """
        Pixel extent of a Text.

        get_window_extent() returns a DEGENERATE zero-area box for axis titles
        (matplotlib repositions them during draw), which silently defeated every
        title-vs-annotation test. Fall back to get_tightbbox, which is laid out
        correctly, and only then to the window extent.
        """
        for fn in (lambda: t.get_tightbbox(r), lambda: t.get_window_extent(renderer=r)):
            try:
                b = fn()
            except Exception:
                continue
            if b is not None and b.width > 0 and b.height > 0:
                return b
        return None

    def clash(t1, b1, t2, b2):
        inter = mpl.transforms.Bbox.intersection(b1, b2)
        if inter is None:
            return False
        a = inter.width * inter.height
        smallest = min(b1.width * b1.height, b2.width * b2.height)
        return smallest > 0 and a / smallest > frac

    # Axes with the axis switched off (spacer/annotation axes, and the upper
    # partner of a sharex pair) still hand back tick-label artists that report
    # visible=True even though nothing is drawn. Testing them yields pure
    # false positives, so drop them up front.
    live = [ax for ax in fig.axes if getattr(ax, "axison", True)]

    # (1) adjacent tick labels within an axis -- the dense-axis failure mode
    for ax in live:
        for grp in (ax.get_xticklabels(), ax.get_yticklabels()):
            grp = [t for t in grp if t.get_visible() and t.get_text().strip()]
            boxes = [(t, box(t)) for t in grp]
            boxes = [(t, b) for t, b in boxes if b is not None]
            for i in range(len(boxes) - 1):
                (t1, b1), (t2, b2) = boxes[i], boxes[i + 1]
                if clash(t1, b1, t2, b2):
                    hits.append((t1.get_text(), t2.get_text()))

    # (2) every free-standing label in the figure against every other, across
    #     axes -- titles, annotations, axis labels and fig.text. This is the
    #     class of collision that dense multi-panel layouts actually produce.
    free = []
    for t in fig.texts:
        if t.get_visible() and t.get_text().strip():
            free.append(t)
    for ax in fig.axes:
        # ax.texts survive axis("off") and are genuinely drawn, so keep those
        # even for hidden axes; tick labels and axis labels are not.
        for t in list(ax.texts):
            if t.get_visible() and t.get_text().strip():
                free.append(t)
        if not getattr(ax, "axison", True):
            continue
        # ax.title is only the CENTRE title. set_title(loc="left"/"right")
        # stores its text in ax._left_title / ax._right_title, so collecting
        # ax.title alone silently skipped every left-aligned title -- which is
        # all of them in this figure set.
        # Legend entries are Text artists too, and they are placed by a rule
        # that knows nothing about annotations sitting in the same corner.
        lg = ax.get_legend()
        if lg is not None and lg.get_visible():
            for t in lg.get_texts():
                if t.get_text().strip():
                    free.append(t)
        titles = [getattr(ax, n, None)
                  for n in ("title", "_left_title", "_right_title")]
        for t in titles:
            if t is not None and t.get_visible() and t.get_text().strip():
                free.append(t)
        if not getattr(ax, "axison", True):
            continue          # tick/axis labels are not drawn; titles are
        for t in [ax.xaxis.label, ax.yaxis.label]:
            if t is not None and t.get_visible() and t.get_text().strip():
                free.append(t)
        # Tick labels are deliberately NOT added here. Within an axis they are
        # already covered by loop (1); across axes, the last tick of one panel
        # legitimately sits close to the first tick of the next, and comparing
        # them yields only false positives.
    fb = [(t, box(t)) for t in free]
    fb = [(t, b) for t, b in fb if b is not None]
    for i in range(len(fb)):
        for j in range(i + 1, len(fb)):
            (t1, b1), (t2, b2) = fb[i], fb[j]
            if clash(t1, b1, t2, b2):
                pair = (t1.get_text()[:28], t2.get_text()[:28])
                if pair not in hits:
                    hits.append(pair)

    # (3) free text sitting ON TOP OF THE DATA. This is the failure the
    #     text-vs-text test cannot see: an annotation dropped into what looked
    #     like empty space but which a bar or point cloud actually occupies.
    #     Bars are checked as filled rectangles; scatter/line marks as points.
    for ax in live:
        notes = [t for t in ax.texts
                 if t.get_visible() and t.get_text().strip()]
        # Legend entries sit on the data just as easily as a stray annotation,
        # and matplotlib's "best" placement only avoids some artists.
        _lg = ax.get_legend()
        if _lg is not None and _lg.get_visible():
            notes += [t for t in _lg.get_texts() if t.get_text().strip()]
        if not notes:
            continue
        bars = []
        for p in ax.patches:
            try:
                if p.get_visible() and p.get_window_extent().width > 0:
                    bars.append(p.get_window_extent())
            except Exception:
                pass
        pts = []
        for ln in ax.lines:
            if not ln.get_visible():
                continue
            try:
                pts.append(ax.transData.transform(ln.get_xydata()))
            except Exception:
                pass
        for t in notes:
            b = box(t)
            if b is None:
                continue
            for rb in bars:
                inter = mpl.transforms.Bbox.intersection(b, rb)
                if inter is None:
                    continue
                if inter.width * inter.height > 0.30 * (b.width * b.height):
                    hits.append((t.get_text()[:28], "<over a bar>"))
                    break
            else:
                for arr in pts:
                    if len(arr) and ((arr[:, 0] > b.x0) & (arr[:, 0] < b.x1) &
                                     (arr[:, 1] > b.y0) & (arr[:, 1] < b.y1)).sum() > 2:
                        hits.append((t.get_text()[:28], "<over plotted points>"))
                        break
    return hits


def save(fig, stem, caption=None, verify=True):
    """Write PNG + PDF + SVG, run the overlap check, and drop a caption sidecar."""
    if verify:
        hits = check_overlap(fig)
        if hits:
            print(f"  !! {stem}: {len(hits)} overlapping label pair(s): "
                  f"{hits[:4]}{' …' if len(hits) > 4 else ''}")
        else:
            print(f"  ok {stem}: no label collisions")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"{stem}.{ext}")
    if caption:
        (FIGURES / f"{stem}.txt").write_text(caption.strip() + "\n")
    plt.close(fig)
    print(f"  wrote {stem}.png/.pdf/.svg")


def tau(v):
    """Tissue-specificity index (Yanai). 0 = ubiquitous, 1 = one cell type only."""
    import numpy as np
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    if v.size < 2 or v.max() <= 0:
        return float("nan")
    return float((1 - v / v.max()).sum() / (v.size - 1))
