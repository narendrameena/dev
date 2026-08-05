# Kidney gene-panel figures — GSE211785 (Abedini et al., *Nat Genet* 2024)

Four Nature Genetics-style figures answering six questions about the 56-gene
panel from `figr_cluster_pipeline/`, built from the CKD dataset only.

```
Figures/
├── scripts/          run in numeric order; 00 must run first
├── figures/          PNG + PDF + SVG + a caption sidecar (.txt) per figure
└── source_data/      the numbers behind every panel, one CSV per analysis
    └── derived/      compact panel matrix extracted from the 17 GB h5ad
```

Run with the `biomni_e1` interpreter (the only env here with `pyBigWig`):

```bash
P=/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python
cd scripts && $P 00_extract_panel.py && for f in 01_* 02_* 03_* 04_*; do $P $f; done
```

## The figures

| Figure | Questions | What it shows |
|---|---|---|
| **1** `Figure1_expression_atlas` | Q1, Q2 | Detection per gene; specificity τ vs abundance; dot plot of all 56 genes × 38 cell types ordered down the nephron |
| **2** `Figure2_disease` | Q3, Q4 | Injured vs healthy proximal tubule; per-donor composition shifts with exact P; Disease−Control pseudobulk per cell type |
| **3** `Figure3_scRNA_vs_snRNA` | Q5 | Assay comparison within matched cell types — concordance, per-gene bias, detection, and where the two disagree |
| **4** `Figure4_regulatory` | Q6 | Promoter accessibility per cell type from snATAC, HNF4A motif scan, and the CUBN/HNF4A coupling across PT states |

## Results worth knowing

**Q1/Q2.** 52/56 genes clear a 1% detection floor. `IL1RL1`, `TDGF1`, `TWIST1`,
`CRB2` fall below it — but `CRB2` has the highest specificity in the panel
(τ = 0.99, podocyte). Low detection reflects a rare cell type, not absence, which
is why Figure 1b plots specificity against abundance. Marker sanity checks land
correctly: CRB2→Podo (0.535 vs 0.173 for the runner-up), HAVCR1→iPT (0.475 vs
0.131), UMOD→M_TAL (2.713 vs 1.818 for C_TAL), ACTA2→MyoFib/VSMC (2.755 vs
1.227), CUBN/HNF4A/AMN/ALDOB→PT_S1–S3.

**One marker does not behave as expected.** `VCAM1` peaks in **PEC** (parietal
epithelial cells, mean 0.6173, n = 598) not in **iPT** (0.6083, n = 19,190) —
a 1.48% margin, i.e. a statistical tie rather than a PEC result. The published
description of this dataset frames VCAM1 as an injured-PT marker, and iPT is
the immediate runner-up, but by strict argmax over cell types with ≥25 cells
the top slot is PEC. Do not cite VCAM1 as a clean iPT marker on the strength of
this dataset; the two populations are indistinguishable here. Independently
re-derived by two different algorithms, identical to 6 decimal places.

**Q3.** Injured PT loses its mature transport programme and gains injury and
progenitor markers — a clean dedifferentiation signature:

| down in iPT | | up in iPT | |
|---|---|---|---|
| ALDOB | −2.59 | BICC1 | +1.50 |
| CUBN | −1.70 | PROM1 (CD133) | +0.80 |
| MME | −1.23 | PAX8 | +0.70 |
| PEPD | −0.92 | VCAM1 | +0.48 |
| ACY1 | −0.82 | HAVCR1 (KIM-1) | +0.43 |

**Q4 — a negative result, stated as one.** DKD and HKD are tested separately
against Control (16 Control, 7 DKD, 13 HKD donors), each Benjamini–Hochberg
corrected within its own contrast. **Neither contrast produces a single
significant gene:** DKD 0 of 1,120 tests (smallest FDR 0.62), HKD 0 of 1,400
(smallest FDR 0.91). With 7 DKD donors this is an underpowered result as much
as a null one, and it is reported as such on the figure.

Composition carries what signal there is, and **splitting the diseases was
worth doing** — several effects are DKD-specific and were averaged away when
the two were pooled:

| cell type | Control | DKD | *P* | HKD | *P* |
|---|---|---|---|---|---|
| PT_S1 | 19.3% | 0.76% | 0.016 | 3.23% | 0.023 |
| PT_S3 | 8.26% | 0.38% | 0.021 | 1.80% | 0.069 |
| Des-Thin_Limb | 0.69% | 0.09% | **0.008** | 0.49% | 0.35 |
| Podo | 0.27% | 0.03% | 0.023 | 0.06% | 0.37 |
| B_Naive | 0.10% | 0.86% | 0.021 | 0.23% | 0.14 |
| CD16_Mono | 0.29% | 3.47% | 0.047 | 0.68% | 0.48 |

Proximal-tubule loss occurs in both diseases but is far more severe in DKD
(19.3% → 0.76% vs → 3.23%). Podocyte loss, descending-thin-limb loss and the
B-cell/CD16-monocyte expansion are nominally DKD-specific — podocyte loss being
the classical hallmark of diabetic nephropathy. Nothing survives BH correction
across 41 cell types (smallest FDR: DKD 0.19, HKD 0.85), so these are nominal
and hypothesis-generating, not established.

This differs from the published analysis, which states it treated each cell as
an independent observation. 282,610 cells from 36 donors are not 282,610
independent observations; per-cell testing inflates significance by orders of
magnitude. Everything here uses donors as the unit of replication.

**Q5.** Overall ρ = 0.54 across 1,848 gene × cell-type pairs — moderate, so the
assays are **not** interchangeable for this panel. The disagreement is
structured, not random:

- higher in **scRNA**: `JUN` (−1.35), `FOS` (−1.17), `VIM` (−0.64) — the
  classic dissociation-induced immediate-early programme. Treat scRNA values
  for these three as artefact.
- higher in **snRNA**: `BICC1` (+1.37), `PAX8` (+0.73), `CUBN` (+0.52) — long,
  intron-rich, nuclear-retained transcripts.

**Q6.** The CUBN/HNF4A argument holds. CUBN promoter accessibility is
PT-restricted (PT_S2 20.6×, PT_S3 19.9×, PT_S1 16.3×, iPT 14.8× over track
background, then a drop to 6.8×); UMOD as a positive control is correctly
TAL-restricted (M_TAL 38.8×, C_TAL 31.1×). The CUBN promoter carries three
HNF4A motif matches, best at chr10:17,129,848–17,129,861 (p = 9.1 × 10⁻⁶,
q = 0.035). And across PT states both fall together — HNF4A −0.46, CUBN −1.72
from PT_S1 to iPT — consistent with HNF4A maintaining CUBN in differentiated PT.

## What these figures do not establish

- **Motif ≠ binding.** 33/56 promoters carry an HNF4A motif at p < 10⁻⁴; the
  motif is degenerate and a 2 kb window yields chance matches. Only the
  conjunction of PT-restricted accessibility *and* a motif is informative.
  Occupancy needs ChIP-seq or footprinting.
- **Coverage ≠ peaks.** The snATAC values are bigWig coverage normalised to each
  track's genome-wide mean, not MACS2 peak calls.
- **Fold-over-track-background is not a neutral depth correction.** The
  denominator is the track's own genome-wide mean, which mixes in real signal
  from every other locus. A cell type with a sparse, sharply peaked
  accessibility landscape therefore gets inflated fold-values relative to one
  with diffuse accessibility at *identical* local signal — the metric is
  sensitive to signal concentration, not depth alone. It also reports the same
  "20×" whether hundreds or a handful of reads support it, and corrects for
  neither GC content, mappability nor copy number. It is defensible across the
  23 control tracks because they share a pipeline and their backgrounds cluster
  tightly (0.041–0.049 despite a ~58× spread in file size). It is *not*
  defensible across tracks of different construction — see the next point.
- **No disease chromatin contrast.** The portal's disease tracks are not merely
  smaller (0.1–7.5 MB vs 3.5–204 MB); they are built at 30–100× coarser
  resolution. Disease C_TAL is a *single 228,743 bp interval of value 0.0*
  spanning the whole UMOD gene and promoter, and the disease tracks carry ~19
  distinct quantised signal levels against 728 for control. Genome-wide
  `sumData` is comparable between the two, so this is the same aggregate signal
  smeared over much larger bins, not less data. The apparent disease CD4T 42.8×
  at UMOD is a 101 bp sliver of one quantised level inside an otherwise
  all-zero ~30–195 kb block — a binning artefact. This invalidates any
  promoter-scale ratio structurally. Control tracks only. A real disease
  contrast needs the raw snATAC FASTQs from SRA (GEO hosts no processed ATAC
  for GSE211785 — all 17 samples have `supplementary_file = NONE`).
- **snATAC gene-activity is not expression.** The `SN_ATAC` rows in the h5ad are
  fragments over gene body + promoter. They are excluded from every expression
  panel; pooling them would be a category error.
- **`DCN` peaks in podocytes** (0.83 vs 0.56 in Fibroblast_1), which is
  biologically odd for a stromal proteoglycan. n = 2,326 podocytes, so not a
  small-sample artefact — most likely ambient RNA from interstitium or stromal
  doublets in that cluster. Flagged rather than silently reported.

## Corrections log

Every number in this figure set was independently re-derived from the raw
sources by three separate verification passes that did not read the derived
files. Three defects were found and fixed; they are recorded here rather than
quietly patched, because two of them changed published numbers.

**1. A donor was counted in both arms of the composition test.**
The test de-duplicated donors on `(sample, status)`, which returns **37 rows for
36 donors**: donor `HK2770` appears once as Control and once as Disease, so it
was compared against itself. Every composition *P*-value in the first version of
Figure 2b was invalid. Root cause is that `obs['Status']` is internally
inconsistent — HK2770's 3,635 scRNA cells are labelled `Disease` while its 4,233
snRNA cells are labelled `Control`. De-duplication is now on `sample` alone and
grouping comes from `obs['group']`, which is consistent for all 36 donors
(asserted at runtime). Two further donors disagree between the columns:
`HK2891` (group HKD, Status Control) and `HK2663` (group Control, Status
Disease).

**2. The pseudobulk library-size offset was not a library size.**
The code used `obs['nCount_RNA']` and a comment asserted it was "the summed
total UMI across all 34,733 genes". It is not: it is integer-valued for only
**5.6%** of RNA cells, equals the true row-sum for 2 cells in 282,610, and
under-counts the true total by ~2% for scRNA but ~11% for snRNA — a
modality-dependent bias applied to a modality-mixed pseudobulk. `00` now streams
the counts layer to compute the true per-cell library size (verified 100%
integer) and every pseudobulk uses it. Independent testing showed this did not
overturn the "nothing survives FDR" conclusion, but it did move individual
*P*-values by up to 0.22 absolute.

**3. VCAM1 was described as an iPT marker. It is not, in this dataset.**
Asserted from reading the heatmap rather than from the computed value. By strict
argmax over cell types with ≥25 cells, VCAM1 peaks in **PEC** (0.6173, n = 598)
above **iPT** (0.6083, n = 19,190) — a 1.48% margin, i.e. a tie, independently
re-derived by two different algorithms to 6 decimal places. The figures always
computed this correctly; only the prose was wrong. See the Q1/Q2 section.

**Also revised:** the disease-track exclusion in Figure 4 was originally
justified by file size. The real defect is resolution — see "What these figures
do not establish". The original justification was correct in conclusion but weak
in evidence.

**Scope of the fixes.** Only Figure 2 uses donor labels (`01`, `03` and `04` use
`celltype` and `tech` only), so corrections 1 and 2 are confined to it. Figures
1, 3 and 4 passed verification unchanged, several values to 4+ significant
figures.

## Conventions

- **Cell-type labels are read from `obs['Cluster_Idents']`**, not assumed. The
  KPMP adaptive/degenerative vocabulary (aPT, aTAL1/2, dPT, aFIB, cycEPI)
  belongs to the *PostSCVI KPMP-merged* object and does **not** exist here.
  Use the PreSCVI object: the PostSCVI files are cut to 3,000 HVGs and lack HNF4A.
- **Colour is semantic and constant** across figures: blue = control, orange =
  disease, and the categorical slots follow a fixed CVD-safe order. Sequential
  scales are single-hue; diverging scales are blue↔red with a neutral grey
  midpoint.
- **Expression heatmaps are dot plots** — size = % of cells expressing, colour =
  z-score across cell types. Scaling each row to its own maximum was rejected:
  it makes a gene averaging 0.05 look like one averaging 5.0.
- **Every figure is checked for label collisions before it is written.**
  `kidneyviz.check_overlap()` tests adjacent tick labels, all free-standing text
  pairwise, and text sitting on top of bars or point clouds. Scripts print
  `ok <figure>` or list the offending pairs.

## Data sources

| What | Where |
|---|---|
| sc/sn/ATAC expression | `data/Kidney_SC/GSE211785_Susztak_SC_SN_ATAC_merged_PreSCVI_final.h5ad` (338,565 cells × 34,733 genes) |
| snATAC coverage | `data/Kidney_scATAC/{control,disease}/*.bw` — 46 tracks, 23 cell types, from `s3.us-east-2.amazonaws.com/ksusztak.genemap/HKSI/` via susztaklab.com/Human_CKD_snATAC |
| motifs | `resources/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme` |
| gene models | `source_data/hg38.refGene.txt.gz` — **not in the repo**, fetch with `curl -o Figures/source_data/hg38.refGene.txt.gz https://susztaklab.com/Human_CKD_snATAC/Annotation/hg38.refGene.txt.gz` |
| promoter sequence | UCSC REST API, hg38 |
