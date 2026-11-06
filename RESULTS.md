# ROD/CADO results summary

Every number and figure the manuscript quotes should come from this page. If the
manuscript and this page disagree, this page wins. If this page and `Analysis.Rmd`
disagree, the Rmd wins, so fix this page.

**Last checked:** 2026-11-06, against commit `5b40739` and the tracked outputs in
`analysis/` and `data/`.

**Status key**
- **Verified:** recomputed from a tracked repo file on the date above.
- **Knit:** printed by `Analysis.Rmd` in the last knit (`Analysis.md`, 2026-02-25).
- **Server:** produced on the server (Oct 2026) and not yet in a committed knit.
  Confirm at the next full re-knit.

---

## 1. Design

| Item | Value | Status |
|---|---|---|
| Animals | Juvenile *C. virginica* from one mass spawn | — |
| Treatments | CONCON, STRCON, CONROD, STRROD (STR = cycling CADO; ROD = *A. crassostreae* challenge) | — |
| Survival trial | 90 animals per treatment, followed to day 22 | Verified (`data/oyster_survival_data.csv`) |
| Sequenced | 80 animals, 20 per treatment, sampled at day 22 | Knit (PERMANOVA df: 3 + 76) |

## 2. Mortality (Figure 1A, 1C)

| Treatment | Deaths / N | Mortality ± SE |
|---|---|---|
| CONCON | 12 / 90 | 13.3 ± 3.6 % |
| STRCON | 19 / 90 | 21.1 ± 4.3 % |
| CONROD | 25 / 90 | 27.8 ± 4.7 % |
| STRROD | 59 / 90 | 65.6 ± 5.0 % |

- **Log-rank test across all four curves:** χ² = 76.9, df = 3, p = 1.4 × 10⁻¹⁶.
  The Rmd prints this p on the KM plot (`add_pvalue`). Status: Verified.
- **Pairwise log-rank, unadjusted:** these are *not in the Rmd yet*; add them before
  quoting.
  - STRROD vs CONCON: p = 1.4 × 10⁻¹³
  - STRROD vs STRCON: p = 3.4 × 10⁻¹⁰
  - STRROD vs CONROD: p = 3.7 × 10⁻⁷
  - CONROD vs CONCON: p = 0.014
  - STRCON vs CONCON: p = 0.17
  - CONROD vs STRCON: p = 0.25

## 3. SNP data set

| Item | Value | Status |
|---|---|---|
| Filtered VCF | `SNP.TRSdp10g1.FIL.vcf.gz`: minDP 10, no missing data, MAF 0.01, biallelic for plink/LFMM | Rmd |
| SNP count after filtering | **TBD.** Fill this in at the re-knit; there's no tracked output. Older docs quote 22,874, which is likely a different filter set, so don't use it. | Pending |
| Diversity windows | 10 kb non-overlapping, mitochondrion excluded | Rmd |

## 4. Genome-wide diversity (Figure 1D; 10 kb windows)

These are one-way ANOVAs with windows as replicates, followed by Tukey HSD. With
about 40,000 windows per treatment, almost every contrast is significant. Report
the effect sizes, not just the p-values. Status for this whole section: Verified.

| Treatment | π (mean ± SE) | H<sub>o</sub> | H<sub>e</sub> | Tajima's *D* |
|---|---|---|---|---|
| CONCON | 0.00800 ± 0.00003 | 0.233 ± 0.0004 | 0.197 | 0.003 ± 0.003 |
| STRCON | 0.00777 ± 0.00003 | 0.224 ± 0.0004 | 0.192 | 0.021 ± 0.003 |
| CONROD | 0.00844 ± 0.00003 | 0.262 ± 0.0004 | 0.207 | 0.131 ± 0.003 |
| STRROD | 0.00773 ± 0.00003 | 0.216 ± 0.0004 | 0.191 | −0.031 ± 0.003 |
| ANOVA | F₃,₁₆₃₃₄₃ = 113, p = 3 × 10⁻⁷³ | F₃,₁₆₃₅₆₈ = 2895, p < 10⁻³⁰⁰ | F = 747, p < 10⁻³⁰⁰ | F₃,₁₆₂₉₁₃ = 471, p = 1 × 10⁻³⁰⁴ |
| Tukey | All pairs differ except STRCON–STRROD (p = 0.80) | All pairs p < 0.001 | All pairs p < 0.001 | All pairs p < 0.001 |

- **Pattern:** CONROD has the highest diversity on every metric, and STRROD the lowest.
- **H<sub>o</sub> exceeds H<sub>e</sub> in every treatment.** This is consistent with
  negative F<sub>IS</sub> genome-wide.
- **Stale numbers to replace:** the older Results doc gives Tajima's *D* as about −25.
  That doesn't match these outputs. Use the values above.
- **Known quirk:** `treat.10kb.txt` is built with `NR>1` on headerless bed files, which
  drops each treatment's first window. These numbers reproduce that quirk (40,893
  windows per treatment). Fixing it changes the means by less than 0.0001.

## 5. Outlier detection

| Item | Value | File | Status |
|---|---|---|---|
| LFMM2 settings | `lfmm_ridge`, K = 4; GIF fixed at 0.98 (CADO, disease), calculated for combined; q < 0.1 | Rmd | Rmd |
| LFMM2 outliers, CADO | 5 | `analysis/outlier_detection/results/outliers_lfmm2_CADO_098.loci.txt` | Verified |
| LFMM2 outliers, disease | 3 | `…_DIS_098.loci.txt` | Verified |
| LFMM2 outliers, combined | 107 | `…_CB.loci.txt` | Verified |
| LFMM2 outliers, total | 115 | `outliers_lfmm2.bed` | Verified |
| Combined-model GIF | 1.41 (calculated) | Outline | **Confirm at re-knit** |
| OutFLANK outliers (4 treatments as pops, q < 0.1) | 126 | `analysis/FST/Outflank_outliers.bed` | Verified |
| Shared LFMM2 ∩ OutFLANK | 91 | `outliers_overlap_lfmm2_outflank.bed` | Verified |
| Union LFMM2 ∪ OutFLANK | 150 SNPs (115 + 126 − 91); 148 merged BED intervals | `outliers_lfmm2_outflank_combined.bed` | Verified |
| LFMM-only SNPs with no FST | 14 (FST set filtered at MAF 0.05, LFMM at MAF 0.01) | — | Server |

**Do not use "230 outliers (115 + 115)".** That figure predates the OutFLANK rerun and
double-counts the overlap.

### Outlier PCA + PERMANOVA (adonis2, 999 permutations, 80 animals)

| Outlier set | R² | F₃,₇₆ | p | Status |
|---|---|---|---|---|
| All LFMM2 (115) | 0.544 | 30.2 | 0.001 | Knit |
| OutFLANK (126) | 0.525 | 28.0 | 0.001 | Knit |
| Union | 0.527 | 28.2 | 0.001 | Knit |
| Overlap (91) | 0.550 | 30.9 | 0.001 | Knit |

**Caveat:** these sets were chosen for treatment association, so separation by
treatment is expected. Treat these tests as descriptive, not as independent evidence.

## 6. Figure 2: LFMM Manhattan + focal regions

- **File:** `figures/Figure2.png`.
- **Panel A:** the combined-model LFMM, truncated at p ≤ 0.1 (`man.pmax`).
- **Panel B:** π, H<sub>o</sub> and Tajima's *D* in 500 bp bins, ±20 kb around each
  focal SNP.
- **Focal SNPs:** 9 SNPs in 7 regions, taken from `analysis/figure2/focal_loci.csv`
  (picked with `exploratory/separation_score.Rmd`):
  - DOCK3-like: chr 8, 62,714,735
  - semaphorin receptor: chr 8, 60,756,955
  - PP2A regulatory subunit: chr 4, 12,033,923
  - GGT1-like: chr 5, 18,844,269 and 18,844,273
  - MsrB-like: chr 6, 32,344,416
  - titin / SRRM2: chr 2, 19,679,570 and 19,679,574
  - UNCH / mucin-like: chr 1, 65,412,698
- **Gene product names** were transcribed from the Discussion draft. They have not
  been checked against NCBI.
- **Truncation hides the null cloud.** Add a supplementary QQ plot (open task).

## 7. Figure 3: meadow plots (pairwise windowed FST × H<sub>o</sub> ratio)

All of section 7 is Server status.

- **Method:** 10 kb windows, each treatment vs CONCON, Weir & Cockerham FST, windows
  with fewer than 10 SNPs dropped. Output: `figures/Figure3.png` (not yet committed).
- **Main selection signal: two narrow peaks on chr 8 (NC_035787.1).**
  - **60.69–60.78 Mb (semaphorin receptor):** STRROD FST up to about 0.18.
  - **62.68–62.84 Mb (DOCK3-like):** STRROD FST up to about 0.20; STRROD H<sub>e</sub>
    falls from 0.26 to 0.12.
  - At both peaks CONROD is weaker and STRCON shows nothing.
- **Read-depth check:** no CNV or paralog signal (all |log2 depth ratio| < 0.05).
  Output: `figures/Figure.S.meadow.depth.png`.

## 8. Supplement: relatedness and the chr 7 block

All of section 8 is Server status.

- **Within-treatment kinship (KING, permutation):** p = 0.40, so relatedness is even
  across treatments. Mean kinship is about 0.07.
- **chr 7 (NC_035786.1), 0.4–14.8 Mb:** this is a long parental haplotype, not an
  inversion.
  - Carriers per treatment (of 20): CONCON 5, STRCON 8, CONROD 13, STRROD 9.
  - CONROD vs the others: p = 0.038 (post hoc). All four treatments: p = 0.10.
  - Supplement only.
- **Output:** `figures/Figure.S.family.chr7.png` (not yet committed).

## 9. Gene Ontology (union outliers)

| Item | Value | File | Status |
|---|---|---|---|
| Candidate genes containing outliers | 97 | `analysis/GO/Sig.gene.ROD.CADO.LOC` | Verified |
| Enriched GO terms (topGO, Fisher p < 0.05) | 68: 29 BP, 13 CC, 26 MF | `analysis/GO/ROD_OUT_GO_en_sig_gene.csv` | Verified |
| Genes carrying an enriched term | 48 | same | Verified |

**Do not use "144 genes / 110 GO terms (47/15/48)".** Those numbers predate the
corrected LFMM run.

**Terms that are no longer supported:**
- bicarbonate transport and Na⁺:HCO₃⁻ symporter
- methylcytosine dioxygenase and 5-mC catabolism
- deubiquitinase terms (GO:0004843, GO:1990380)
- p53 binding
- phosphatidylinositol dephosphorylation
- regulation of GTPase activity

The section-by-section rewrite guidance is in "Changes that NEED TO BE MADE" on
SharePoint.

## 10. Figures

| Figure | File | In repo? |
|---|---|---|
| 1: mortality, PCA, KM, π/H<sub>o</sub> | `figures/Figure1.png` | Yes |
| 2: LFMM + focal regions | `figures/Figure2.png` | Yes |
| 3: meadow plots | `figures/Figure3.png` | **No.** Commit after the re-knit |
| S: FST Manhattan | `figures/Figure.X.FST.man.png` | Yes |
| S: depth check | `figures/Figure.S.meadow.depth.png` | No |
| S: kinship + chr 7 | `figures/Figure.S.family.chr7.png` | No |
| S: GO enrichment (union) | `figures/GO_enrichment_plot.png` | Yes |
| S: GO overlap / comparison | `figures/GO_enrichment_overlap_plot.png`, `…_comparison_plot.png` | No |
| S: LFMM QQ plot | — | Not built yet |

## Open items before the numbers are final

- [ ] Do a full re-knit on the server, then commit `Analysis.md` and the figures
      marked "No" above.
- [ ] Record the filtered SNP count and the combined-model GIF.
- [ ] Add the pairwise log-rank tests to the Mortality chunk.
- [ ] Add the LFMM QQ plot.
- [ ] Check the Figure 2 gene names against NCBI.
