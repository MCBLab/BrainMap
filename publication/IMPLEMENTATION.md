# VEGABrain manuscript — implementation notes

How the *Neuroinformatics* submission in this folder is put together: what each
file is for, where every number in the text comes from, what is still a
placeholder, and how to build the PDF.

- **Target journal:** *Neuroinformatics* (Springer), <https://link.springer.com/journal/12021>
- **Article type:** **Software Original Article** — the journal's type for a paper
  describing a software product. It asks for more than a description: a short
  account of how to install and run the tool, a validation of its correctness
  against accepted standards, and an evaluation against comparable packages
  including benchmarking where applicable. Section 4 of `main.tex` exists
  specifically to satisfy that requirement (`Kennedy2009` in the bibliography is
  the journal's own statement of the bar for software papers).
- **Working title:** VEGABrain: a Visualizer of Enriched Gene Analysis for the
  developing human brain
- **Status:** full text drafted; figures 2–5 are placeholders and the measured
  values in Section 4 are marked `[TODO]`. See [Second pass](#second-pass).

---

## 1. Files

```
publication/
├── main.tex                  the manuscript
├── main.pdf                  built output (17 pp, A4)
├── references.bib            34 entries
├── build.sh                  build / clean
├── sn-jnl.cls                Springer Nature class (copied from template/)
├── sn-basic.bst              bibliography style used by main.tex
├── bst/                      the other eight SN styles, unused
├── template/                 pristine unpacked template, incl. user-manual.pdf
├── figures/
│   ├── make_placeholders.R   generates the placeholder PDFs currently in use
│   ├── make_figures.R        generates the real figures (not yet run)
│   └── fig{2,3,4,5}*.pdf     placeholders
├── references/               20 PDFs of cited works
└── IMPLEMENTATION.md         this file
```

The template package is the official Springer Nature journal article template,
v3.1 (December 2024), from
<https://www.springernature.com/gp/authors/campaigns/latex-author-support>.
`template/` holds the untouched extraction, including `user-manual.pdf` (the
authoritative reference for the class) and `sn-article.pdf` (a rendered
example). `sn-jnl.cls` and `sn-basic.bst` are copied up one level because BibTeX
does not search subdirectories.

**Document class:** `\documentclass[pdflatex,sn-basic]{sn-jnl}`. SN Basic is the
author–year reference style *Neuroinformatics* uses. Add the `lineno` option for
the reviewers' copy:

```latex
\documentclass[lineno,pdflatex,sn-basic]{sn-jnl}
```

`\usepackage[T1]{fontenc}` is deliberately **absent**. This TinyTeX has no Type 1
EC fonts, so T1 triggers Metafont bitmap generation, which bloats the PDF and
breaks `latexmk`'s error detection. The SN template leaves fontenc off for the
same reason; UTF-8 accents in author names compose correctly without it.

## 2. Building

```bash
bash build.sh          # -> main.pdf, keeps main.log / main.blg
bash build.sh clean    # remove LaTeX intermediates
```

LaTeX comes from TinyTeX in `~/.TinyTeX` (TeX Live 2026), installed with
`Rscript -e 'tinytex::install_tinytex()'`. The build goes through
`tinytex::latexmk()` rather than calling `latexmk` directly because it
auto-installs any LaTeX package the template pulls in.

Verifying a build:

```bash
grep -c 'Citation .* undefined'  main.log   # expect 0
grep -c 'Reference .* undefined' main.log   # expect 0
grep -c '^!'                     main.log   # expect 0
grep -c 'Overfull \\hbox'        main.log   # expect 0
grep -oP "You've used \K\d+"     main.blg   # expect 34, = entries in references.bib
```

Last verified build: 17 pages, A4, 0 undefined citations, 0 undefined
references, 0 errors, 0 overfull boxes, 34/34 bibliography entries cited.

## 3. Where the numbers come from

Every quantitative claim in `main.tex` was read off the code or the data in this
repository, not estimated. Keep this table in sync if the pipeline changes.

| Claim in the text | Value | Verified from |
|---|---|---|
| Expression source | BrainSpan RNA-Seq, Gencode v10 RPKM, gene-summarised | `ETL.qmd`, `genes_matrix_csv/` |
| Matrix dimensions | 52,376 genes × 524 samples | `wc -l genes_matrix_csv/{rows,columns}_metadata.csv` |
| Transformation | `log2(RPKM + 1)`, duplicate symbols dropped | `preparaDados.R` |
| Developmental windows | 5; n = 5, 10, 5, 8, 14 donors | `preparaDados.R` (`age_mapping`) |
| Atlas regions (micro) | 42 | `length(unique(ontologia_micro$region))` |
| Anatomical groups (macro) | 6 | `mapeamento_regioes.csv`, column `macro_regions` |
| Views per map | 4 (DK lateral, DK medial, aseg sagittal, aseg coronal_1) | `plumber.R:build_brain_grid` |
| Gene sets scored | 18,650 | `length(readRDS("genesets_list.rds"))` |
| — GOBP / HP / GOMF / REACTOME / GOCC / BIOCARTA / KEGG / HALLMARK | 7,538 / 5,793 / 1,872 / 1,839 / 1,080 / 292 / 186 / 50 | prefix table over `names(genesets_list)` |
| Aggregated score rows | 3,188,466 | `nrow(readRDS("ontologia_micro.rds"))` |
| Sample-level scores before aggregation | ≈9.8 M | 18,650 × 524 |
| Cloud Run config | 4 GiB, 2 vCPU, concurrency 4, timeout 900 s, min-instances 0 | `deploy.sh` |
| Cold start | 24.8 s | timed `GET /list_ontologies` against the live API |
| Live front end | <https://mcblab.github.io/BrainMap/> | `.github/workflows/pages.yml` |
| Live API | <https://brainmap-api-5kofy56taq-rj.a.run.app> | extracted from the deployed JS bundle |

Regenerate the gene-set breakdown with:

```r
gs <- readRDS("genesets_list.rds")
sort(table(sub("_.*", "", names(gs))), decreasing = TRUE)
```

## 4. Section plan

| § | Section | ~words | Content source |
|---|---|---|---|
| — | Abstract | 200 | — |
| 1 | Introduction | 850 | The three gaps: portals are non-anatomical and gene-at-a-time; programmatic tools assume coding and local data; nothing takes a *gene set* to a developmental map. Contributions listed explicitly. |
| 2.1 | Expression data | 150 | `preparaDados.R`, `ETL.qmd` |
| 2.2 | Developmental windows | 120 | `preparaDados.R` |
| 2.3 | Atlas mapping | 300 | `mapeamento_regioes.csv`, `plumber.R` |
| 2.4 | Composite enrichment scores | 450 | `plumber.R`, `preparaDados.R`; Eq. 1 from `Barbie2009` |
| 2.5 | Implementation & deployment | 300 | `plumber.R`, `Dockerfile.api`, `deploy.sh`, `pages.yml` |
| 3 | Results / tool description | 700 | `front/src/componentes/Home.jsx`, `plumber.R` endpoints |
| 4 | Validation & comparison | 700 | **partly TODO** — see below |
| 5 | Discussion & limitations | 700 | Donor n, bulk tissue, mapping approximation, rank-based scores, cold start |
| 6 | Conclusion | 130 | — |
| — | Back matter | 250 | Declarations; several `[TODO]` |

## 5. Figures and tables

| Item | Label | State | Regenerate with |
|---|---|---|---|
| Fig. 1 Architecture | `fig:arch` | **real** — TikZ, inline in `main.tex` | rebuild |
| Fig. 2 Interface | `fig:ui` | placeholder | screen capture of the live app |
| Fig. 3 Single gene (*SOX10*) | `fig:gene` | placeholder | `Rscript figures/make_figures.R fig3` |
| Fig. 4 Ontology composite score | `fig:ontology` | placeholder | `Rscript figures/make_figures.R fig4` |
| Fig. 5 User gene list | `fig:genelist` | placeholder | `Rscript figures/make_figures.R fig5` |
| Table 1 Data summary | `tab:data` | **real** | — |
| Table 2 API endpoints | `tab:api` | **real** | `plumber.R` |
| Table 3 Response times | `tab:perf` | cold start real, warm `[TODO]` | benchmark script (to write) |
| Table 4 Tool comparison | `tab:compare` | **real** | — |

`figures/make_figures.R` is written and parses, but has **not been run**. It
reproduces `plumber.R:build_brain_grid` exactly — the same four views, the same
facetting, the same colour ramp — so the figures in the paper are the figures the
tool serves. It writes PDF directly rather than converting the API's SVG, because
no `rsvg-convert`/`inkscape`/ImageMagick is installed here. It needs
`dados_otimizados.rds`, `ontologia_micro.rds` and `ontologia_macro.rds` in the
repository root.

`figures/make_placeholders.R` regenerates the four placeholder PDFs currently
referenced by `main.tex`. Filenames are identical to the real outputs, so
swapping them in requires no edit to `main.tex`.

## 6. Per-structure enrichment testing

The maps paint magnitude. `figures/enrichment_stats.R` adds inference: which
structures are *significantly* enriched for a gene set, in which window.

### Unit of inference

**The BrainSpan structure, never the ggseg parcel.** The atlas draws 29 parcels
but BrainSpan sampled only 22 structures, and several parcels are byte-identical
(`Caudate`/`Putamen`/`accumbens area`/`insula` are all one "striatum"
dissection). Testing parcels would report one measurement as up to four
independent findings. Results are painted outward for display only, and a parcel
inherits its structure's verdict — scores can be averaged across contributing
structures, p-values cannot.

That gives **87 populated (structure, window) cells, 70 with >= 3 donors**
(`cells_from_meta()`, straight from `dados_otimizados.rds`). Cells below that
are not tested (`tested = FALSE`) and render grey. `p.adjust` is fed `NA` for
them, so the BH family is 70, not 87.

### One competitive test and a donor guard

| | asks | null |
|---|---|---|
| competitive | in this cell, is the set's deviation from the window mean higher than random genes would give? | 4,000 random gene sets of the same size, each query gene swapped for one from the same mean-expression bin |
| guard | does it survive dropping any single donor? | the same statistic, recomputed without each donor in turn |

The competitive p is a rank against the null sets in that cell, BH over the ~70
tested cells, q < 0.05. A cell that passes is then re-ranked with each of its
donors left out (query and null alike, window mean adjusted); it is `enriched`
only if every leave-one-out p stays below 0.05 (`p_loo` reports the worst one).
The guard runs at nominal alpha on cells that already passed BH; it is not part
of the multiplicity family.

**Scores without GSVA.** GSVA's ssGSEA statistic is a sum of three per-gene
terms (`GSVA:::.fastRndWalk`): with per-sample ranks R, weights R^0.25 and
C = n − position + 1, `ES = Σ_S R^α·C / Σ_S R^α − (n(n+1)/2 − Σ_S C)/(n − |S|)`.
`rank_tables()` computes the three tables once (~6 s, ~300 MB) and
`ssgsea_from_ranks()` scores any number of sets as sparse-matrix products.
Check 0 of `validate_stats.R` asserts identity with `GSVA::gsva(normalize =
FALSE)` (relative error ≤ 1e-5). A user's list therefore costs ~1 s, not the
~45 s of a GSVA call, and there is no precomputed reference to build or ship.

### Why the earlier two-test design was replaced

The first version ranked the query against the 4,000 nearest-size **MSigDB
sets** and gated it with a **within-donor label permutation**. Diagnostics on the
real data (200 sampled MSigDB sets) showed both halves answering the wrong
question:

- Ranks within a population of curated sets are uniform by construction, so at
  most ~1% of sets can be called per cell whatever their biology, and a set had
  to be more regionally specific than 99.9% of MSigDB — itself brain-structured.
  Only **4 of 200** real sets produced any enriched cell.
- The permutation compared a structure with its donor's average over the other
  structures, which is dominated by cortex (11 of 16). Subcortical tissue differs
  from that average for every gene set: it passed **32% of subcortical cells and
  4% of cortical ones** on arbitrary sets. It filtered cortex out instead of
  checking donor consistency, and pooled random-gene nulls inherit the same
  bias (subcortical cells have heavier tails). Only a per-cell rank against
  size- and expression-matched random sets is tissue-neutral.

### Three things that are load-bearing

**Within-window centring.** Centre globally instead and any gene set with a
developmental trend lights up every structure of its peak window — a true
statement about time, a false one about place.

**Expression matching of the null.** Roughly half the matrix is barely
detected; a null drawn uniformly from it is noise, narrower in every cell than
any real signature, and calls almost anything enriched. Each null gene comes
from the query gene's own mean-expression bin (50 bins).

**The leave-one-donor-out guard.** Cells with 3–5 donors (1st and 3rd
trimester) can pass BH on one donor's tissue; the guard removes exactly those
(HALLMARK_OXIDATIVE_PHOSPHORYLATION: see table).

### Calibration, power and controls (`validate_stats.R`)

Run with `Rscript figures/validate_stats.R 300` (k = 4,000). Null lists and
null sets are drawn from the 17,445-symbol MSigDB universe.

| check | result | nominal |
|---|---|---|
| score identity, `ssgsea_from_ranks` vs `GSVA::gsva` | max rel. diff 1e-5 (6 sets, sizes 15–200) | 0 |
| cell-level, 300 random lists, size 15 / 50 / 200 | p ≤ 0.05: **0.051 / 0.046 / 0.048**; p ≤ 0.01: 0.012 / 0.010 / 0.010 | 0.05 / 0.01 |
| family level (≥ 1 asterisk on a null list), size 15 / 50 / 200 | **0.043 / 0.027 / 0.027** | ≤ 0.05 |
| weighted path, 100 random lists × N(0,1) weights | cell 0.043 at p ≤ 0.05; family 0.020 | 0.05 / ≤ 0.05 |
| power, 200 real MSigDB sets (size 15–300) | **64 / 200 (32%)** with ≥ 1 enriched cell; median 4 cells per hit set | — |
| donor guard | 384 cells pass BH, 308 survive (removes 20%) | — |
| hits by class (tested cells 50 cortical : 20 subcortical) | 118 cortical : 190 subcortical | — |
| `HALLMARK_OXIDATIVE_PHOSPHORYLATION` | 12 cells pass BH, 5 survive the guard; the 7 removed are 3–5-donor 1st/3rd-trimester cortical cells with p_loo 0.4–1.0 | — |

Controls (cells, structures):

| gene set | cells | where |
|---|---|---|
| `GOBP_CEREBELLAR_CORTEX_DEVELOPMENT` | 4 | cerebellar cortex ×4 (2nd, 3rd trimester, Infant, Adult) |
| `GOBP_STRIATUM_DEVELOPMENT` | 4 | striatum ×4 |
| `GOBP_THALAMUS_DEVELOPMENT` | 7 | mediodorsal thalamus ×3, dorsal thalamus ×1, striatum ×2, amygdala ×1 |
| `GOBP_HIPPOCAMPUS_DEVELOPMENT` | 6 | hippocampus ×2, amygdala ×2, anterior cingulate ×2 |
| `GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION` | 35 | 27 cortical cells across 9 areas, plus amygdala, hippocampus, striatum |
| `GOBP_FOREBRAIN_GENERATION_OF_NEURONS` | 19 | amygdala ×4, hippocampus ×4, 10 cortical cells |
| `GOBP_MYELIN_ASSEMBLY` | 1 | mediodorsal thalamus (2nd trimester) |

Under the previous design the same 200-set sample gave 4 sets with any hit,
and the cerebellar control was the only anatomically named set that found its
own structure.

The honest limitations: the test is calibrated per cell, but hits are not
evenly spread over cells. Hippocampus (77 of 308), amygdala (49) and cerebellar
cortex (39) — the tissues most distinct from the cortical mean — each collect
more hits than any cortical area (7–15), so a subcortical cell is called about
four times as often as a cortical one. That is a property of the tissue under a
within-window contrast, not of the null (a random list is equally likely to be
called anywhere), and it is printed by `validate_stats.R` so it is never a
surprise. The p floor is
1/(k+1) = 2.5e-4 at k = 4,000, so a lone hit must beat all but two null sets to
survive BH over 70 cells; Fig. 6 uses k = 10,000.

### Power is uneven and must be shown

Donors per window are 5 / 10 / 5 / 14 / 8. The **3rd trimester** is the weak one:
its donors contribute 1, 1, 3, 15, 15 samples and its median structure rests on
two people. It is flagged (`underpowered_window`), not suppressed. The 1st
trimester is *not* weak despite also having 5 donors — 67 samples across 21
structures.

### In the app

`plumber.R` sources `figures/enrichment_stats.R`; cells come from
`dados_otimizados.rds` at boot, and the rank tables and the z-score matrix for
weighted signatures are built on first use (not at boot: Cloud Run runs at
min-instances = 0 and times the cold start). The Gene List endpoints
(`/plot_genelist`, `/plot_genelist_png`, `/data_genelist`) take
`significancia=true`; the front-end sends it from a checkbox that is on by
default. Enriched parcels get the same white-on-dark asterisk as Fig. 6 and the
caption stamped on the map says what it means. The CSV gains `tested / enriched
/ q_competitive / p_loo` per region, painted outward with the same rule as the
map. Fewer than 5 genes found: no test, and the caption says so.

The unweighted path scores the list with `ssgsea_from_ranks()` (identical to
GSVA, so the displayed map is unchanged) and divides by the set's own range for
display, exactly what GSVA's normalisation does. **A log2FC-weighted signature
is not an ssGSEA score**; `enrich_weighted()` uses the same null family in the
weighted score's own currency (random expression-matched sets carrying the
query's weights, z-scored expression) and the same guard.

### Files and runtimes

```
figures/enrichment_stats.R    statistics; pure, no I/O
figures/validate_stats.R      the table above, ~15 min
figures/fig6_enrichment.R     -> fig6_enrichment.pdf, ~15 s
```

A gene list costs ~1 s warm; the first request after a cold start adds ~6 s for
the rank tables.

### Rendering

Fill is the **within-window deviation**, not the raw score. Painting the raw
score produces a figure that argues with itself: the cerebellum carries its
highest absolute ssGSEA prenatally but is only *distinctive* among structures in
infancy and adulthood, so a raw-score fill shows a blue parcel starred as
significant.

Significance is an asterisk drawn on each enriched parcel (`geom_text`, white
glyph over a thin dark halo so it reads on dark and light fills alike). The point is
the area-weighted centroid of the parcel's largest polygon, taken from ggseg's
own repositioned polygons for that view so it lands where the parcel is drawn.
Where the two hemispheres are stacked on top of each other (aseg's sagittal
cerebellum) they collapse to one asterisk. The text layer sets no `fill`, so it
does not trigger the discrete `scale_fill_manual` ggseg appends to fill-less
brain layers, which is what made an outline overlay unusable.

ggseg >= 2.2 no longer inherits the plot's data: `geom_brain` is given the data
explicitly, grouped by `broad_age`, so each facet receives a full atlas.

## 7. Reference ledger

34 entries in `references.bib`; 20 have the full text archived in
`references/<Key>.pdf`. Bib key ⇄ filename is 1:1 (`Zhu2019.pdf` was renamed to
`Zhu2018.pdf` and `mowinckel2020.pdf` to `Mowinckel2020.pdf` to keep it so).

**Archived locally (20):** `Chen2013`, `Chopra2023`, `GOConsortium2023`,
`Gillespie2022`, `Hanzelmann2013`, `Kanehisa2000`, `Kohler2021`, `Kuleshov2016`,
`Lau2008`, `Li2018`, `Liberzon2011`, `Liberzon2015`, `Markello2021`,
`Markello2022`, `Mowinckel2020`, `Subramanian2005`, `Sunkin2013`,
`Wilkinson2016`, `Yi2020`, `Zhu2018`.

**No local PDF (14)** — cited on verified Crossref metadata; no open-access copy
was retrievable:

| Key | Why no PDF |
|---|---|
| `Kang2011`, `Miller2014`, `Barbie2009`, `Ashburner2000` | PMC author manuscripts; PMC and Europe PMC both refuse automated PDF retrieval |
| `Hawrylycz2012`, `Desikan2006`, `Fischl2002`, `Kennedy2009`, `Fleck2021` | paywalled, no OA copy |
| `BrainSpan2011` | web resource, not an article |
| `RCoreTeam2026`, `Schloerke2024plumber`, `Chang2025shiny`, `Dolgalev2025msigdbr` | software citations |

Bibliographic metadata for every DOI-registered entry came from the Crossref
REST API, so titles, journals, volumes and pages are publisher records rather
than hand-typed. Author lists longer than 12 names were truncated with
`and others` (14 entries; the consortium papers would otherwise run to hundreds
of names).

To add a reference:

```bash
curl -s "https://api.crossref.org/works/<DOI>?mailto=<you>" | python3 -c '...'   # or by hand
# then archive the PDF if open access:
curl -sL -o references/<Key>.pdf "https://europepmc.org/articles/<pmcid lowercase>?pdf=render"
```

The lowercase `pmcid` matters — `?pdf=render` returns HTTP 500 for an uppercase
`PMCID`.

## 8. Second pass

Everything below is deliberately unfinished in the current draft.

### Figures and measurements

- [ ] Capture Fig. 2 from the live app: the three modes, the micro/macro toggle,
      the reference-map modal.
- [ ] `Rscript figures/make_figures.R` for Figs. 3–5; check the placeholders are
      overwritten and rebuild.
- [ ] §4.1 correctness: re-run `GSVA::gsva`/`ssgseaParam` on a random sample of
      gene sets, aggregate as `preparaDados.R` does, and compare against the
      served values. Fill in *n*, Pearson *r*, max |Δ|. Do the same for the
      on-demand path by submitting a precomputed set's member genes.
- [ ] §4.2 face validity: replace the qualitative statements with per-window
      means from the CSV exports for the myelination, GABAergic and neurogenesis
      cases.
- [ ] §4.3 / Table 3: script a warm-latency benchmark of the four endpoints.

### Author and submission metadata

- [ ] Confirm author **order**; add ORCIDs (`\orcid{}` inside `\author`).
- [ ] Replace the two `TODO@ufrn.br` addresses.
- [ ] Full postal affiliation (`\orgdiv`, `\street`, `\postcode`) for the MCB Lab.
- [ ] Funding statement and grant numbers.
- [ ] Author contributions in CRediT terms.

### Repository

- [ ] Choose and add a licence to `MCBLab/BrainMap` — a software paper needs one.
- [ ] Archive the release on Zenodo and cite the DOI in *Code availability*.
- [ ] Merge `VEGAbrain` into `main` and repoint `.github/workflows/pages.yml`
      (its comment already notes this).

### Before submission

- [ ] Re-read the journal's submission guidelines at
      <https://link.springer.com/journal/12021/submission-guidelines> — this
      draft was structured from the Software Original Article definition and the
      SN template, and the live page should be checked for length limits and
      declaration wording, which could not be retrieved automatically (the page
      is behind a bot challenge).
- [ ] Build the reviewers' copy with the `lineno` option.
- [ ] Springer wants a single `.tex`; figures are uploaded separately.
