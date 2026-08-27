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

## 6. Reference ledger

37 entries in `references.bib`; 23 have the full text archived in
`references/<Key>.pdf`. Bib key ⇄ filename is 1:1 (`Zhu2019.pdf` was renamed to
`Zhu2018.pdf` and `mowinckel2020.pdf` to `Mowinckel2020.pdf` to keep it so).

**Archived locally (23):** `Bahl2017`, `Cao2023`, `Chen2013`, `Chopra2023`,
`GOConsortium2023`, `Gillespie2022`, `Hanzelmann2013`, `Huisman2017`,
`Kanehisa2000`, `Kohler2021`, `Kuleshov2016`, `Lau2008`, `Li2018`,
`Liberzon2011`, `Liberzon2015`, `Markello2021`, `Markello2022`,
`Mowinckel2020`, `Subramanian2005`, `Sunkin2013`, `Wilkinson2016`, `Yi2020`,
`Zhu2018`.

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

## 7. Second pass

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
