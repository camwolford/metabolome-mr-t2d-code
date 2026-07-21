# Metabolome-wide Mendelian randomisation identifies metabolites and pathways causally associated with type 2 diabetes

This repository contains R analysis/figure code showing how the manuscript analysis was conducted.

## What this repository contains

The code documents the two-sample Mendelian randomisation workflow from plasma metabolites to type 2 diabetes, fasting glucose and HbA1c. It covers the dual conservative and liberal metabolite quantitative trait locus (mQTL) instrument designs, forward Mendelian randomisation, sensitivity analyses, reverse Mendelian randomisation, colocalisation, and pathway/prioritisation summaries.

## Read the analysis

Read the eight analysis stages in sequence:

1. [Stage 01 — GWAS preparation](analysis/README.md#stage-01-gwas-preparation)
2. [Stage 02 — instrument selection](analysis/README.md#stage-02-instrument-selection)
3. [Stage 03 — forward Mendelian randomisation](analysis/README.md#stage-03-forward-mendelian-randomisation)
4. [Stage 04 — significance filtering](analysis/README.md#stage-04-significance-filtering)
5. [Stage 05 — sensitivity analysis](analysis/README.md#stage-05-sensitivity-analysis)
6. [Stage 06 — Steiger directionality and reverse Mendelian randomisation](analysis/README.md#stage-06-steiger-directionality-and-reverse-mendelian-randomisation)
7. [Stage 07 — colocalisation and downstream integration](analysis/README.md#stage-07-colocalisation-and-downstream-integration)
8. [Stage 08 — pathways and prioritisation](analysis/README.md#stage-08-pathways-and-prioritisation)

## Figure code

The data-figure scripts are linked directly:

1. [Figure 2](figures/R/fig2_volcano.R) — volcano plot.
2. [Figure 3](figures/R/fig3_reverse_heatmap.R) — reverse-Mendelian-randomisation heatmap.
3. [Figure 4](figures/R/fig4_forest.R) — forest plot.
4. [Supplementary Figure 1](figures/R/supp_fig1_prep_data.R) — locus data preparation, with [locus plotting](figures/R/supp_fig1_locuszoom.R).

The [shared figure theme](figures/R/theme_cvdiab.R) and [figure README](figures/README.md) describe common rendering choices, inputs and dependencies.

## Data and scope

This repository contains no GWAS summary statistics, derived tables, rendered figures, binaries, or source assets.

You can read the code without configuring anything. To run a script, its external data, working, output, reference and executable locations need to be supplied separately.

[`config/environment.R`](config/environment.R) can be used to find the location of environment variables (if present) on local computers without including them in the repository. Anyone running a script from this repo with the required data should set the relevant paths first, then load the helper. The following generic example shows the pattern; replace each placeholder with a location on your own computer:

```r
Sys.setenv(
  METABOLOME_MR_INPUT_DIR = "/path/to/authorised-inputs",
  METABOLOME_MR_FULL_METABOLITE_GWAS_DIR = "/path/to/complete-per-metabolite-gwas",
  METABOLOME_MR_WORK_DIR = "/path/to/working-files",
  METABOLOME_MR_OUTPUT_DIR = "/path/to/output-files"
)
source("config/environment.R")
```

`METABOLOME_MR_FULL_METABOLITE_GWAS_DIR` must point to the folder holding the complete per-metabolite summary statistics used by Stages 6 and 7. This is distinct from the sparse Stage 1 input: Stage 1 keeps only validated mQTL rows for instrument selection, whereas the complete files retain the genome-wide rows needed to match outcome instruments and extract regional windows.

Some stages also require an LD reference or an external executable; their file paths are listed in the relevant stage's Inputs table.

## Repository guide

| Repository path | Description |
| --- | --- |
| [`analysis/`](analysis/) | Eight analysis stages in order, with each stage's purpose, scope, inputs and outputs. |
| [`figures/R/`](figures/R/) | Data-figure R scripts. |
| [`config/environment.R`](config/environment.R) | Small R helper that reads the named path settings and checks that they point outside this repository. |
| [`CITATION.cff`](CITATION.cff) | Reserved for citation details after acceptance; currently `TBD`. |
| [`LICENSE`](LICENSE) | MIT licence. |

## External data and tools

Obtain each resource from its provider and comply with its access and use terms.

| Resource | Access route |
| --- | --- |
| Surendran metabolite GWAS (not yet publicly available; academic-use access restrictions apply) | [Omicscience](https://omicscience.org/apps/mgwas/) |
| Suzuki type 2 diabetes GWAS | [DIAGRAM](https://www.diagram-consortium.org/downloads.html) |
| Chen fasting glucose and HbA1c GWAS | [GWAS Catalog: fasting glucose](https://www.ebi.ac.uk/gwas/studies/GCST90002232) and [HbA1c](https://www.ebi.ac.uk/gwas/studies/GCST90002244) |
| LD references: 1000G and UKB | [1000 Genomes](https://mrcieu.github.io/ieugwasr/articles/local_ld.html#:~:text=http%3A//fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz) and [UK Biobank access](https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/) |
| PLINK and PwCoCo | [PLINK](https://www.cog-genomics.org/plink/) and [PwCoCo](https://github.com/jwr-git/pwcoco) |

## Software

The analyses were conducted with R 4.4.1. Package imports appear at the top of the individual R scripts; figure-specific dependencies are listed in the figure README.

## Citation and licence

Citation details are not yet available; [CITATION.cff](CITATION.cff) is currently `TBD`. The original code is available under the [MIT licence](LICENSE).
