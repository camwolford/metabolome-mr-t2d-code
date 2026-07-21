# Figure code

| Figure | Code | What it shows |
| --- | --- | --- |
| Figure 2 | [Volcano plot](R/fig2_volcano.R) | Forward-Mendelian-randomisation significance across metabolites and outcomes. |
| Figure 3 | [Reverse-MR heatmap](R/fig3_reverse_heatmap.R) | Reverse Mendelian randomisation estimates across metabolites and outcomes. |
| Figure 4 | [Forest plot](R/fig4_forest.R) | Prioritised metabolite effect estimates for type 2 diabetes, fasting glucose and HbA1c. |
| Supplementary Figure 1 | [Locus-data preparation](R/supp_fig1_prep_data.R) and [locus plotting](R/supp_fig1_locuszoom.R) | Regional association and colocalisation evidence at selected loci. |

## External inputs and outputs

The scripts read authorised external inputs from configured locations and write vector PDFs only to `METABOLOME_MR_OUTPUT_DIR/figures/`. No input datasets, derived result tables or rendered figures are included here.

Figures 2 and 4 use an external, author-supplied Table 2 annotation seed. They read and check its confidence labels but do not calculate confidence tiers. Supplementary Figure 1 requires authorised regional association and colocalisation inputs. Figure 1 and Supplementary Figures 2–5 are non-R figures and are not included.

## Dependencies

Use R 4.4.1. The packages and external tools required by each script are:

| Code | Packages and tools |
| --- | --- |
| Figure 2 volcano plot | `readr`, `dplyr`, `stringr`, `ggplot2`, `ggrepel`, `patchwork`, `cowplot` |
| Figure 3 reverse-MR heatmap | `readr`, `dplyr`, `ggplot2`, `scales` |
| Figure 4 forest plot | `readr`, `dplyr`, `stringr`, `purrr`, `ggplot2` |
| Supplementary Figure 1 preparation | `readr`, `dplyr`; external `PLINK` executable |
| Supplementary Figure 1 locus plot | `readr`, `dplyr`, `stringr`, `locuszoomr`, `cowplot`, `ggplot2`, `EnsDb.Hsapiens.v75` |
| Shared figure theme | `ggplot2` |
