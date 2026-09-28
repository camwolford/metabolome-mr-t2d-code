# Metabolome-wide Mendelian randomisation prioritises candidate causal metabolites and pathways for type 2 diabetes

## Analysis steps

| Folder | Step |
| --- | --- |
| `analysis/01_gwas_preparation` | Standardise the metabolite and outcome GWAS |
| `analysis/02_instrument_selection` | Select, clump, match and harmonise the conservative and liberal instrument sets |
| `analysis/03_forward_mr` | MR of metabolites on type 2 diabetes, fasting glucose and HbA1c |
| `analysis/04_significance_filtering` | Bonferroni and robust-method filters |
| `analysis/05_sensitivity_analysis` | Heterogeneity and pleiotropy filters |
| `analysis/06_reverse_mr` | Steiger filtering, reverse MR and the reverse-direction screen |
| `analysis/07_colocalisation` | PwCoCo colocalisation and MR re-runs for partially colocalising associations |
| `analysis/08_pathways_and_prioritisation` | Pathway summaries and instrument sharing |
| `figures/R` | Figures 2–4 and Supplementary Figure 1 |

## Data and tools

| Resource | Source |
| --- | --- |
| Metabolite GWAS (Surendran et al. 2022; free for academic use) | [Omicscience](https://omicscience.org/apps/mgwas/) |
| Type 2 diabetes GWAS (Suzuki et al. 2024, European ancestry) | [DIAGRAM](https://www.diagram-consortium.org/downloads.html) |
| Fasting glucose and HbA1c GWAS (Chen et al. 2021) | GWAS Catalog [GCST90002232](https://www.ebi.ac.uk/gwas/studies/GCST90002232), [GCST90002244](https://www.ebi.ac.uk/gwas/studies/GCST90002244) |
| LD reference panels | [1000 Genomes EUR](https://mrcieu.github.io/ieugwasr/articles/local_ld.html) |
| External tools | [PLINK](https://www.cog-genomics.org/plink/), [PwCoCo](https://github.com/jwr-git/pwcoco) |

Analyses were run in R 4.4.1 with TwoSampleMR and MendelianRandomization.

## Licence

[MIT](LICENSE)
