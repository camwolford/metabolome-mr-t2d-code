# Analysis code

This guide follows the analysis in order, from GWAS preparation to pathway summaries and prioritisation.

## Stage 01: GWAS preparation

### Purpose

Prepare retained metabolite and glycaemic GWAS inputs for the subsequent analysis stages.

### Scope

Standardise the retained association tables, create one metabolite-association file per metabolite, and count SNP occurrence across the metabolite table.

### Inputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_INPUT_DIR/metabolite_sup_associations.txt` | Retained metabolite association table. |
| `METABOLOME_MR_INPUT_DIR/fasting_glucose_gwas_cleaned.csv` | Retained fasting-glucose association table. |
| `METABOLOME_MR_INPUT_DIR/hbA1c_gwas_cleaned.csv` | Retained HbA1c association table. |
| `METABOLOME_MR_INPUT_DIR/t2dm_gwas_cleaned.tsv` | Type 2 diabetes association table retained for later matching. |

### Outputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_WORK_DIR/metabolite_gwas_associations_cleaned.tsv` | Standardised metabolite association table. |
| `METABOLOME_MR_WORK_DIR/fasting_glucose_gwas_cleaned.tsv` | Standardised fasting-glucose association table. |
| `METABOLOME_MR_WORK_DIR/hbA1c_gwas_cleaned.tsv` | Standardised HbA1c association table. |
| `METABOLOME_MR_WORK_DIR/Individual_Metabolite_GWAS/<metabolite>_GWAS.tsv` | Per-metabolite association tables. |
| `METABOLOME_MR_WORK_DIR/Repeated_SNPS.tsv` | SNP occurrence counts across the metabolite table. |

## Stage 02: Instrument selection

### Purpose

Select conservative and liberal mQTL instruments and prepare them for forward Mendelian randomisation.

### Scope

Create, LD-clump, match and harmonise the conservative and liberal instrument sets with type 2 diabetes, fasting glucose and HbA1c outcomes.

### Inputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_WORK_DIR/Individual_Metabolite_GWAS/<metabolite>_GWAS.tsv` | Per-metabolite association tables from Stage 01. |
| `METABOLOME_MR_WORK_DIR/Repeated_SNPS.tsv` | SNP occurrence counts from Stage 01. |
| `METABOLOME_MR_INPUT_DIR/t2dm_gwas_cleaned.tsv` | Type 2 diabetes outcome association table. |
| `METABOLOME_MR_WORK_DIR/fasting_glucose_gwas_cleaned.tsv` | Fasting-glucose outcome association table. |
| `METABOLOME_MR_WORK_DIR/hbA1c_gwas_cleaned.tsv` | HbA1c outcome association table. |
| `METABOLOME_MR_EUR_PANEL_DIR/EUR.{bed,bim,fam}` | European LD-reference panel files. |
| `METABOLOME_MR_PLINK` | PLINK executable used for LD operations. |

### Outputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_WORK_DIR/02_instrument_selection/conservative/Filtered_IVs/<metabolite>_IVs.tsv` | Conservative eligible instruments. |
| `METABOLOME_MR_WORK_DIR/02_instrument_selection/conservative/Clumped_IVs/<metabolite>_Clumped_IVs.tsv` | LD-clumped conservative instruments. |
| `METABOLOME_MR_WORK_DIR/02_instrument_selection/conservative/Matched_IVs/<metabolite><OUTCOME>_Matched_IVs.tsv` | Conservative instruments matched to an outcome. |
| `METABOLOME_MR_WORK_DIR/02_instrument_selection/liberal/Filtered_IVs/<metabolite>_IVs_Liberal.tsv` | Liberal eligible instruments. |
| `METABOLOME_MR_WORK_DIR/02_instrument_selection/liberal/Clumped_IVs/<metabolite>_Clumped_IVs_Liberal.tsv` | LD-clumped liberal instruments. |
| `METABOLOME_MR_WORK_DIR/02_instrument_selection/liberal/Matched_IVs/<metabolite><OUTCOME>_Matched_IVs_Liberal.tsv` | Liberal instruments matched to an outcome. |
| `METABOLOME_MR_OUTPUT_DIR/02_instrument_selection/conservative/Harmonised_<OUTCOME>_IVs/<metabolite><OUTCOME>_Harmonised_IVs.tsv` | Harmonised conservative instruments for forward MR. |
| `METABOLOME_MR_OUTPUT_DIR/02_instrument_selection/liberal/Harmonised_<OUTCOME>_IVs_Liberal/<metabolite><OUTCOME>_Harmonised_IVs_Liberal.tsv` | Harmonised liberal instruments for forward MR. |

## Stage 03: Forward Mendelian randomisation

### Purpose

Estimate forward metabolite-to-outcome associations for both instrument designs.

### Scope

Calculate forward MR summaries for type 2 diabetes, fasting glucose and HbA1c, then compare instrument counts between the conservative and liberal designs.

### Inputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/02_instrument_selection/conservative/Harmonised_<OUTCOME>_IVs/<metabolite><OUTCOME>_Harmonised_IVs.tsv` | Harmonised conservative instruments from Stage 02. |
| `METABOLOME_MR_OUTPUT_DIR/02_instrument_selection/liberal/Harmonised_<OUTCOME>_IVs_Liberal/<metabolite><OUTCOME>_Harmonised_IVs_Liberal.tsv` | Harmonised liberal instruments from Stage 02. |

### Outputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/03_forward_mr/conservative/<OUTCOME>_MR_Results.tsv` | Conservative forward-MR result for one outcome. |
| `METABOLOME_MR_OUTPUT_DIR/03_forward_mr/conservative/Full_MR_Results.tsv` | Combined conservative forward-MR results. |
| `METABOLOME_MR_OUTPUT_DIR/03_forward_mr/liberal/<OUTCOME>_MR_Results_Liberal.tsv` | Liberal forward-MR result for one outcome. |
| `METABOLOME_MR_OUTPUT_DIR/03_forward_mr/liberal/Full_MR_Results_Liberal.tsv` | Combined liberal forward-MR results. |
| `METABOLOME_MR_OUTPUT_DIR/03_forward_mr/instrument_comparison/metabolites_with_fewer_IVs_<OUTCOME>.tsv` | Metabolites with fewer conservative than liberal instruments for one outcome. |

## Stage 04: Significance filtering

### Purpose

Identify forward associations that pass the manuscript significance-filtering criteria.

### Scope

Filter the conservative and liberal forward-MR summaries independently, then combine qualifying associations with conservative-design precedence.

### Inputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/03_forward_mr/conservative/<OUTCOME>_MR_Results.tsv` | Conservative forward-MR results from Stage 03. |
| `METABOLOME_MR_OUTPUT_DIR/03_forward_mr/liberal/<OUTCOME>_MR_Results_Liberal.tsv` | Liberal forward-MR results from Stage 03. |

### Outputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/04_significance_filtering/conservative/significant_<OUTCOME>_results.tsv` | Filtered conservative result for one outcome. |
| `METABOLOME_MR_OUTPUT_DIR/04_significance_filtering/liberal/significant_<OUTCOME>_results_liberal.tsv` | Filtered liberal result for one outcome. |
| `METABOLOME_MR_OUTPUT_DIR/04_significance_filtering/combined/Full_Significant_Results_Manuscript.tsv` | Combined significant-association summary. |

## Stage 05: Sensitivity analysis

### Purpose

Retain significant metabolite--outcome associations that pass the instrument-strength, heterogeneity and horizontal-pleiotropy filters.

### Scope

Identify conservative-significant membership first, remove those associations from the liberal arm, then apply the F statistic, Cochran's Q and MR-Egger-intercept filters to the conservative and remaining liberal-only rows. A liberal copy cannot rescue a conservative-significant association that fails sensitivity filtering.

### Inputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/04_significance_filtering/combined/Full_Significant_Results_Manuscript.tsv` | Combined significant-association summary from Stage 04. |
| `METABOLOME_MR_OUTPUT_DIR/04_significance_filtering/conservative/significant_<OUTCOME>_results.tsv` | Conservative significant results from Stage 04. |
| `METABOLOME_MR_OUTPUT_DIR/04_significance_filtering/liberal/significant_<OUTCOME>_results_liberal.tsv` | Liberal significant results from Stage 04. |

### Outputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/05_sensitivity_analysis/conservative/Filtered_<OUTCOME>_Results.tsv` | Retained conservative associations for one outcome. |
| `METABOLOME_MR_OUTPUT_DIR/05_sensitivity_analysis/liberal/Filtered_<OUTCOME>_Results_Liberal.tsv` | Retained liberal-only associations for one outcome after conservative-significant membership takes precedence. |
| `METABOLOME_MR_OUTPUT_DIR/05_sensitivity_analysis/combined/Full_Filtered_Results_Manuscript.tsv` | Combined post-sensitivity result summary. |
| `METABOLOME_MR_OUTPUT_DIR/05_sensitivity_analysis/combined/sensitivity_filter_exclusions.tsv` | Sensitivity-filter exclusion record. |

## Stage 06: Steiger directionality and reverse Mendelian randomisation

### Purpose

Evaluate directionality and reverse metabolite--outcome associations after sensitivity analysis.

### Scope

Create Steiger assessments, select outcome instruments, match and harmonise reverse-MR inputs, and calculate reverse-MR summaries.

### Inputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/05_sensitivity_analysis/combined/Full_Filtered_Results_Manuscript.tsv` | Combined post-sensitivity result summary from Stage 05. |
| `METABOLOME_MR_OUTPUT_DIR/05_sensitivity_analysis/{conservative,liberal}/Filtered_<OUTCOME>_Results*.tsv` | Per-design retained associations from Stage 05. |
| `METABOLOME_MR_OUTPUT_DIR/02_instrument_selection/{conservative,liberal}/Harmonised_<OUTCOME>_IVs*/` | Harmonised forward-instrument files from Stage 02. |
| `METABOLOME_MR_INPUT_DIR/t2dm_gwas_cleaned.tsv` | Type 2 diabetes association table. |
| `METABOLOME_MR_WORK_DIR/fasting_glucose_gwas_cleaned.tsv` | Fasting-glucose association table. |
| `METABOLOME_MR_WORK_DIR/hbA1c_gwas_cleaned.tsv` | HbA1c association table. |
| `METABOLOME_MR_FULL_METABOLITE_GWAS_DIR/<metabolite>_GWAS.tsv` | Complete per-metabolite summary statistics for reverse-MR outcome matching, distinct from the sparse Stage 01 tables. |
| `METABOLOME_MR_INPUT_DIR/t2dm_sig_snps_ancestry.csv` | Type 2 diabetes outcome-instrument source. |
| `METABOLOME_MR_INPUT_DIR/fg_hba1c_sig_snps_ancestry.csv` | Fasting-glucose and HbA1c outcome-instrument source. |
| `METABOLOME_MR_INPUT_DIR/metabolite_sample_sizes.tsv` | Metabolite SNP sample-size mapping. |
| `METABOLOME_MR_EUR_PANEL_DIR/EUR.{bed,bim,fam}` | European LD-reference panel files. |
| `METABOLOME_MR_PLINK` | PLINK executable used for outcome-instrument clumping. |

### Outputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/steiger/steiger_results.tsv` | Steiger directionality results. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/steiger/steiger_exclusions.tsv` | Steiger exclusion record. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/steiger/post_steiger_candidates.tsv` | All 54 post-sensitivity candidates with their Steiger-retention status. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/outcome_instruments/<OUTCOME>_{selected,clumped}.tsv` | Selected and clumped reverse-MR outcome instruments. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/matched/<OUTCOME>/` | Matched outcome--metabolite association files. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/matching/{matching_summary,unmatched_instruments}.tsv` | Reverse-MR matching records. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/harmonised/<OUTCOME>/` | Harmonised reverse-MR association files. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/harmonised/harmonisation_summary.tsv` | Reverse-MR harmonisation record. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/results/{t2dm,fg,hba1c}_reverse_mr_raw.tsv` | Outcome-specific raw reverse-MR summaries. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/results/reverse_mr_raw.tsv` | Combined raw reverse-MR summary. |

## Stage 07: Colocalisation and downstream integration

### Purpose

Assess regional colocalisation and produce the final retained-candidate set.

### Scope

Prepare regional association windows, run and classify PwCoCo results, re-run selected forward associations, and integrate downstream evidence.

### Inputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/05_sensitivity_analysis/combined/Full_Filtered_Results_Manuscript.tsv` | Combined post-sensitivity result summary from Stage 05. |
| `METABOLOME_MR_OUTPUT_DIR/05_sensitivity_analysis/{conservative,liberal}/Filtered_<OUTCOME>_Results*.tsv` | Per-design retained associations from Stage 05. |
| `METABOLOME_MR_OUTPUT_DIR/02_instrument_selection/{conservative,liberal}/Harmonised_<OUTCOME>_IVs*/` | Harmonised forward-instrument files from Stage 02. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/steiger/{steiger_results,post_steiger_candidates}.tsv` | Steiger evidence from Stage 06. |
| `METABOLOME_MR_OUTPUT_DIR/06_reverse_mr/results/reverse_mr_raw.tsv` | Reverse-MR evidence from Stage 06. |
| `METABOLOME_MR_INPUT_DIR/t2dm_gwas_cleaned.tsv` | Type 2 diabetes regional association source. |
| `METABOLOME_MR_WORK_DIR/{fasting_glucose_gwas_cleaned,hbA1c_gwas_cleaned}.tsv` | Fasting-glucose and HbA1c regional association sources. |
| `METABOLOME_MR_FULL_METABOLITE_GWAS_DIR/<metabolite>_GWAS.tsv` | Complete per-metabolite summary statistics used as regional association sources, distinct from the sparse Stage 01 tables. |
| `METABOLOME_MR_INPUT_DIR/metabolite_sample_sizes.tsv` | Metabolite SNP sample-size mapping. |
| `METABOLOME_MR_PWCOCO` | PwCoCo executable. |
| `METABOLOME_MR_UKB_EUR_BFILE` | UK Biobank European binary-reference template. |

### Outputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_WORK_DIR/07_colocalisation/candidate_status_manifest.tsv` | Candidate-status record used for downstream integration. |
| `METABOLOME_MR_WORK_DIR/07_colocalisation/regions/locus_eligibility.tsv` | Locus eligibility record. |
| `METABOLOME_MR_WORK_DIR/07_colocalisation/regions/*_{metabolite,outcome}_region.tsv` | Extracted regional association windows. |
| `METABOLOME_MR_WORK_DIR/07_colocalisation/pwcoco/pwcoco_ready_manifest.tsv` | PwCoCo-ready locus manifest. |
| `METABOLOME_MR_WORK_DIR/07_colocalisation/pwcoco/ready/*_{metabolite,outcome}_pwcoco.tsv` | PwCoCo input files. |
| `METABOLOME_MR_WORK_DIR/07_colocalisation/pwcoco/runs/*.coloc` | Raw PwCoCo results. |
| `METABOLOME_MR_WORK_DIR/07_colocalisation/pwcoco/pwcoco_run_manifest.tsv` | PwCoCo run manifest. |
| `METABOLOME_MR_OUTPUT_DIR/07_colocalisation/classification/{raw_pwcoco_rows,locus_classification,association_colocalisation,candidate_colocalisation}.tsv` | Colocalisation classification records. |
| `METABOLOME_MR_OUTPUT_DIR/07_colocalisation/{rerun_mr_results,t2dm_rerun_mr_results,fg_rerun_mr_results,hba1c_rerun_mr_results}.tsv` | Forward-MR re-run results. |
| `METABOLOME_MR_OUTPUT_DIR/07_colocalisation/{candidate_rerun_status,candidate_downstream_status,final_retained_candidates}.tsv` | Final downstream candidate-status records. |

## Stage 08: Pathways and prioritisation

### Purpose

Summarise pathways and outcome-specific instrument sharing for significant metabolites.

### Scope

Create pathway and shared-instrument summaries; confidence-tier annotation is outside this stage.

### Inputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_INPUT_DIR/annotations/responses_combined.tsv` | Metabolite pathway annotations. |
| `METABOLOME_MR_OUTPUT_DIR/04_significance_filtering/combined/Full_Significant_Results_Manuscript.tsv` | Significant-metabolite summary from Stage 04. |
| `METABOLOME_MR_OUTPUT_DIR/02_instrument_selection/conservative/Harmonised_T2DM_IVs/` | Conservative type 2 diabetes harmonised-instrument files. |
| `METABOLOME_MR_OUTPUT_DIR/02_instrument_selection/conservative/Harmonised_FG_IVs/` | Conservative fasting-glucose harmonised-instrument files. |
| `METABOLOME_MR_OUTPUT_DIR/02_instrument_selection/conservative/Harmonised_HBA1C_IVs/` | Conservative HbA1c harmonised-instrument files. |

### Outputs

| File path | What it is |
| --- | --- |
| `METABOLOME_MR_OUTPUT_DIR/08_pathways_and_prioritisation/superpathway_counts.tsv` | Superpathway counts. |
| `METABOLOME_MR_OUTPUT_DIR/08_pathways_and_prioritisation/subpathway_counts.tsv` | Subpathway counts. |
| `METABOLOME_MR_OUTPUT_DIR/08_pathways_and_prioritisation/T2DM_snp_counts.tsv` | Type 2 diabetes shared-SNP counts. |
| `METABOLOME_MR_OUTPUT_DIR/08_pathways_and_prioritisation/FG_snp_counts.tsv` | Fasting-glucose shared-SNP counts. |
| `METABOLOME_MR_OUTPUT_DIR/08_pathways_and_prioritisation/HBA1C_snp_counts.tsv` | HbA1c shared-SNP counts. |
| `METABOLOME_MR_OUTPUT_DIR/08_pathways_and_prioritisation/significant_metabolites_data.tsv` | Significant-metabolite pathway summary. |
