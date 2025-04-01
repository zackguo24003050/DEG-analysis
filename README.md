# DEG-analysis

This project compares gene expression between early-stage (stage I–II) and advanced-stage (stage III–IV) small cell lung cancer (SCLC) tumors.

## Data

- Expression data: `expr_normalized_median_polish.rds`
- Clinical data: `pheno_sclc_ucologne_2015.tsv`

## Method

- Samples are grouped by tumor stage using `uicc_tumor_stage`
- Differential expression is analyzed using the limma package
- Linear modeling and empirical Bayes moderation are applied
- Volcano plot is used to visualize results

## Results

- DEGs (p < 0.01 and |log2FC| > 1): 4

## Output

- DEG results table: `DEG_results_limma.csv`
- Volcano plot (optional)

## Requirements

- R (version ≥ 4.4)
- Packages: `readr`, `limma`, `ggplot2`
