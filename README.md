## Promoter Activity Analysis for CPAbP and nonCPAbP Genes

This repository contains a collection of R scripts and R Markdown notebooks used to analyze **promoter activities of CPAbP (Cleavage and Polyadenylation site–Associated Promoter Binding) and nonCPAbP genes** for the accompanying publication.

The analyses include gene classification, read quantification near promoters, differential expression testing, and metagene visualization. All code is designed to be reproducible and modular.

---

## Repository Structure

- **`Define_CPAbP_genes.Rmd`**
  - Identifies CPAbP and nonCPAbP genes from annotated gene models.
  - Outputs lists of gene IDs and associated promoter coordinates.

- **`HITstat_analysis.R`**
  - Analyzes high-throughput interaction data (e.g., HIT statistics).
  - Generates gene-level summary statistics and visualizations for promoter-linked interactions.

- **`differential_reads_analysis.Rmd`**
  - Counts sequencing reads near transcription start sites (eg. TSS ±200 bp) for CPAbP and nonCPAbP genes.
  - Performs differential analysis between control and U1-depleted conditions using `DESeq2`.
  - Produces statistical comparisons and boxplots of log2 fold changes.

- **`metagene_plot.R`**
  - Generates metagene plots of read enrichment across TSS regions of different gene groups.
  - Useful for visualizing global promoter behavior.

---

## Author

GyeungYun Kim

---

## Cite

If you use this script in your research, please cite: Kim et al., Mol Cell, 2025.

