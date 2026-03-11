# Multi-Omics-Analysis-of-Young-CSF
Overview
This repository contains a comprehensive bioinformatics pipeline built in R to analyze the effects of Young Cerebrospinal Fluid (CSF) on brain tissue. It processes bulk RNA-seq data to identify transcriptomic signatures, deconvolute cellular proportions, and uncover the epigenetic and transcriptional master regulators driving brain repair (specifically oligodendrocyte maturation and myelination).

Furthermore, the pipeline cross-references the "Young CSF Rescue Signature" against public Alzheimer's Disease datasets to quantify disease reversal and therapeutic potential.

Pipeline Architecture
The analysis is structured into modular phases, allowing for end-to-end execution or standalone analysis once the base DESeq2 object is created.

Phase 1 & 2: Pre-processing & Quality Control

Imports lightweight pseudo-alignments (quant.sf) from Salmon using tximport.

Maps Ensembl IDs to standard Gene Symbols via biomaRt.

Generates global PCA plots for strict quality control.

Phase 3: Tissue Deconvolution

Utilizes BisqueRNA and the Tabula Muris Senis single-cell aging atlas to estimate the physical proportions of Oligodendrocytes, Microglia, and Astrocytes from the bulk tissue.

Phase 5: Differential Gene Expression

Runs DESeq2 to identify significant gene changes (Young CSF vs. Control).

Outputs publication-ready Volcano Plots and hierarchically clustered Heatmaps (pheatmap).

Phase 6: Pathway Enrichment (GSEA)

Maps DESeq2 statistics to ENTREZ IDs and runs Gene Set Enrichment Analysis (GSEA) on the Gene Ontology (GO) Biological Process database via clusterProfiler.

Phase 7: Phenotype Scoring

Calculates a composite z-score based on a literature-derived 11-gene signature to quantify Astrocyte A1 Toxicity/Neuroinflammation.

Phase 8: Alzheimer's Disease Reversal Mapping

Downloads a public 5xFAD vs. WT RNA-seq dataset (GEO: GSE161848).

Calculates an AD disease signature using limma.

Generates a 4-quadrant identity plot to isolate genes suppressed by AD but successfully restored by Young CSF.

Phase 9: Master Regulator Analysis

Queries the Enrichr API (ChEA, ENCODE, TRRUST databases).

Identifies upstream transcription factors (e.g., Klf6, Tcf7) and epigenetic complexes driving the transcriptomic rescue.
