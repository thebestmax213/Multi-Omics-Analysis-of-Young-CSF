suppressPackageStartupMessages({
  library(tximport)
  library(biomaRt)
  
  library(BisqueRNA)
  library(Biobase)
  library(SingleCellExperiment)
  library(TabulaMurisSenisData)
  
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  
  library(GEOquery)
  library(limma)
  library(dplyr)
  
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(gridExtra)
  library(pheatmap)
})

options(timeout = 3600)

## -----------------------------
## USER SETTINGS
## -----------------------------
dir <- "~/Downloads/raw_data"
my_conditions <- c(rep("Treatment", 7), rep("Control", 8))
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")


## -----------------------------
## HELPERS
## -----------------------------
strip_version <- function(x) sub("\\.[0-9]+$", "", x)

get_salmon_files <- function(dir, pattern = "_quant\\.sf$") {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No Salmon quant files found. Check dir/pattern.")
  names(files) <- gsub("_quant\\.sf$", "", basename(files))
  files
}

make_metadata <- function(sample_names, conditions) {
  if (length(sample_names) != length(conditions)) {
    stop("Condition vector length does not match number of samples.")
  }
  data.frame(row.names = sample_names, Condition = factor(conditions))
}

get_tx_maps <- function(mart) {
  bm <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
    mart = mart
  )
  
  tx2gene_symbol <- bm[, c("ensembl_transcript_id", "external_gene_name")]
  colnames(tx2gene_symbol) <- c("TXNAME", "GENEID")
  
  tx2gene_ens <- bm[, c("ensembl_transcript_id", "ensembl_gene_id")]
  colnames(tx2gene_ens) <- c("TXNAME", "GENEID")
  
  id_to_symbol <- unique(bm[, c("ensembl_gene_id", "external_gene_name")])
  colnames(id_to_symbol) <- c("Ensembl_ID", "Symbol")
  id_to_symbol$Ensembl_ID <- strip_version(id_to_symbol$Ensembl_ID)
  
  list(tx2gene_symbol = tx2gene_symbol, tx2gene_ens = tx2gene_ens, id_to_symbol = id_to_symbol)
}

plot_pca_qc <- function(counts, metadata, title = "Bulk RNA-seq PCA (QC Check)", subtitle = NULL) {
  log_counts <- log2(counts + 1)
  gene_vars <- apply(log_counts, 1, var)
  log_counts <- log_counts[gene_vars > 0, , drop = FALSE]
  pca <- prcomp(t(log_counts), scale. = TRUE)
  var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
  
  plot_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Condition = metadata$Condition
  )
  
  ggplot(plot_df, aes(PC1, PC2, color = Condition)) +
    geom_point(size = 4, alpha = 0.85) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)")
    ) +
    theme_bw()
}

plot_cell_type <- function(df, cell_type, color_palette) {
  if (!cell_type %in% colnames(df)) {
    warning(paste(cell_type, "not found in results. Skipping."))
    return(NULL)
  }
  ggplot(df, aes(x = Condition, y = .data[[cell_type]], fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
    stat_compare_means(method = "t.test", label = "p.format",
                       label.x.npc = "center", vjust = -1) +
    scale_fill_manual(values = color_palette) +
    labs(title = cell_type, y = "Proportion") +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

merge_gene_mats_union <- function(mat_a, mat_b) {
  all_genes <- union(rownames(mat_a), rownames(mat_b))
  fill_zeros <- function(mat, all_genes) {
    missing <- setdiff(all_genes, rownames(mat))
    if (length(missing) > 0) {
      zero_mat <- matrix(0, nrow = length(missing), ncol = ncol(mat),
                         dimnames = list(missing, colnames(mat)))
      mat <- rbind(mat, zero_mat)
    }
    mat[all_genes, , drop = FALSE]
  }
  cbind(fill_zeros(mat_a, all_genes), fill_zeros(mat_b, all_genes))
}

# ---- SAFELY CREATE A NAMED STAT VECTOR FOR GSEA ----
make_named_stat <- function(res_obj, dds_obj = NULL) {
  stats <- res_obj$stat
  stats <- stats[!is.na(stats)]
  
  # If names are missing, restore them
  if (is.null(names(stats)) || all(is.na(names(stats))) || all(names(stats) == "")) {
    if (!is.null(rownames(res_obj)) && length(rownames(res_obj)) >= length(res_obj$stat)) {
      # rownames(res_obj) corresponds to full res; subset to non-NA stat
      rn <- rownames(res_obj)[!is.na(res_obj$stat)]
      names(stats) <- rn
    } else if (!is.null(dds_obj)) {
      # fallback: use rownames(dds) aligned to res$stat length (rarely needed)
      if (length(rownames(dds_obj)) == length(res_obj$stat)) {
        rn <- rownames(dds_obj)[!is.na(res_obj$stat)]
        names(stats) <- rn
      }
    }
  }
  
  # Final check
  if (is.null(names(stats)) || all(names(stats) == "")) {
    stop("Could not recover names for stat vector. Check rownames(res) and rownames(dds).")
  }
  
  # strip Ensembl version
  names(stats) <- strip_version(names(stats))
  stats
}

# ---- MAP IDs TO ENTREZ (returns named numeric vector or NULL) ----
map_to_entrez <- function(stat_named, fromType) {
  ids <- unique(names(stat_named))
  mapped <- bitr(ids, fromType = fromType, toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  if (is.null(mapped) || nrow(mapped) == 0) return(NULL)
  
  colnames(mapped)[1] <- "ID"
  df <- data.frame(ID = names(stat_named), stat = as.numeric(stat_named), stringsAsFactors = FALSE)
  merged <- merge(df, mapped, by = "ID")
  if (nrow(merged) == 0) return(NULL)
  
  # Keep strongest stat per ENTREZ
  merged <- merged[order(abs(merged$stat), decreasing = TRUE), ]
  merged <- merged[!duplicated(merged$ENTREZID), ]
  
  out <- merged$stat
  names(out) <- merged$ENTREZID
  sort(out, decreasing = TRUE)
}

run_gsea_bp_safe <- function(res_obj, dds_obj) {
  stat_named <- make_named_stat(res_obj, dds_obj)
  
  # Try ENSEMBL first (most likely)
  gene_list <- map_to_entrez(stat_named, "ENSEMBL")
  
  # If ENSEMBL fails, try SYMBOL
  if (is.null(gene_list)) {
    message("ENSEMBL->ENTREZ mapping returned 0. Trying SYMBOL->ENTREZ...")
    gene_list <- map_to_entrez(stat_named, "SYMBOL")
  }
  
  if (is.null(gene_list) || length(gene_list) < 100) {
    message("GSEA skipped: mapping failed or too few genes after mapping.")
    message("Debug: print head(names(stat_named)) to see the ID format DESeq2 produced.")
    return(NULL)
  }
  
  message(paste("GSEA genes after mapping:", length(gene_list)))
  
  gseGO(
    geneList     = gene_list,
    OrgDb        = org.Mm.eg.db,
    ont          = "BP",
    keyType      = "ENTREZID",
    minGSSize    = 15,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    verbose      = FALSE,
    eps          = 0
  )
}

attach_symbols_from_gse <- function(gset, ad_results) {
  feature_data <- fData(gset)
  ids <- rownames(ad_results)
  
  sym_col <- grep("Symbol", colnames(feature_data), ignore.case = TRUE)[1]
  if (is.na(sym_col)) sym_col <- grep("Gene", colnames(feature_data), ignore.case = TRUE)[1]
  
  if (!is.na(sym_col)) {
    sym <- as.character(feature_data[ids, sym_col])
    sym <- sapply(strsplit(sym, " /// "), `[`, 1)
    ad_results$Symbol <- sym
  } else {
    ad_results$Symbol <- ids
  }
  ad_results
}

find_samples_anycol <- function(df, keyword) {
  apply(df, 1, function(row) any(grepl(keyword, row, ignore.case = TRUE)))
}


## -----------------------------
## PHASE 1: Salmon (SYMBOL) + PCA
## -----------------------------
files <- get_salmon_files(dir)
maps  <- get_tx_maps(mart)

txi_symbol <- tximport(files, type = "salmon", tx2gene = maps$tx2gene_symbol, ignoreTxVersion = TRUE)
raw_counts <- txi_symbol$counts
write.csv(raw_counts, "Phase1_Raw_Count_Matrix.csv")

preview_table <- data.frame(FileName = names(files), Assigned_Condition = my_conditions)
print(preview_table)

metadata <- make_metadata(colnames(raw_counts), my_conditions)
print(plot_pca_qc(raw_counts, metadata, subtitle = "Prior Manipulation"))


## -----------------------------
## PHASE 3: Bisque Deconvolution
## -----------------------------
sce_myeloid <- TabulaMurisSenisFACS(tissues = "Brain_Myeloid")
sce_nonmyel <- TabulaMurisSenisFACS(tissues = "Brain_Non-Myeloid")

sce_m_raw  <- if (is.list(sce_myeloid)) sce_myeloid[[1]] else sce_myeloid
sce_nm_raw <- if (is.list(sce_nonmyel)) sce_nonmyel[[1]] else sce_nonmyel

sce_m_hip  <- sce_m_raw[,  grepl("Hippocampus", sce_m_raw$subtissue,  ignore.case = TRUE)]
sce_nm_hip <- sce_nm_raw[, grepl("Hippocampus", sce_nm_raw$subtissue, ignore.case = TRUE)]

print(paste("Myeloid Cells:", ncol(sce_m_hip)))
print(paste("Non-Myeloid Cells:", ncol(sce_nm_hip)))

counts_m  <- as.matrix(counts(sce_m_hip))
counts_nm <- as.matrix(counts(sce_nm_hip))
final_counts <- merge_gene_mats_union(counts_m, counts_nm)

cols_keep <- c("subtissue", "cell_ontology_class", "mouse.id")
meta_m  <- colData(sce_m_hip)[, cols_keep, drop = FALSE]
meta_nm <- colData(sce_nm_hip)[, cols_keep, drop = FALSE]
final_meta <- rbind(meta_m, meta_nm)

valid_cells <- !is.na(final_meta$cell_ontology_class)
final_counts <- final_counts[, valid_cells, drop = FALSE]
final_meta   <- final_meta[valid_cells, , drop = FALSE]

ref_eset <- ExpressionSet(
  assayData = final_counts,
  phenoData = AnnotatedDataFrame(as.data.frame(final_meta))
)

overlap <- intersect(rownames(raw_counts), rownames(ref_eset))
print(paste("Overlapping genes found:", length(overlap)))
if (length(overlap) < 2000) stop("Low overlap: check ID types (symbol vs ensembl).")

bulk_eset <- ExpressionSet(assayData = raw_counts[overlap, , drop = FALSE])

res_bisque <- ReferenceBasedDecomposition(
  bulk.eset = bulk_eset,
  sc.eset = ref_eset,
  cell.types = "cell_ontology_class",
  subject.names = "mouse.id",
  use.overlap = FALSE
)

estimated_props <- t(res_bisque$bulk.props)
colnames(estimated_props) <- gsub("microglial cell", "Microglia", colnames(estimated_props))
colnames(estimated_props) <- gsub("oligodendrocyte", "Oligodendrocytes", colnames(estimated_props))
colnames(estimated_props) <- gsub("astrocyte", "Astrocytes", colnames(estimated_props))
colnames(estimated_props) <- gsub("endothelial cell", "Endothelial", colnames(estimated_props))

plot_props <- as.data.frame(estimated_props)
plot_props$Condition <- metadata[rownames(plot_props), "Condition"]
write.csv(plot_props, "Phase3_Final_Cell_Proportions.csv")
print("Deconvolution complete. Saved to Phase3_Final_Cell_Proportions.csv")

my_colors <- c("Control" = "#F8766D", "Treatment" = "#00BFC4")
grid.arrange(
  plot_cell_type(plot_props, "Oligodendrocytes", my_colors),
  plot_cell_type(plot_props, "Microglia", my_colors),
  plot_cell_type(plot_props, "Astrocytes", my_colors),
  ncol = 3
)


## -----------------------------
## PHASE 5: DESeq2 (ENSEMBL) + annotate
## -----------------------------
txi_ens <- tximport(files, type = "salmon", tx2gene = maps$tx2gene_ens, ignoreTxVersion = TRUE)
metadata_ens <- make_metadata(colnames(txi_ens$counts), my_conditions)

dds <- DESeqDataSetFromMatrix(
  countData = round(txi_ens$counts),
  colData = metadata_ens,
  design = ~ Condition
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds$Condition <- relevel(dds$Condition, ref = "Control")

dds <- DESeq(dds)
res <- results(dds, contrast = c("Condition", "Treatment", "Control"))

# Save DEGs with symbols
res_df <- as.data.frame(res)
res_df$Ensembl_ID <- strip_version(rownames(res_df))
res_final <- merge(res_df, maps$id_to_symbol, by = "Ensembl_ID", all.x = TRUE)
res_final <- res_final[order(res_final$padj), ]
write.csv(res_final, "Phase5_DEGs_Final.csv", row.names = FALSE)


## -----------------------------
## Volcano + Heatmap
## -----------------------------
hero_genes <- c("Mbp","Plp1","Mag","Mog","Olig1","Olig2","Sox10","Cnp","Myrf")

res_final$Significant <- ifelse(!is.na(res_final$padj) &
                                  res_final$padj < 0.05 &
                                  abs(res_final$log2FoldChange) > 0.5, "Sig", "NS")
res_final$Label <- ifelse(res_final$Symbol %in% hero_genes & res_final$Significant == "Sig",
                          res_final$Symbol, NA)

p_volcano <- ggplot(res_final, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = Significant), alpha = 0.6) +
  scale_color_manual(values = c("NS" = "grey", "Sig" = "firebrick")) +
  geom_text_repel(aes(label = Label), box.padding = 0.5, max.overlaps = Inf, fontface = "bold") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot: Young CSF Effect",
       subtitle = "Highlighted genes are Oligodendrocyte markers") +
  theme_bw()
print(p_volcano)

vsd <- vst(dds, blind = FALSE)
top_genes <- head(res_final$Ensembl_ID[order(res_final$padj)], 30)

mat <- assay(vsd)
rownames(mat) <- strip_version(rownames(mat))
mat_top <- mat[top_genes, , drop = FALSE]

sym_lookup <- maps$id_to_symbol$Symbol[match(rownames(mat_top), maps$id_to_symbol$Ensembl_ID)]
rownames(mat_top) <- ifelse(is.na(sym_lookup), rownames(mat_top), sym_lookup)

mat_scaled <- t(scale(t(mat_top)))
ann_col <- data.frame(Condition = colData(dds)$Condition)
rownames(ann_col) <- colnames(dds)

pheatmap(
  mat_scaled,
  cluster_rows = TRUE, cluster_cols = TRUE,
  annotation_col = ann_col,
  main = "Top DEGs (Young CSF vs Control)",
  color = colorRampPalette(c("navy","white","firebrick"))(50)
)


## -----------------------------
## PHASE 6: GSEA (FIXED FOR NAMELESS STAT VECTOR)
## -----------------------------
gse <- run_gsea_bp_safe(res, dds)

if (!is.null(gse) && nrow(gse@result) > 0) {
  p_gsea <- dotplot(gse, showCategory = 10, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle("GSEA: Young CSF Pathways")
  print(p_gsea)
  print(head(gse@result[, c("Description", "NES", "p.adjust")]))
} else {
  print("No significant gseGO pathways at adjusted p < 0.05 (or mapping failed).")
}


## -----------------------------
## PHASE 7: Astrocyte A1 Toxicity Score
## -----------------------------
a1_genes_symbol <- c("C3","H2-D1","Serpina3n","Ggta1","Gbp2","Fbln5",
                     "Fkbp5","Cryab","Srgn","Amigo2","Lcn2")

mat_vsd <- assay(vsd)
rownames(mat_vsd) <- strip_version(rownames(mat_vsd))

sym <- maps$id_to_symbol$Symbol[match(rownames(mat_vsd), maps$id_to_symbol$Ensembl_ID)]
mat_sym <- mat_vsd
rownames(mat_sym) <- sym
mat_sym <- mat_sym[!is.na(rownames(mat_sym)), , drop = FALSE]

valid_a1 <- intersect(a1_genes_symbol, rownames(mat_sym))
print(paste("Found", length(valid_a1), "A1 marker genes."))

if (length(valid_a1) > 2) {
  a1_mat <- mat_sym[valid_a1, , drop = FALSE]
  a1_scaled <- t(scale(t(a1_mat)))
  toxicity_scores <- colMeans(a1_scaled, na.rm = TRUE)
  
  tox_df <- data.frame(
    Sample = names(toxicity_scores),
    Condition = colData(dds)$Condition,
    A1_Score = toxicity_scores
  )
  
  pval <- t.test(A1_Score ~ Condition, data = tox_df)$p.value
  
  p_tox <- ggplot(tox_df, aes(Condition, A1_Score, fill = Condition)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3) +
    scale_fill_manual(values = my_colors) +
    labs(
      title = "Astrocyte A1 Toxicity Score",
      subtitle = paste0("t-test p = ", signif(pval, 3)),
      y = "Composite Z-score"
    ) +
    theme_bw()
  print(p_tox)
} else {
  print("Not enough A1 genes found to calculate score.")
}


## -----------------------------
## PHASE 8: AD Reversal (GSE161848)
## -----------------------------
my_res <- read.csv("Phase5_DEGs_Final.csv")
my_res_clean <- my_res %>%
  filter(!is.na(padj), !is.na(Symbol)) %>%
  transmute(Symbol, My_LogFC = log2FoldChange, My_Padj = padj)

print("Downloading GSE161848...")
gset <- getGEO("GSE161848", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
pdata <- pData(gset)

print("Identifying WT Vehicle vs 5xFAD Vehicle...")
is_vehicle <- find_samples_anycol(pdata, "vehicle")
is_wt <- find_samples_anycol(pdata, "N-tg") | find_samples_anycol(pdata, "Non-transgenic")
is_ad <- find_samples_anycol(pdata, "5xFAD")

group_wt <- which(is_wt & is_vehicle)
group_ad <- which(is_ad & is_vehicle)
print(paste("WT Vehicle:", length(group_wt), "| 5xFAD Vehicle:", length(group_ad)))
if (length(group_wt) == 0 || length(group_ad) == 0) stop("Could not identify groups from metadata.")

ex <- exprs(gset)
if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)

idx <- c(group_ad, group_wt)
ex_subset <- ex[, idx, drop = FALSE]

groups <- factor(c(rep("AD", length(group_ad)), rep("WT", length(group_wt))))
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(ex_subset, design)
fit2 <- contrasts.fit(fit, makeContrasts(AD - WT, levels = design))
fit2 <- eBayes(fit2)

ad_results <- topTable(fit2, number = Inf, adjust.method = "BH")
ad_results <- attach_symbols_from_gse(gset, ad_results)

ad_sig <- ad_results %>%
  transmute(Symbol, AD_LogFC = logFC, AD_Padj = adj.P.Val)

merged_df <- inner_join(my_res_clean, ad_sig, by = "Symbol")

sig_genes <- merged_df %>%
  filter(My_Padj < 0.05, abs(My_LogFC) > 0.25) %>%
  mutate(Category = case_when(
    My_LogFC > 0 & AD_LogFC < 0 ~ "Rescue (Restored)",
    My_LogFC < 0 & AD_LogFC > 0 ~ "Rescue (Suppressed)",
    My_LogFC > 0 & AD_LogFC > 0 ~ "Concordant Up",
    My_LogFC < 0 & AD_LogFC < 0 ~ "Concordant Down",
    TRUE ~ "Other"
  ))

p_identity <- ggplot(sig_genes, aes(AD_LogFC, My_LogFC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = Category), alpha = 0.85, size = 3) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  geom_text_repel(aes(label = Symbol), box.padding = 0.5, point.padding = 0.3,
                  max.overlaps = 15, size = 3.5, fontface = "bold") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
  scale_color_manual(values = c(
    "Rescue (Restored)" = "dodgerblue",
    "Rescue (Suppressed)" = "forestgreen",
    "Concordant Up" = "firebrick",
    "Concordant Down" = "orange",
    "Other" = "grey"
  )) +
  labs(
    title = "Identity of Rescue Genes (Young CSF vs 5xFAD)",
    subtitle = "Quadrants show reversal vs concordance",
    x = "AD effect (5xFAD logFC)",
    y = "Young CSF effect (Treatment log2FC)"
  ) +
  theme_bw(base_size = 14)
print(p_identity)

ad_down <- ad_sig$Symbol[ad_sig$AD_LogFC < -0.25 & ad_sig$AD_Padj < 0.05]
csf_up  <- my_res_clean$Symbol[my_res_clean$My_LogFC > 0.25 & my_res_clean$My_Padj < 0.05]
ov <- intersect(ad_down, csf_up)
print(paste("Restored genes (AD down, CSF up):", length(ov)))

universe <- nrow(merged_df)
mat_f <- matrix(
  c(length(ov),
    length(csf_up) - length(ov),
    length(ad_down) - length(ov),
    universe - length(ad_down) - length(csf_up)),
  nrow = 2
)
print(paste("Fisher enrichment p:", format.pval(fisher.test(mat_f, alternative = "greater")$p.value)))


## -----------------------------
## PHASE 9: Upstream TF Enrichment (Enrichr)
## -----------------------------
if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
library(enrichR)

setEnrichrSite("Enrichr")
dbs <- c("ChEA_2016", "TRRUST_Transcription_Factors_2019", "ENCODE_TF_ChIP-seq_2015")

up_genes <- my_res_clean$Symbol[my_res_clean$My_LogFC > 0.25 & my_res_clean$My_Padj < 0.05]
print(paste("Upregulated rescue genes for TF enrichment:", length(up_genes)))

if (length(up_genes) > 0) {
  enriched <- enrichr(up_genes, dbs)
  tf_tbl <- enriched[["TRRUST_Transcription_Factors_2019"]]
  tf_sig <- tf_tbl %>% filter(Adjusted.P.value < 0.05)
  
  if (nrow(tf_sig) > 0) {
    top_tfs <- head(tf_sig, 10)
    p_tf <- ggplot(top_tfs, aes(reorder(Term, -log10(Adjusted.P.value)), -log10(Adjusted.P.value))) +
      geom_bar(stat = "identity", alpha = 0.75) +
      coord_flip() +
      labs(title = "Top Upstream Regulators (TRRUST)",
           subtitle = "Predicted TFs for CSF-up rescue genes",
           x = "Transcription factor", y = "-log10 adj p") +
      theme_bw()
    print(p_tf)
    print(head(tf_sig[, c("Term", "Adjusted.P.value", "Genes")]))
  } else {
    message("No TRRUST TFs at adj p < 0.05; showing ChEA top hits.")
    chea_tbl <- enriched[["ChEA_2016"]]
    top_chea <- head(chea_tbl[order(chea_tbl$Adjusted.P.value), ], 10)
    p_chea <- ggplot(top_chea, aes(reorder(Term, -log10(Adjusted.P.value)), -log10(Adjusted.P.value))) +
      geom_bar(stat = "identity", alpha = 0.75) +
      coord_flip() +
      labs(title = "Top Upstream Regulators (ChEA)",
           subtitle = "Fallback if TRRUST has no hits",
           x = "Transcription factor", y = "-log10 adj p") +
      theme_bw()
    print(p_chea)
  }
} else {
  print("Not enough upregulated genes for Enrichr TF enrichment.")
}
