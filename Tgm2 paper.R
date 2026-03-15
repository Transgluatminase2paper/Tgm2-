#Figure 2B PCA plot
library(ggplot2)
library(factoextra)
data <- read.csv('Proteomics_data.csv', row.names = 1)
data_t <- t(data)
pca_res <- prcomp(data_t, scale = TRUE)
explained_var <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100
pca_df <- as.data.frame(pca_res$x)  # PC scores for each sample
pca_df$Group <- factor(gsub("_\\d+", "", rownames(pca_df)))

pdf("verify_PCA_Protein_Expression.pdf", width = 8, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = rownames(pca_df))) +
  geom_point(size = 4, alpha = 0.8) + 
  stat_ellipse(aes(group = Group), type = "norm",  linetype = 2) +  
  labs(title = "PCA Plot of Protein Expression",
       x = paste0("PC1 (", round(explained_var[1], 2), "% Variance)"),
       y = paste0("PC2 (", round(explained_var[2], 2), "% Variance)")) +
  theme_minimal() +
  theme(legend.title = element_blank()) 
dev.off()


#Figure 2C Heatmap using pheatmap for all the protein expression data
library(grid)
library(pheatmap)
library(dplyr)
data <- read.csv('Proteomics_data.csv', row.names = 1)
data <- na.omit(data)
rownames(data) <- NULL
annotation_col <- data.frame(
  Group = c(rep("LFD", 8), rep("HFD", 8), rep("HFDinf", 8)) 
)
rownames(annotation_col) <- colnames(data)
pdf("Protein expression heatmap.pdf", width = 8, height = 8)
pheatmap(as.matrix(data), 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         scale = "row", # Normalize rows to show Z-scores
         show_rownames = FALSE)
dev.off()



#Figure 2D, Pie split up
Before the pie distribution of the (non) significantly altered proteins, differential expression analysis was done using Limma.
library(limma)
data <- read.csv('Proteomics_data.csv', sep = ";", check.names = FALSE)
colnames(data)
data <- data[ , 1:25]
data$Name <- make.unique(as.character(data$Name))
rownames(data) <- data$Name
data$Name <- NULL
data[] <- lapply(data, function(x) as.numeric(as.character(x)))
mat <- as.matrix(data)
expr_data <- mat
group <- factor(c(rep("LFD", 8), rep("HFD", 8), rep("HFD_inf", 8)))
design <- model.matrix(~ 0 + group)  
colnames(design) <- levels(group)
# Fit linear model
fit <- lmFit(expr_data, design)
# For contrast definition 
contrast_matrix <- makeContrasts(
  HFD_vs_LFD = HFD - LFD,
  HFD_inf_vs_HFD = HFD_inf - HFD,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
# Extract differential proteins
# Disease effect: HFD vs LFD
results_HFD_vs_LFD <- topTable(fit2, coef = "HFD_vs_LFD", adjust = "fdr", number = Inf)
# Treatment effect: HFD_inf vs HFD
results_HFD_inf_vs_HFD <- topTable(fit2, coef = "HFD_inf_vs_HFD", adjust = "fdr", number = Inf)
write.csv(results_HFD_vs_LFD, "limma_HFD_vs_LFD_results.csv")
write.csv(results_HFD_inf_vs_HFD, "limma_HFD_inf_vs_HFD_results.csv")
#Significance threshold definition
adj_pval_threshold <- 0.05
logFC_threshold <- 0.5
significant_HFD_vs_LFD <- subset(results_HFD_vs_LFD,
                                 adj.P.Val < adj_pval_threshold & abs(logFC) > logFC_threshold)
significant_HFD_inf_vs_HFD <- subset(results_HFD_inf_vs_HFD,
                                     adj.P.Val < adj_pval_threshold & abs(logFC) > logFC_threshold)
write.csv(significant_HFD_vs_LFD, "significant_HFD_vs_LFD.csv", row.names = TRUE)
write.csv(significant_HFD_inf_vs_HFD, "significant_HFD_inf_vs_HFD.csv", row.names = TRUE)

#Pie split up function
library(ggplot2)
library(cowplot)
df <- data.frame(
  group = c("Altered", "Not significantly altered"),
  n     = c(250, 3591)
)
df$fraction <- df$n / sum(df$n)
df$ymax <- cumsum(df$fraction)
df$ymin <- c(0, head(df$ymax, -1))
df$y    <- (df$ymin + df$ymax) / 2
df$label <- c("Altered\n(n=250)", "Not significantly altered\n(n=3591)")
cols <- c("Altered" = "#FBB042", "Not significantly altered" = "#C3996C")
p <- ggplot(df) +
  geom_rect(
    aes(ymin = ymin, ymax = ymax, xmin = 1, xmax = 2, fill = group),
    color = NA, linewidth = 2
  ) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.2) +                  
  scale_fill_manual(values = cols) +
  geom_text(aes(x = 1.45, y = y, label = label), size = 3.6) +
  theme_void() +
  theme(legend.position = "none")
ggsave("donut_plotproteins.png",
       plot = p,
       width = 6, height = 6, dpi = 600)
ggsave("donut_plot_proteins.pdf",
       plot = p,
       width = 6, height = 6,
       device = cairo_pdf)  

Figure 2E, combined volcano plot
The results from the differential expression analysis is thus visualized using the volcano plot including both contrasts
results_HFD_vs_LFD <- topTable(fit2, coef = "HFD_vs_LFD", adjust = "fdr", number = Inf)
results_HFD_inf_vs_HFD <- topTable(fit2, coef = "HFD_inf_vs_HFD", adjust = "fdr", number = Inf)
results_HFD_vs_LFD$Feature <- rownames(results_HFD_vs_LFD)
results_HFD_inf_vs_HFD$Feature <- rownames(results_HFD_inf_vs_HFD)
adj_pval_threshold <- 0.05
logFC_threshold <- 0.5
results_HFD_vs_LFD$Comparison <- "Disease (HFD vs LFD)"
results_HFD_inf_vs_HFD$Comparison <- "Treatment (HFDinf vs HFD)"

volcano_data <- rbind(results_HFD_vs_LFD, results_HFD_inf_vs_HFD)
volcano_data$negLog10AdjP <- -log10(volcano_data$adj.P.Val)
volcano_data$Significant <- ifelse(
  volcano_data$adj.P.Val < adj_pval_threshold & abs(volcano_data$logFC) > logFC_threshold,
  "TRUE", "FALSE"
)
volcano_data$Comparison <- factor(
  volcano_data$Comparison,
  levels = c("Disease (HFD vs LFD)", "Treatment (HFDinf vs HFD)")
)
label_features <- c("Fabp5", "Gstp1", "Ighg1", "Gm5629", "Vps4a",
                    "Tgm2", "Cyp2b9", "Pex11a", "Acaa1b", "S100a9", "Ear2")
label_data <- subset(volcano_data, Feature %in% label_features)
p <- ggplot(volcano_data, aes(x = logFC, y = negLog10AdjP)) +
  geom_point(aes(color = Comparison, alpha = Significant), size = 2) +
  scale_color_manual(values = c(
    "Disease (HFD vs LFD)" = "blue",
    "Treatment (HFDinf vs HFD)" = "red"
  )) +
  scale_alpha_manual(values = c("FALSE" = 0.25, "TRUE" = 0.95)) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(adj_pval_threshold),
             linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = label_data,
    aes(label = Feature),
    size = 5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "black",
    show.legend = FALSE
  ) +
  labs(
    title = "Volcano Plot: Disease and Treatment Effects",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Comparison",
    alpha = "Significant"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )
ggsave("combined_volcano_plot.png", p, width = 15, height = 10, dpi = 300)

#Figure 3F. Enrichment analysis including GO and KEGG
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(readr)
keys(org.Mm.eg.db, keytype = "SYMBOL")[1:10]
dep_data <- read.csv("DEPs_conditions_combined.csv", header = FALSE)
dep_genes <- dep_data[[1]]
dep_genes <- dep_genes[dep_genes != "Protein name"]

background_data <- read.csv("Proteomics_data.csv", header = FALSE)
background_genes <- background_data[[1]]
background_genes <- background_genes[background_genes != "Protein name"]

dep_entrez <- bitr(dep_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
background_entrez <- bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

ego <- enrichGO(
  gene = dep_entrez$ENTREZID,
  universe = background_entrez$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",  # Options: "BP", "MF", "CC", or "ALL"
  pAdjustMethod = "fdr",
  pvalueCutoff = 1.0,
  qvalueCutoff = 1.0,
  readable = TRUE  
)
go_results <- as.data.frame(ego)
write.csv(go_results, "GO_Enrichment_Results_Full2.csv", row.names = FALSE)

#KEGG analysis
kegg_results <- enrichKEGG(
  gene = dep_entrez$ENTREZID,
  universe = background_entrez$ENTREZID,
  organism = "mmu",           # mmu = Mus musculus
  pvalueCutoff = 1.0,
  qvalueCutoff = 1.0,
  pAdjustMethod = "fdr"
)
all_entrez_ids <- unique(unlist(strsplit(kegg_results@result$geneID, "/")))
mapping <- bitr(all_entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)

kegg_results@result$geneSymbols <- sapply(kegg_results@result$geneID, function(x) {
  ids <- unlist(strsplit(x, "/"))
  syms <- mapping$SYMBOL[match(ids, mapping$ENTREZID)]
  paste(na.omit(syms), collapse = "/")
})
kegg_df <- as.data.frame(kegg_results)
write.csv(kegg_df, "KEGG_Enrichment_Results_Full.csv", row.names = FALSE)


#Figure 3B
library(mixOmics)
library(ropls)
library(ggplot2)
data <- read.csv('Metabolites_clean.csv', sep = ";", check.names = FALSE)
colnames(data)
data <- data[ , 1:25]
data$Name <- make.unique(as.character(data$Name))
rownames(data) <- data$Name
data$Name <- NULL
data[] <- lapply(data, function(x) as.numeric(as.character(x)))
mat <- as.matrix(data)
X <- t(mat)  # Samples as rows
Y <- factor(c(rep("LFD", 8), rep("HFD", 8), rep("HFD_inf", 8)))  # Group labels
#  Run PLS-DA
plsda_result_all <- plsda(X, Y, ncomp = 2)
scores_df <- as.data.frame(plsda_result_all$variates$X)
scores_df$Group <- Y
png("PLSDA_allmetabolites_metabolomics.png", width = 800, height = 600)
plotIndiv(plsda_result_all, comp = c(1,2),
          group = Y,
          legend = TRUE,
          title = "PLS-DA: All Metabolites",
          ellipse = TRUE)
pdf("PLSDA_allmetabolites_metabolomics.pdf", width = 8, height = 6)
plotIndiv(plsda_result_all, comp = c(1,2),
          group = Y,
          legend = TRUE,
          title = "PLS-DA: All Metabolites",
          ellipse = TRUE)
dev.off()
var_exp_x1 <- round(plsda_result_all$explained_variance$X[1] * 100, 1)
var_exp_x2 <- round(plsda_result_all$explained_variance$X[2] * 100, 1)
p <- ggplot(scores_df, aes(x = comp1, y = comp2, color = Group)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(type = "norm", size = 1, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "PLS-DA: All Metabolites",
       x = paste0("X-variate 1 (", var_exp_x1, "31.2%)"),
       y = paste0("X-variate 2 (", var_exp_x2, "14.5%)")) +
  scale_color_manual(values = c("blue", "orange", "red")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("PLSDA _Xvariate_metabolomics.png", plot = p, width = 8, height = 6, dpi = 300)
ggsave("PLSDA _Xvariate_metabolomics.pdf", plot = p, width = 8, height = 6)

















#Figure 3C: Metabolites, heatmap data
library(grid)
library(pheatmap)
library(dplyr)
data <- read.csv('Metabolites_clean.csv', sep = ";", check.names = FALSE)
colnames(data)
data <- data[ , 1:25]
data$Name <- make.unique(as.character(data$Name))
rownames(data) <- data$Name
data$Name <- NULL
data[] <- lapply(data, function(x) as.numeric(as.character(x)))
mat <- as.matrix(data)
#Define group vectors before the heatmap
groups <- c(rep("LFD", 8),     
            rep("HFD", 8),     
            rep("HFD_inf", 8)) 
annotation_col <- data.frame(
  Group = factor(groups)
)
rownames(annotation_col) <- colnames(mat)
annotation_colors <- list(
  Group = c(LFD = "red", 
            HFD = "green", 
            HFD_inf = "orange")
)
pheatmap(mat_clean,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         filename = "heatmap_metabolites.png")


#Figure 3D, Pie split up
#Before the pie distribution of the (non) significantly altered metabolites, differential expression analysis was done using Limma.
library(limma)
data <- read.csv('Metabolites.csv', sep = ";", check.names = FALSE)
colnames(data)
data <- data[ , 1:25]
data$Name <- make.unique(as.character(data$Name))
rownames(data) <- data$Name
data$Name <- NULL
data[] <- lapply(data, function(x) as.numeric(as.character(x)))
mat <- as.matrix(data)
expr_data <- mat
group <- factor(c(rep("LFD", 8), rep("HFD", 8), rep("HFD_inf", 8)))
design <- model.matrix(~ 0 + group)  
colnames(design) <- levels(group)
# Fit linear model
fit <- lmFit(expr_data, design)
# For contrast definition 
contrast_matrix <- makeContrasts(
  HFD_vs_LFD = HFD - LFD,
  HFD_inf_vs_HFD = HFD_inf - HFD,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
# Extract differential metabolites
# Disease effect: HFD vs LFD
results_HFD_vs_LFD <- topTable(fit2, coef = "HFD_vs_LFD", adjust = "fdr", number = Inf)
# Treatment effect: HFD_inf vs HFD
results_HFD_inf_vs_HFD <- topTable(fit2, coef = "HFD_inf_vs_HFD", adjust = "fdr", number = Inf)
write.csv(results_HFD_vs_LFD, "limma_HFD_vs_LFD_results.csv")
write.csv(results_HFD_inf_vs_HFD, "limma_HFD_inf_vs_HFD_results.csv")
#Significance threshold definition
adj_pval_threshold <- 0.05
logFC_threshold <- 0.5
significant_HFD_vs_LFD <- subset(results_HFD_vs_LFD,
                                 adj.P.Val < adj_pval_threshold & abs(logFC) > logFC_threshold)
significant_HFD_inf_vs_HFD <- subset(results_HFD_inf_vs_HFD,
                                     adj.P.Val < adj_pval_threshold & abs(logFC) > logFC_threshold)
write.csv(significant_HFD_vs_LFD, "significant_HFD_vs_LFD.csv", row.names = TRUE)
write.csv(significant_HFD_inf_vs_HFD, "significant_HFD_inf_vs_HFD.csv", row.names = TRUE)
#Pie split up function
library(ggplot2)
library(cowplot)
df <- data.frame(
  group = c("Altered", "Not significantly altered"),
  n     = c(150, 350)
)
df$fraction <- df$n / sum(df$n)
df$ymax <- cumsum(df$fraction)
df$ymin <- c(0, head(df$ymax, -1))
df$y    <- (df$ymin + df$ymax) / 2
df$label <- c("Altered\n(n=150)", "Not significantly altered\n(n=350)")
cols <- c("Altered" = "#FBB042", "Not significantly altered" = "#C3996C")
p <- ggplot(df) +
  geom_rect(
    aes(ymin = ymin, ymax = ymax, xmin = 1, xmax = 2, fill = group),
    color = NA, linewidth = 2
  ) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.2) +                  
  scale_fill_manual(values = cols) +
  geom_text(aes(x = 1.45, y = y, label = label), size = 3.6) +
  theme_void() +
  theme(legend.position = "none")
ggsave("donut_plotmetabolites.png",
       plot = p,
       width = 6, height = 6, dpi = 600)
ggsave("donut_plot_metabolites.pdf",
       plot = p,
       width = 6, height = 6,
       device = cairo_pdf)  
















#Figure 3F
#Protein vs Metabolite correlation
#“For easy annotation of the metabolites, the whole universe of the significant metabolites was subset to those that can be annotated using the RefMet database.
library(RefMet)
library(readr)
library(dplyr)
metab_data <- read_csv("Significant metabolites for correlation.csv")
metab_names <- metab_data$Metabolites
# Map to RefMet
refmet_mapped <- refmet_map_df(metab_names)
write.csv(refmet_mapped, "Significant metabolites_matched RefMet_for correlation.csv", row.names = FALSE)
#Correlation analysis”
metabolites <- read.csv("Significant metabolites_matched RefMet_ for correlation.csv", sep = ";", check.names = FALSE)
proteins    <- read.csv("Significant protein for correlation.csv", sep = ";", check.names = FALSE)
sample_cols <- grep("LFD_|HFD", colnames(metabolites), value = TRUE)
met_mat <- metabolites[, c("Metabolites", "Main.class", sample_cols)]
prot_mat <- proteins[, c("Protein name", sample_cols)]
#To transpose data
met_expr <- as.data.frame(t(met_mat[, -(1:2)]))
colnames(met_expr) <- met_mat$Metabolites
rownames(met_expr) <- sample_cols

prot_expr <- as.data.frame(t(prot_mat[, -1]))
colnames(prot_expr) <- prot_mat$`Protein name`
rownames(prot_expr) <- sample_cols
# Pearson correlation 
results <- expand.grid(
  Metabolite = colnames(met_expr),
  Protein = colnames(prot_expr),
  stringsAsFactors = FALSE
)
cor_fun <- function(met, prot) {
  test <- cor.test(met_expr[[met]], prot_expr[[prot]], method = "pearson")
  c(r = unname(test$estimate),
    p = test$p.value)
}
cor_vals <- mapply(cor_fun, results$Metabolite, results$Protein)
results$r <- cor_vals["r", ]
results$p <- cor_vals["p", ]

# Multiple testing correction 
results$q <- p.adjust(results$p, method = "BH")
#Add metabolite Main.class
met_classes <- metabolites[, c("Metabolites", "Main.class")]
results <- results %>%
  left_join(met_classes, by = c("Metabolite" = "Metabolites"))
# Filter significant correlations
sig_results <- results %>%
  filter(abs(r) >= 0.6, q <= 0.05)
write.csv(results, "Metabolite_Protein_Correlations_all.csv", row.names = FALSE)
write.csv(sig_results, "Metabolite_Protein_Correlations_significant.csv", row.names = FALSE)

#Heatmap to visualize the metabolites and proteins
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
all_csv <- "Metabolite_Protein_Correlations_all.csv"
all_df <- read.csv(all_csv, check.names = FALSE)
r_mat <- all_df %>%
  select(Metabolite, Protein, r) %>%
  pivot_wider(names_from = Protein, values_from = r) %>%
  as.data.frame()
rownames(r_mat) <- r_mat$Metabolite
r_mat$Metabolite <- NULL
r_mat <- as.matrix(r_mat)
met2class <- all_df %>%
  select(Metabolite, Main.class) %>%
  distinct()
class_vec <- met2class$Main.class[ match(rownames(r_mat), met2class$Metabolite) ]
class_counts <- table(class_vec)
class_levels <- names(sort(class_counts, decreasing = TRUE))
class_fac <- factor(class_vec, levels = class_levels)
col_fun <- colorRamp2(c(-0.7, 0, 0.7), c("#2166AC", "#FFFFFF", "#B2182B"))
# Alternating strip colors for each Main.class block 
block_colors <- rep(c("black", "grey70"), length.out = length(class_levels))
block_col_vec <- setNames(block_colors, class_levels)
block_fac_col <- block_col_vec[as.character(class_fac)]
row_ha <- rowAnnotation(
  Block = block_fac_col,
  col = list(Block = c("black" = "black", "grey70" = "grey70")),
  show_annotation_name = FALSE,
  width = unit(2, "mm")
)
# Heatmap function 
ht <- Heatmap(
  r_mat,
  name = "Pearson r",
  col = col_fun,
  cluster_rows = FALSE,         
  cluster_columns = FALSE,
  show_row_names = FALSE,       
  show_column_names = FALSE,    
  row_split = class_fac,
  gap = unit(0.25, "mm"),         
  row_title_rot = 0,            
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_title = "Proteins",
  left_annotation = row_ha,     
  heatmap_legend_param = list(title = "Pearson r")
)
pdf("Heatmap_grouped_by_Mainclass.pdf", width = 10, height = 12)
draw(ht, heatmap_legend_side = "right")
dev.off()

png("Heatmap_grouped_by_Mainclass.png", width = 2000, height = 2600, res = 250)
draw(ht, heatmap_legend_side = "right")
dev.off()


#Figure 3G Tgm2 and significantly correlated metabolites
library(tidyverse)
library(ggplot2)
df <- tribble(
  ~Metabolite, ~Protein, ~r, ~p, ~q, ~Main.class,
  "L-(+)-Valine", "Tgm2", 0.772099, 9.89e-06, 0.001686, "Amino acids and peptides",
  "N-acetyl-L-2-aminoadipic acid", "Tgm2", 0.737281, 3.95e-05, 0.002989, "Amino acids and peptides",
  "L-Hexanoylcarnitine", "Tgm2", -0.769786, 1.09e-05, 0.001720, "Fatty esters",
  "Troxipide", "Tgm2", 0.767122, 1.22e-05, 0.001752, "Phenols",
  "Nicotinamide ribotide", "Tgm2", 0.672858, 3.15e-04, 0.009551, "Pyridine alkaloids",
  "alpha-Ketoglutaric acid", "Tgm2", 0.736298, 4.10e-05, 0.003042, "TCA acids"
)
df <- df %>% arrange(r)
df$Metabolite <- factor(df$Metabolite, levels = df$Metabolite)
p <- ggplot(df, aes(x = r, y = Metabolite)) +
  geom_point(aes(size = -log10(q), fill = Main.class), shape = 21, color = "black") +
  scale_fill_brewer(palette = "Set2") +
  scale_size_continuous(name = "-log10(q)", range = c(3,8)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "Pearson correlation (r)",
    y = "Metabolites",
    title = "Correlation of Tgm2 with metabolites"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10)
  )
ggsave("Tgm2_metabolite_dotplot.pdf", p, width = 7, height = 5)
ggsave("Tgm2_metabolite_dotplot.png", p, width = 7, height = 5, dpi = 30
       
