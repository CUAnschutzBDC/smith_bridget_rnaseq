library(here)
library(DESeq2)
library(tidyverse)
library(MetBrewer)
library(pheatmap)
source(here("src/scripts/functions.R"))

save_dir <- here("results/R_analysis/images/figures")
ifelse(!dir.exists(save_dir), dir.create(save_dir, recursive = TRUE), FALSE)

dds <- readRDS(here("results/R_analysis/objs/dds_obj.rda"))

vsd <- readRDS(here("results/R_analysis/objs/vsd_obj.rda"))

de_res <- read.table(here("results/R_analysis/DE_files/KO_Ctl_all.csv"),
                     sep = ",", header = TRUE)

blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")

genotype_colors <- met.brewer(palette_name = "Derain",
                              n = length(unique(dds$genotype)),
                              type = "continuous")
names(genotype_colors) <- unique(dds$genotype)

# Heatmap of all DE genes ------------------------------------------------------
sig_genes <- de_res %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 0)

sig_genes$plot_name <- paste(sig_genes$gene_name, sig_genes$ens_id,
                             sep = "_")

heatmap <- make_heatmap(dds = dds,
                        vsd = vsd,
                        de_genes = sig_genes$plot_name,
                        treatment = "KO",
                        control = "Control",
                        group = "genotype",
                        print_genenames = FALSE,
                        cluster_cols = FALSE,
                        save_heatmap = TRUE,
                        output_dir = here("results/R_analysis/images"),
                        color_test = genotype_colors,
                        plot_groups = c("KO", "Control"),
                        save_name = "DE_heatmap")

graphics.off()
pdf(file.path(save_dir, "DE_heatmap.pdf"))
print(heatmap)
dev.off()

# Volcano plot of selected genes -----------------------------------------------
selected_genes <- c(
  "Socs3", "Tbx21", "Socs1", "Xylt1", "Gbp7", "Ighg2c", "Ddx25", "Plekhg3",
  "Iigp1", "Zdhhc23", "Ly6c2", "Siglecg", "Socs2", "Zeb2", "Itgb1", "Rgs1",
  "Ighg2b", "Stat1", "Ighm", "Zbtb32", "Cxcr3", "Fas", "Fcrl5", "Irf8", "Irf1",
  "Stat3"
)

de_res$color <- ifelse(de_res$padj < 0.05 & abs(de_res$log2FoldChange) > 0.5,
                       "red", "black")

de_res$label <- ifelse(de_res$gene %in% selected_genes, 
                       de_res$gene_name, 
                       "") 
p_all <- ggplot2::ggplot(data = de_res, mapping = ggplot2::aes(x = log2FoldChange,
                                                      y = -log(padj, 10),
                                                      color = color)) +
  ggplot2::geom_point() +
  ggplot2::xlim(c(-4,4)) +
  ggplot2::scale_color_manual(values = c("black" = "black", "red" = "red")) +
  ggplot2::geom_vline(xintercept = 0.5, linetype = "dotted") +
  ggplot2::geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  ggplot2::geom_vline(xintercept = -0.5, linetype = "dotted") +
  ggrepel::geom_text_repel(ggplot2::aes(label = label),
                           max.overlaps = Inf, box.padding = 0.2) +
  ggplot2::theme_classic()
  

sub_selected_genes <- selected_genes <- c(
  "Socs3", "Tbx21", "Socs1", "Ighg2c", "Ddx25",
  "Iigp1", "Ly6c2", "Siglecg", "Socs2", "Rgs1",
  "Ighg2b", "Stat1", "Ighm", "Cxcr3", "Fas", "Fcrl5", "Irf8", "Irf1",
  "Stat3"
)

de_res$label <- ifelse(de_res$gene %in% sub_selected_genes, 
                       de_res$gene_name, 
                       "") 

p_sub <- ggplot2::ggplot(data = de_res, mapping = ggplot2::aes(x = log2FoldChange,
                                                      y = -log(padj, 10),
                                                      color = color)) +
  ggplot2::geom_point() +
  ggplot2::xlim(c(-4,4)) +
  ggplot2::scale_color_manual(values = c("black" = "black", "red" = "red")) +
  ggplot2::geom_vline(xintercept = 0.5, linetype = "dotted") +
  ggplot2::geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  ggplot2::geom_vline(xintercept = -0.5, linetype = "dotted") +
  ggrepel::geom_text_repel(ggplot2::aes(label = label),
                           max.overlaps = Inf, box.padding = 0.2) +
  ggplot2::theme_classic()

# Recommend removing Siglecg, Fcrl5, Irf8, Stat3, Stat1, Ighm, Irf1

de_res$label <- ifelse(de_res$gene %in% sub_selected_genes &
                         abs(de_res$log2FoldChange) > 0.5, 
                       de_res$gene_name, 
                       "") 

p_sig <- ggplot2::ggplot(data = de_res, mapping = ggplot2::aes(x = log2FoldChange,
                                                      y = -log(padj, 10),
                                                      color = color)) +
  ggplot2::geom_point() +
  ggplot2::xlim(c(-4,4)) +
  ggplot2::scale_color_manual(values = c("black" = "black", "red" = "red")) +
  ggplot2::geom_vline(xintercept = 0.5, linetype = "dotted") +
  ggplot2::geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  ggplot2::geom_vline(xintercept = -0.5, linetype = "dotted") +
  ggrepel::geom_text_repel(ggplot2::aes(label = label),
                           max.overlaps = Inf, box.padding = 0.3) +
  ggplot2::theme_classic()

pdf(file.path(save_dir, "DE_volcano_all.pdf"), width = 6, height = 4)
print(p_all)
dev.off()

pdf(file.path(save_dir, "DE_volcano_subset.pdf"), width = 6, height = 4)
print(p_sub)
dev.off()

pdf(file.path(save_dir, "DE_volcano_sig.pdf"), width = 6, height = 4)
print(p_sig)
dev.off()

# GO plots ---------------------------------------------------------------------
go_res <- read.table(here("results/R_analysis/GSE_files/KO_Ctl.csv"),
                     sep = ",",
                     header = TRUE)

terms <- c(
  "Metabolic pathways", 
  "MAPK signaling pathway",
  "Growth hormone synthesis, secretion and action",
  "Toxoplasmosis",
  "JAK-STAT signaling pathway",
  "Chemokine signaling pathway"
)

full_terms <- c(
  "Th1 and Th2 cell differentiation",
  "Th17 cell differentiation",
  "T cell receptor signaling pathway",
  "B cell receptor signaling pathway"
)

go_res$log_padj <- -log10(go_res$p_value)

short_res <- go_res %>%
  dplyr::filter(source == "KEGG", term_name %in% terms)

short_res <- short_res %>%
  dplyr::arrange(precision)

short_res$term_name <- factor(short_res$term_name,
                              levels = short_res$term_name)

plot_short <- ggplot2::ggplot(short_res, ggplot2::aes(x = precision,
                                                   y = term_name,
                                                   color = log_padj,
                                                   size = intersection_size)) +
  ggplot2::geom_point() + theme_classic() +
  viridis::scale_color_viridis() +
  theme(text = ggplot2::element_text(size = 10)) +
  labs(color = "-log10(p-value)") +
  xlab("Precision (proportion of genes)") +
  ylab("Term")

full_res <- go_res %>%
  dplyr::filter(source == "KEGG", term_name %in% c(terms, full_terms))

full_res <- full_res %>%
  dplyr::arrange(precision)

full_res$term_name <- factor(full_res$term_name,
                              levels = full_res$term_name)

plot_full <- ggplot2::ggplot(full_res, ggplot2::aes(x = precision,
                                                      y = term_name,
                                                      color = log_padj,
                                                      size = intersection_size)) +
  ggplot2::geom_point() + theme_classic() +
  viridis::scale_color_viridis() +
  theme(text = ggplot2::element_text(size = 10)) +
  labs(color = "-log10(p-value)") +
  xlab("Precision (proportion of genes)") +
  ylab("Term")

pdf(file.path(save_dir, "kegg_short.pdf"), width = 10, height = 6)
print(plot_short)
dev.off()

pdf(file.path(save_dir, "kegg_full.pdf"), width = 10, height = 6)
print(plot_full)
dev.off()
