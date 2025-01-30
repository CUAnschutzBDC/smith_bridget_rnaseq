library(CIBERSORT)
library(here)


human_mouse <- read.csv(file.path("/pl/active/Anschutz_BDC/", 
                                  "resources/single_cell_references", 
                                  "gene_lists", 
                                  "human_mouse_gene_mapping_20221003.csv"))

counts <- read.table(here("results/bridget_rnaseq_countTable_cutadapt_trim.txt"),
                     header = TRUE)

sample_info <- read.table(here("files/sample_info.csv"), header = TRUE,
                          sep = ",")

colnames(counts) <- gsub(".*featureCount_cutadapt_trim.", "", colnames(counts))
counts$gene_id <- gsub("_ENS.*", "", counts$gene)
rownames(counts) <- make.unique(counts$gene_id)
counts$gene <- NULL
counts$gene_id <- NULL

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ genotype)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

dds <- DESeq(dds)

norm_data <- counts(dds, normalized = TRUE)

log_norm_data <- log1p(norm_data)

sig_matrix <- read.table(here("files/LM22.txt"), sep = "\t", header = TRUE,
                         row.names = 1, check.names = FALSE)

#Run CIBERSORT abs 
#The number of permutation
cibersort_perm = 100
#Quantile normalization of input mixture, default = FALSE for RNA-Seq data
cibersort_qn = FALSE

res_ciber <- cibersort(data.matrix(sig_matrix), data.matrix(log_norm_data), perm = 10, 
                       QN = cibersort_qn)


final_res <- res_ciber %>%
  data.frame() %>%
  dplyr::select(!c(P.value, Correlation, RMSE)) %>%
  tibble::rownames_to_column("sample") %>%
  tidyr::pivot_longer(cols = colnames(.)[2:ncol(.)],
                      names_to = "cell_type",
                      values_to = "fraction") %>%
  dplyr::mutate(Sample_name = gsub("X", "", sample)) %>%
  merge(sample_info, by = "Sample_name")


colors <- colorRampPalette(
  colors = RColorBrewer::brewer.pal(
    n = 9, name = "Set1"
  )
)(length(unique(final_res$cell_type)))

p1 <- ggplot2::ggplot(final_res, ggplot2::aes(x = Sample_name, y = fraction,
                                        fill = cell_type)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = colors) +
  ggplot2::theme_classic()

pdf(file.path("results/R_analysis/images/cibersort_cell_types.pdf"),
    width = 10, height = 10)
print(p1)
dev.off()
