library(tidyverse)
library(GOstats)
library(AnnotationDbi)
library(Seurat)
library(viridis)
library(VennDiagram)
library(ggrepel)

richter_l <- readRDS("/Volumes/Massoni_external/B_cell_atlas/RICHTER/current/results/R_objects/richter_seurat_patients_list_clustered.rds")
richter_l <- richter_l[c("ICGC_012", "ICGC_019", "ICGC_3299")]

Idents(richter_l$ICGC_012) <- "annotation"
richter_l$ICGC_012 <- subset(richter_l$ICGC_012, idents = c("CXCR4+", "CXCR4-"))
Idents(richter_l$ICGC_019) <- "annotation"
richter_l$ICGC_019 <- subset(richter_l$ICGC_019, idents = c("CXCR4+", "CXCR4-"))
Idents(richter_l$ICGC_3299) <- "annotation"
richter_l$ICGC_3299 <- subset(richter_l$ICGC_3299, idents = c("CXCR4+", "CXCR4-"))
richter_l <- purrr::map(richter_l, function(seurat) {
  seurat <- seurat %>%
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()
})
cxcr4_expression <- purrr::map2(richter_l, names(richter_l), function(seurat, title) {
  cxcr4 <- as.vector(seurat[["RNA"]]@data["CXCR4", ])
  cxcr4 <- data.frame(cxcr4 = cxcr4)
  ggplot(cxcr4, aes(cxcr4)) +
    geom_histogram() +
    geom_vline(xintercept = 1.5, color = "darkblue", linetype = "dashed") +
    geom_vline(xintercept = 3.5, color = "darkblue", linetype = "dashed") +
    labs(title = title, x = "CXCR4 expression") +
    theme_classic() +
    theme(plot.title = element_text(size = 13, hjust = 0.5))
})
ggpubr::ggarrange(plotlist = cxcr4_expression, ncol = 3)
richter_l <- purrr::map(richter_l, function(seurat) {
  cxcr4 <- as.vector(seurat[["RNA"]]@data["CXCR4", ])
  cxcr4_status <- case_when(
    cxcr4 < 1.5 ~ "low",
    cxcr4 >= 1.5 & cxcr4 <= 3.5 ~ "mid",
    cxcr4 > 3.5 ~ "high"
  )
  seurat$cxcr4_status <- factor(cxcr4_status, levels = c("low", "mid", "high"))
  Idents(seurat) <- "cxcr4_status"
  seurat <- subset(seurat, subset = cxcr4_status %in% c("low", "high"))
  seurat
})
cxcr4_expression2 <- purrr::map2(richter_l, names(richter_l), function(seurat, title) {
  cxcr4 <- as.vector(seurat[["RNA"]]@data["CXCR4", ])
  cxcr4 <- data.frame(cxcr4 = cxcr4)
  ggplot(cxcr4, aes(cxcr4)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = 1.5, color = "darkblue", linetype = "dashed") +
    geom_vline(xintercept = 3.5, color = "darkblue", linetype = "dashed") +
    labs(title = title, x = "CXCR4 expression") +
    theme_classic() +
    theme(plot.title = element_text(size = 13, hjust = 0.5))
})
ggpubr::ggarrange(plotlist = cxcr4_expression2, ncol = 3)
dea_cxcr4_status <- purrr::map(richter_l, function(seurat) {
  Idents(seurat) <- "cxcr4_status"
  df <- FindMarkers(
    seurat,
    ident.1 = "high",
    ident.2 = "low",
    test.use = "wilcox"
  )
  df
})

volcano_plots <- purrr::map2(dea_cxcr4_status, names(dea_cxcr4_status), function(df, title) {
  df <- df %>%
    rownames_to_column(var = "gene") %>% 
    dplyr::mutate(
      log_10_adj_p = -1 * log10(p_val_adj),
      significance = ifelse(p_val_adj < 0.001, TRUE, FALSE)
    )
  subset_df <- df %>%
    filter(abs(avg_logFC) > 0.5, significance == TRUE)
  p <- df %>% 
    ggplot(aes(avg_logFC, log_10_adj_p, color = significance)) +
      geom_point() +
      geom_text_repel(data = subset_df, aes(label = gene), color = "black", size = 3) +
      scale_color_manual("significance", values = c("darkgrey", "green")) +
      labs(
        title = title,
        x = "log (CXCR4+ / CXCR4-)",
        y = "-log10 adjusted p-value"
      ) +
      theme_classic() +
      theme(plot.title = element_text(size = 13, hjust = 0.5))
  p
})    
ggpubr::ggarrange(plotlist = volcano_plots, ncol = 3)



up_in_high_cxcr4 <- purrr::map(dea_cxcr4_status, ~ rownames(.x)[.x$avg_logFC > 0])
down_in_high_cxcr4 <- purrr::map(dea_cxcr4_status, ~ rownames(.x)[.x$avg_logFC < 0])

vd_high <- VennDiagram::venn.diagram(up_in_high_cxcr4, filename = NULL, fill = 2:4)
grid.draw(vd_high)

Reduce(intersect, up_in_high_cxcr4)


# GO
library(org.Hs.eg.db)
up_in_high_cxcr4_inters <- Reduce(intersect, up_in_high_cxcr4[c("ICGC_012", "ICGC_019")])
target_entrez <- AnnotationDbi::select(
  x = org.Hs.eg.db, 
  keys = up_in_high_cxcr4_inters, 
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID
target_entrez <- target_entrez[!is.na(target_entrez)]
target_entrez

universe_entrez <- AnnotationDbi::select(
  x = org.Hs.eg.db, 
  keys = rownames(richter_l$ICGC_012), 
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID
universe_entrez <- universe_entrez[!is.na(universe_entrez)]
universe_entrez
library(GOstats)
params <- new("GOHyperGParams", geneIds = target_entrez, 
              universeGeneIds = universe_entrez, annotation = "org.Hs.eg.db",
              ontology = "BP", pvalueCutoff = 1, 
              conditional = TRUE, testDirection = "over")
hgOver_gt <- hyperGTest(params)
go_results <- summary(hgOver_gt)

selection <- go_results$Size >= 3 & go_results$Size <= 600 & go_results$Count >= 4 & go_results$OddsRatio > 3 & go_results$Pvalue < 0.05
go_results_filt <- go_results[selection, ]
go_results_filt_gg <- go_results_filt %>% 
  ggplot(aes(fct_reorder(Term, OddsRatio), OddsRatio)) +
  geom_col() +
  labs(x = "", y = "Odds Ratio") +
  theme_classic() +
  coord_flip()

go_lst <- mapIds(org.Hs.eg.db, go_results_filt$GOBPID, "ENTREZID", "GOALL", multiVals = "list")
go_lst <- purrr::map(go_lst, function(x) x[x %in% target_entrez])
go_lst_symb <- purrr::map_chr(go_lst, function(x) {
  symbs <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  names(symbs) <- NULL
  symbs <- unique(symbs)
  symbs <- str_c(symbs, collapse = "/")
  symbs
})
go_lst_symb <- go_lst_symb[go_results_filt$GOBPID]
go_lst_df <- data.frame(
  term = go_results_filt$Term, 
  enriched_genes = go_lst_symb
)
write_tsv(go_lst_df, path = "/Volumes/Massoni_external/B_cell_atlas/RICHTER/current/results/tables/gene_ontology_cxcr4_hi_vs_low.tsv", col_names = TRUE, quote_escape = FALSE)



# CIRBP
seurat <- richter_l$ICGC_012
title <- "ICGC_012"
cirbp_l <- purrr::map2(richter_l, names(richter_l), function(seurat, title) {
  p <- VlnPlot(seurat, features = "CIRBP", group.by = "cxcr4_status")
  p +
    scale_x_discrete("", labels = c("CXCR4 low", "CXCR4 high")) +
    labs(y = "CIRBP Expression Level", title = title) +
    theme(legend.position = "none")
})
ggpubr::ggarrange(plotlist = cirbp_l, ncol = 3)


# Save excel files
dea_cxcr4_status <- purrr::map(dea_cxcr4_status, rownames_to_column, var = "gene")
openxlsx::write.xlsx(dea_cxcr4_status, file = "/Volumes/Massoni_external/B_cell_atlas/RICHTER/current/results/tables/differential_expression_cxcr4_hi_vs_low.xlsx")


vd_low <- VennDiagram::venn.diagram(down_in_high_cxcr4, filename = NULL, fill = 2:4)
grid.draw(vd_low)
Reduce(intersect, down_in_high_cxcr4)






# pca_gg <- purrr::map(richter_l, function(seurat) {
#   plots <- purrr::map(c("CXCR4", "AC007952.4", "AC245014.3"), function(x) {
#     df <- seurat@reductions$pca@cell.embeddings %>%
#       as.data.frame() %>% 
#       dplyr::mutate(expression = seurat[["RNA"]]@data[x, ])
#     p <- df %>% 
#       ggplot(aes(PC_1, PC_2, color = expression)) +
#         geom_point() +
#         scale_color_viridis("Expression") +
#         labs(title = x, x = "PC1", "PC2") +
#         theme_classic() +
#         theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
#     p
#   })
#   plots
# })
# pca_gg_arranged <- purrr::map(pca_gg, function(x) {
#   ggpubr::ggarrange(plotlist = x, ncol = 1, nrow = 3)
# })
# VizDimLoadings(richter_l$ICGC_012, dims = 1, balanced = TRUE)
# 
# 
# 
# c("AC007952.4", "AC245014.3") %in% richter_l$ICGC_019[["RNA"]]@var.features
# FeaturePlot(richter_l$ICGC_019, features = c("AC007952.4", "AC245014.3"))
# 
# ################################################
# ################### ICGC_012 ###################
# ################################################
# mat_012 <- richter_l$ICGC_012[["RNA"]]@scale.data
# distance_mat_012 <- dist(mat_012, method = "euclidean")
# h_clust_012 <- hclust(distance_mat_012, method = "ward.D")
# plot(h_clust_012, labels = FALSE)
# clusters_012 <- cutree(h_clust_012, k = 17)
# table(clusters_012)
# clusters_012[c("AC007952.4", "AC245014.3")]
# genes_012 <- names(clusters_012[clusters_012 == 11])
# 
# target <- names(clusters_012[clusters_012 == 11])
# universe <- rownames(mat_012)
# library(org.Hs.eg.db)
# target_entrez <- AnnotationDbi::select(
#   x = org.Hs.eg.db, 
#   keys = target, 
#   keytype = "SYMBOL",
#   columns = "ENTREZID"
# )$ENTREZID
# target_entrez <- target_entrez[!is.na(target_entrez)]
# target_entrez
# 
# universe_entrez <- AnnotationDbi::select(
#   x = org.Hs.eg.db, 
#   keys = universe, 
#   keytype = "SYMBOL",
#   columns = "ENTREZID"
# )$ENTREZID
# universe_entrez <- universe_entrez[!is.na(universe_entrez)]
# universe_entrez
# 
# params <- new("GOHyperGParams", geneIds = target_entrez, 
#               universeGeneIds = universe_entrez, annotation = "org.Hs.eg.db",
#               ontology = "BP", pvalueCutoff = 1, 
#               conditional = TRUE, testDirection = "over")
# hgOver_gt <- hyperGTest(params)
# 
# go_results <- summary(hgOver_gt)
# library(tidyverse)
# 
# selection <- go_results$Size >= 3 & go_results$Size <= 600 & go_results$Count >= 5 & go_results$OddsRatio > 3 & go_results$Pvalue < 0.05
# go_results_filt <- go_results[selection, ]
# go_results_filt_gg <- go_results_filt %>% 
#   ggplot(aes(fct_reorder(Term, OddsRatio), OddsRatio)) +
#   geom_col() +
#   labs(x = "", y = "Odds Ratio") +
#   theme_classic() +
#   coord_flip()
# ggsave(
#   filename = "current/results/plots/go_lncrna_012.png",
#   plot = go_results_filt_gg,
#   height = 14,
#   width = 19,
#   units = "cm"
# )
# 
# go_012_lst <- mapIds(org.Hs.eg.db, go_results_filt$GOBPID, "ENTREZID", "GOALL", multiVals = "list")
# go_012_lst <- purrr::map(go_012_lst, function(x) x[x %in% target_entrez])
# go_012_lst_symb <- purrr::map_chr(go_012_lst, function(x) {
#   symbs <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
#   names(symbs) <- NULL
#   symbs <- unique(symbs)
#   symbs <- str_c(symbs, collapse = "/")
#   symbs
# })
# go_012_lst_symb <- go_012_lst_symb[go_results_filt$GOBPID]
# go_012_lst_df <- data.frame(
#   term = go_results_filt$Term, 
#   enriched_genes = go_012_lst_symb
# )
# write_tsv(go_012_lst_df, path = "current/results/tables/gene_ontology_012_genes", col_names = TRUE, quote_escape = FALSE)
# 
# ################################################
# ################### ICGC_019 ###################
# ################################################
# # mat_019 <- richter_l$ICGC_019[["RNA"]]@scale.data
# # distance_mat_019 <- dist(mat_019, method = "manhattan")
# # h_clust_019 <- hclust(distance_mat_019, method = "ward.D")
# # plot(h_clust_019, labels = FALSE)
# # clusters_019 <- cutree(h_clust_019, k = 23)
# # table(clusters_019)
# # clusters_019[c("AC007952.4", "AC245014.3")]
# # genes_019 <- names(clusters_019[clusters_019 == 19])
# # 
# # target_019 <- names(clusters_019[clusters_019 == 5])
# # universe_019 <- rownames(mat_019)
# # library(org.Hs.eg.db)
# # target_entrez_019 <- AnnotationDbi::select(
# #   x = org.Hs.eg.db, 
# #   keys = target_019, 
# #   keytype = "SYMBOL",
# #   columns = "ENTREZID"
# # )$ENTREZID
# # target_entrez_019 <- target_entrez_019[!is.na(target_entrez_019)]
# # target_entrez_019
# # 
# # universe_entrez_019 <- AnnotationDbi::select(
# #   x = org.Hs.eg.db, 
# #   keys = universe_019, 
# #   keytype = "SYMBOL",
# #   columns = "ENTREZID"
# # )$ENTREZID
# # universe_entrez_019 <- universe_entrez_019[!is.na(universe_entrez_019)]
# # universe_entrez
# # 
# # params_019 <- new("GOHyperGParams", geneIds = target_entrez_019, 
# #                   universeGeneIds = universe_entrez_019, annotation = "org.Hs.eg.db",
# #                   ontology = "BP", pvalueCutoff = 1, 
# #                   conditional = TRUE, testDirection = "over")
# # hgOver_gt_019 <- hyperGTest(params_019)
# # go_results_019 <- summary(hgOver_gt_019)
# # 
# # selection <- go_results_019$Size >= 3 & go_results_019$Size <= 600 & go_results_019$Count >= 5 & go_results_019$OddsRatio > 3 & go_results_019$Pvalue < 0.05
# # go_results_019_filt <- go_results_019[selection, ]
# # go_results_filt_gg_019 <- go_results_019_filt %>% 
# #   ggplot(aes(fct_reorder(Term, OddsRatio), OddsRatio)) +
# #   geom_col() +
# #   labs(x = "", y = "Odds Ratio") +
# #   theme_classic() +
# #   coord_flip()
# # ggsave(
# #   filename = "current/results/plots/go_lncrna_019.png",
# #   plot = go_results_filt_gg_019,
# #   height = 14,
# #   width = 19,
# #   units = "cm"
# # )
# # 
# # go_019_lst <- mapIds(org.Hs.eg.db, go_results_019_filt$GOBPID, "ENTREZID", "GOALL", multiVals = "list")
# # go_019_lst <- purrr::map(go_019_lst, function(x) x[x %in% target_entrez_019])
# # go_019_lst_symb <- purrr::map_chr(go_019_lst, function(x) {
# #   symbs <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
# #   names(symbs) <- NULL
# #   symbs <- unique(symbs)
# #   symbs <- str_c(symbs, collapse = "/")
# #   symbs
# # })
# # go_019_lst_symb <- go_019_lst_symb[go_results_019_filt$GOBPID]
# # go_019_lst_df <- data.frame(
# #   term = go_results_019_filt$Term, 
# #   enriched_genes = go_019_lst_symb
# # )
# # write_tsv(go_019_lst_df, path = "current/results/tables/gene_ontology_019_genes", col_names = TRUE, quote_escape = FALSE)
# # 
# # 
# # 
# # pca_gg <- purrr::map(richter_l, function(seurat) {
# #   plots <- purrr::map(c("CXCR4", "AC007952.4", "AC245014.3"), function(x) {
# #     df <- seurat@reductions$pca@cell.embeddings %>%
# #       as.data.frame() %>% 
# #       dplyr::mutate(expression = seurat@data[x, ])
# #     p <- df %>% 
# #       ggplot(aes(PC_1, PC_2, color = expression)) +
# #       geom_point() +
# #       scale_color_viridis("Expression") +
# #       labs(title = x, x = "PC1", "PC2") +
# #       theme_classic() +
# #       theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
# #     p
# #   })
# #   plots
# # })
# # 
# # 
# # 
# # # pca_gg <- purrr::map(richter_l, function(seurat) {
# # seurat <- richter_l[[1]]
# # # plots <- purrr::map(c("CXCR4", "AC007952.4", "AC245014.3"), function(x) {
# # x <- "CXCR4"
# # df <- seurat@reductions$pca@cell.embeddings %>%
# #   as.data.frame() %>% 
# #   dplyr::mutate(expression = seurat[["RNA"]]@data[x, ])
# # p <- df %>% 
# #   ggplot(aes(PC_1, PC_2, color = expression)) +
# #   geom_point() +
# #   scale_color_viridis("Expression") +
# #   labs(title = x, x = "PC1", "PC2") +
# #   theme_classic() +
# #   theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
# # p
# # # })
# # plots
# # # })
# 
# 
