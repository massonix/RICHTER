library(tidyverse)
library(GOstats)
library(AnnotationDbi)
library(Seurat)

richter_l <- readRDS("/Volumes/Massoni_external/B_cell_atlas/RICHTER/current/results/R_objects/richter_seurat_patients_list_clustered.rds")
richter_l <- richter_l[c("ICGC_012", "ICGC_019")]
VlnPlot(richter_l$ICGC_012, features = c("AC007952.4", "AC245014.3"), group.by = "status")

c("AC007952.4", "AC245014.3") %in% richter_l$ICGC_019[["RNA"]]@var.features

FeaturePlot(richter_l$ICGC_019, features = c("AC007952.4", "AC245014.3"))

################################################
################### ICGC_012 ###################
################################################
mat_012 <- richter_l$ICGC_012[["RNA"]]@scale.data
distance_mat_012 <- dist(mat_012, method = "euclidean")
h_clust_012 <- hclust(distance_mat_012, method = "ward.D")
plot(h_clust_012, labels = FALSE)
clusters_012 <- cutree(h_clust_012, k = 17)
table(clusters_012)
clusters_012[c("AC007952.4", "AC245014.3")]
genes_012 <- names(clusters_012[clusters_012 == 11])

target <- names(clusters_012[clusters_012 == 11])
universe <- rownames(mat_012)
library(org.Hs.eg.db)
target_entrez <- AnnotationDbi::select(
  x = org.Hs.eg.db, 
  keys = target, 
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID
target_entrez <- target_entrez[!is.na(target_entrez)]
target_entrez

universe_entrez <- AnnotationDbi::select(
  x = org.Hs.eg.db, 
  keys = universe, 
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID
universe_entrez <- universe_entrez[!is.na(universe_entrez)]
universe_entrez

params <- new("GOHyperGParams", geneIds = target_entrez, 
              universeGeneIds = universe_entrez, annotation = "org.Hs.eg.db",
              ontology = "BP", pvalueCutoff = 1, 
              conditional = TRUE, testDirection = "over")
hgOver_gt <- hyperGTest(params)

go_results <- summary(hgOver_gt)
library(tidyverse)

selection <- go_results$Size >= 3 & go_results$Size <= 600 & go_results$Count >= 5 & go_results$OddsRatio > 3 & go_results$Pvalue < 0.05
go_results_filt <- go_results[selection, ]
go_results_filt_gg <- go_results_filt %>% 
  ggplot(aes(fct_reorder(Term, OddsRatio), OddsRatio)) +
  geom_col() +
  labs(x = "", y = "Odds Ratio") +
  theme_classic() +
  coord_flip()
ggsave(
  filename = "current/results/plots/go_lncrna_012.png",
  plot = go_results_filt_gg,
  height = 14,
  width = 19,
  units = "cm"
)

go_012_lst <- mapIds(org.Hs.eg.db, go_results_filt$GOBPID, "ENTREZID", "GOALL", multiVals = "list")
go_012_lst <- purrr::map(go_012_lst, function(x) x[x %in% target_entrez])
go_012_lst_symb <- purrr::map_chr(go_012_lst, function(x) {
  symbs <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  names(symbs) <- NULL
  symbs <- unique(symbs)
  symbs <- str_c(symbs, collapse = "/")
  symbs
})
go_012_lst_symb <- go_012_lst_symb[go_results_filt$GOBPID]
go_012_lst_df <- data.frame(
  term = go_results_filt$Term, 
  enriched_genes = go_012_lst_symb
)
write_tsv(go_012_lst_df, path = "current/results/tables/gene_ontology_012_genes", col_names = TRUE, quote_escape = FALSE)

################################################
################### ICGC_019 ###################
################################################
mat_019 <- richter_l$ICGC_019[["RNA"]]@scale.data
distance_mat_019 <- dist(mat_019, method = "euclidean")
h_clust_019 <- hclust(distance_mat_019, method = "ward.D")
plot(h_clust_019, labels = FALSE)
clusters_019 <- cutree(h_clust_019, k = 22)
table(clusters_019)
clusters_019[c("AC007952.4", "AC245014.3")]
genes_019 <- names(clusters_019[clusters_019 == 5])

target_019 <- names(clusters_019[clusters_019 == 5])
universe_019 <- rownames(mat_019)
library(org.Hs.eg.db)
target_entrez_019 <- AnnotationDbi::select(
  x = org.Hs.eg.db, 
  keys = target_019, 
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID
target_entrez_019 <- target_entrez_019[!is.na(target_entrez_019)]
target_entrez_019

universe_entrez_019 <- AnnotationDbi::select(
  x = org.Hs.eg.db, 
  keys = universe_019, 
  keytype = "SYMBOL",
  columns = "ENTREZID"
)$ENTREZID
universe_entrez_019 <- universe_entrez_019[!is.na(universe_entrez_019)]
universe_entrez

params_019 <- new("GOHyperGParams", geneIds = target_entrez_019, 
              universeGeneIds = universe_entrez_019, annotation = "org.Hs.eg.db",
              ontology = "BP", pvalueCutoff = 1, 
              conditional = TRUE, testDirection = "over")
hgOver_gt_019 <- hyperGTest(params_019)
go_results_019 <- summary(hgOver_gt_019)

selection <- go_results_019$Size >= 3 & go_results_019$Size <= 600 & go_results_019$Count >= 5 & go_results_019$OddsRatio > 3 & go_results_019$Pvalue < 0.05
go_results_019_filt <- go_results_019[selection, ]
go_results_filt_gg_019 <- go_results_019_filt %>% 
  ggplot(aes(fct_reorder(Term, OddsRatio), OddsRatio)) +
  geom_col() +
  labs(x = "", y = "Odds Ratio") +
  theme_classic() +
  coord_flip()
ggsave(
  filename = "current/results/plots/go_lncrna_019.png",
  plot = go_results_filt_gg_019,
  height = 14,
  width = 19,
  units = "cm"
)

go_019_lst <- mapIds(org.Hs.eg.db, go_results_019_filt$GOBPID, "ENTREZID", "GOALL", multiVals = "list")
go_019_lst <- purrr::map(go_019_lst, function(x) x[x %in% target_entrez_019])
go_019_lst_symb <- purrr::map_chr(go_019_lst, function(x) {
  symbs <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  names(symbs) <- NULL
  symbs <- unique(symbs)
  symbs <- str_c(symbs, collapse = "/")
  symbs
})
go_019_lst_symb <- go_019_lst_symb[go_results_019_filt$GOBPID]
go_019_lst_df <- data.frame(
  term = go_results_019_filt$Term, 
  enriched_genes = go_019_lst_symb
)
write_tsv(go_019_lst_df, path = "current/results/tables/gene_ontology_019_genes", col_names = TRUE, quote_escape = FALSE)
