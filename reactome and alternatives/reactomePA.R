
library(tidyverse)
library(dplyr)
library(ReactomePA)
library(org.Hs.eg.db)



import_data <- read.csv("aim3_allprots_protein_dataset_DE_results.csv", check.names = F)



### pathway enrichment analysis:

# get sig genes only
sig_genes <- import_data %>%
  filter(`Buttock-Forearm adj p-value` < 0.05) %>%
  dplyr::select("Gene") %>%
  pull(Gene)

# convert to entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = sig_genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# pathway enrichment analysis
enrich_results <- enrichPathway(gene = entrez_ids,
                                organism = "human",
                                pvalueCutoff = 0.05,
                                readable = TRUE)

# Bar plot
barplot(enrich_results, showCategory = 20)

# Dot plot
dotplot(enrich_results, showCategory = 20)

# Enrichment map
emapplot(enrich_results)

# Pathway-gene relationship
cnetplot(enrich_results)






### gene set enrichment analysis GSEA (uses all genes ranked by something e.g. logFC)

ranked_genes <- import_data %>%
  dplyr::select(Gene, `Buttock-Forearm log2FC`)
gene_list <- ranked_genes$`Buttock-Forearm log2FC`
names(gene_list) <- ranked_genes$Gene

# Sort descending
gene_list <- sort(gene_list, decreasing = TRUE)


# convert
entrez_ids_ranked_genes <- mapIds(org.Hs.eg.db,
                     keys = names(gene_list),
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
valid_idx <- !is.na(entrez_ids_ranked_genes)
gene_list <- gene_list[valid_idx]
names(gene_list) <- entrez_ids_ranked_genes[valid_idx]

entrez_gene_list <- data.frame(
  GeneID = names(gene_list),
  Value = as.numeric(gene_list),
  row.names = NULL,
  stringsAsFactors = FALSE
)
# write.csv(entrez_gene_list, "entrez_gene_list.csv")


# GSEA
gsea_results <- gsePathway(geneList = gene_list,
                           organism = "human",
                           pvalueCutoff = 1,
                           verbose = FALSE)


# Dot plot of enriched pathways
dotplot(gsea_results, showCategory = 10)













