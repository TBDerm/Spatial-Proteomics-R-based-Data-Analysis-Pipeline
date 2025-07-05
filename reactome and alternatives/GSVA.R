
library(tidyverse)
library(GSVA)
library(dplyr)
library(ggfortify)
library(limma)
library(pheatmap)
library(ggrepel)



ECM_data <- read.csv("aim3_ECMprots_protein_dataset_DE_results.csv", check.names = F)


# import_data <- read.csv("aim3_allprots_protein_dataset_DE_results.csv", check.names = F)
# panther_class <- read.csv("aim3_panther_classifications.csv", check.names = F)
# import_data <- import_data %>%
#   left_join(panther_class, by = "Gene") %>%
#   filter(!duplicated(Gene))

# make a matrix of genes x samples, with genes as rownames
gene_sample_matrix <- ECM_data %>% 
  dplyr::select("0841_Buttock", "0842_Buttock", "0863_Buttock", "0841_Forearm", "0842_Forearm", "0863_Forearm") %>%
  as.matrix()
rownames(gene_sample_matrix) <- ECM_data$Gene

# 
# # # create gene sets from category data - I've put all categories into one col (if a protein is present in >1 cat, that is preserved)
#  gene_sets <- import_data %>%
#    select(Gene, Matrisome_Category, BM_Status, panther_class) %>%
#    mutate(panther_class = str_to_title(panther_class)) %>%
#    separate(panther_class, into = c("panther_class", "other"), sep = "\\(") %>%
#    pivot_longer(cols = c(Matrisome_Category, BM_Status, panther_class),
#                 names_to = "cat",
#                 values_to = "protein_category") %>%
#    mutate(
#      protein_category = case_when(
#        protein_category == "Core" ~ "BM Core Protein",
#        protein_category == "Confirmed" ~ "BM Confirmed Protein",
#        TRUE ~ protein_category
#      )
#    ) %>%
#    filter(!is.na(protein_category)) %>%
#    select(-cat, -other)

# write.csv(gene_sets, "all_categories.csv", row.names = F)


gene_sets <- ECM_data %>%
  select(Gene, Matrisome_Category, BM_Status) %>%
  pivot_longer(cols = c(Matrisome_Category, BM_Status),
               names_to = "cat",
               values_to = "protein_category") %>%
  mutate(
    protein_category = case_when(
      protein_category == "Core" ~ "BM Core Protein",
      protein_category == "Confirmed" ~ "BM Confirmed Protein",
      TRUE ~ protein_category
    )
  ) %>%
  filter(!is.na(protein_category)) %>%
  dplyr::select(-cat)





# # (not using) create gene sets from category data - I've put all categories into one col (only one cat per gene, prioritising BM cats)
# gene_sets2 <- import_data %>%
#   select(Gene, Matrisome_Category, BM_Status) %>%
#   mutate(protein_category = ifelse(is.na(BM_Status), Matrisome_Category, BM_Status)) %>%  # if BM is NA, takes matrisome cat value. Otherwise takes BM value
#   mutate(protein_category = replace_na(protein_category, "Intracellular Protein")) %>%
#   mutate(
#     protein_category = case_when(
#       protein_category == "Core" ~ "BM Core Protein",
#       protein_category == "Confirmed" ~ "BM Confirmed Protein",
#       TRUE ~ protein_category
#     )
#   ) %>%
#   select(Gene, protein_category)



# makes it into a list where it counts and divides genes according to category
gene_sets_list <- split(gene_sets$Gene, gene_sets$protein_category)


# pass the matrix and gene sets into the param function
gsvaPar <- gsvaParam(gene_sample_matrix, gene_sets_list, minSize = 6)   # minimum num of genes in category to include that cat in analysis
gsvaPar

# then generate GSVA results using the param object - makes a matrix of pathway scores ('pathways' being the gene sets)
gsva_results <- gsva(gsvaPar, verbose=FALSE)
dim(gsva_results)




### then stat tests / vis to compare experimental groups

# DE analysis on GVSA results
group <- factor(c(rep("Buttock", 3), rep("Forearm", 3)))
design <- model.matrix(~ group)
# Fit model to GSVA results
fit <- lmFit(gsva_results, design)
fit <- eBayes(fit)
# View differentially enriched pathways
topTable(fit, coef = 2)
# makeing dfs of results
protCategory_limma_res_top10 <- topTable(fit, coef = 2)
protCategory_limma_res_top10 <- protCategory_limma_res_top10 %>%
  mutate(Gene_Set = rownames(protCategory_limma_res_top10)) %>%
  select(Gene_Set, everything())
protCategory_limma_res_all <- topTable(fit, coef = 2, number = Inf) 
protCategory_limma_res_all <- protCategory_limma_res_all %>%
  mutate(Gene_Set = rownames(protCategory_limma_res_all)) %>%
  select(Gene_Set, everything())

# write.csv(protCategory_limma_res_all, "GSVA results ECM-only.csv", row.names = F)


# view top 10 DE results in a barchart
protCategory_limma_res_all %>%
  mutate(Name = rownames(.)) %>%
  ggplot(aes(x = reorder(Name, logFC), y = logFC, fill = adj.P.Val)) +
  geom_col() +
  geom_col(data = . %>% filter(adj.P.Val < 0.1), 
           aes(x = reorder(Name, logFC), y = logFC),
           fill = NA, colour = "red3", linewidth = 2) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 18, colour = "black"), 
        axis.text.x = element_text(size = 17),
        axis.title = element_text(size = 19, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.75, face = "bold"),
        legend.text = element_text(size = 17),
        legend.title = element_text(size = 17, face = "bold")) +
  labs(title = "Differentially Affected ECM Protein Categories: Photoaged DEJ",
       x = "Protein Category", y = "logFC") +
  scale_fill_gradientn(colours = colorRampPalette(c("green3", "black"))(100), 
                       name = "P-adj", na.value = "white")



pdf("aim3_GSVA_larger.pdf", width = 10, height = 5)

dev.off()



# or dotplot
ggplot(protCategory_limma_res_top10, aes(x = reorder(Gene_Set, logFC), y = logFC, size = logFC, color = adj.P.Val)) +
  geom_point() +
  coord_flip() +
  theme_bw() +
  labs(title = "Differentially Enriched Protein Categories: Photoaged DEJ",
       y = "log2 ( Fold Change )",
       x = "Protein Category") +
  scale_color_gradient(low = "green", high = "black") +
  theme(axis.text.y = element_text(size = 10))


# another dotplot
protCategory_limma_res_top10 %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), size = abs(logFC), color = adj.P.Val)) +
  geom_point() +
  scale_color_gradient(low = "green3", high = "black") +
  theme_bw() +
  labs(title = "Differentially Enriched Protein Categories: Photoaged DEJ",
       x = "logFC",
       y = "-log10(P-adj)", color = "P-adj", size = "logFC") +
  # geom_text(data = protCategory_limma_res_top10,
  #           aes(label = Gene_Set), vjust = 1.5, size = 3, check_overlap = TRUE) +
  geom_text_repel(data = protCategory_limma_res_top10,
                  aes(label = Gene_Set), colour = "black", size = 3,
                  max.overlaps = 20)




# PCA on GVSA results
pca <- prcomp(t(gsva_results), scale. = TRUE)

autoplot(pca, x = 1, y = 2,
         data = data.frame(Group = group), 
         colour = "Group",
         size = 4) +
  scale_colour_manual(values = c("magenta", "blue")) +
  labs(title = "PCA of GSVA Enrichment Scores") +
  theme_bw()



# heatmap? idk
annotation <- data.frame(Group = group)
rownames(annotation) <- colnames(gsva_results)  # 👈 match sample names
pheatmap(gsva_results,
         annotation_col = annotation,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         main = "GSVA Enrichment Heatmap")






dev.off()








library(igraph)

# Correlation matrix of GSVA pathways
cor_mat <- cor(t(gsva_results), method = "pearson")

# Build edges for high-correlation pairs
cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
cor_df <- as.data.frame(as.table(cor_mat)) %>%
  filter(!is.na(Freq) & abs(Freq) > 0.8)  # Adjust threshold

colnames(cor_df) <- c("from", "to", "weight")
g <- graph_from_data_frame(cor_df, directed = FALSE)

plot(g, vertex.label.cex = 0.8, edge.width = abs(cor_df$weight) * 5)
















