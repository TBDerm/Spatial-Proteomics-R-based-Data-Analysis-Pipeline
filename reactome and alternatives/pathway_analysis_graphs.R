
library(ggplot2)
library(dplyr)
library(ReactomeGSA)
library(stringr)


# if do gene expression analysis on reactome website, it gives u a code to use in R to create own visualisations:

all_prots_input <- "9a53a523-f3c5-11ef-8ab5-2677535b140e"

reactomeGSA_result <- get_reactome_analysis_result(analysis_id = all_prots_input, reactome_url = "https://gsa.reactome.org")
reactomeGSA_result <- reactomeGSA_result@results
reactomeGSA_result_df <- reactomeGSA_result[[1]]$pathways

top_20_pathways <- reactomeGSA_result_df %>%
  arrange(FDR) %>%  # Sort by FDR (ascending order, smaller FDR is more significant)
  head(20) 
top_20_pathways$Name_wrapped <- str_wrap(top_20_pathways$Name, width = 30)



ggplot(top_20_pathways, aes(x = reorder(Name_wrapped, NGenes), y = NGenes, size = NGenes, color = PValue)) +
  geom_point() +
  coord_flip() +  # Flip axes for better readability
  theme_bw() +
  labs(title = "Reactome Gene Expression Dotplot",
       y = "Number of genes in pathway",
       x = "Pathway",
       size = "Number of Genes") +  scale_color_gradient(low = "green", high = "black") +
  theme(axis.text.y = element_text(size = 10))



