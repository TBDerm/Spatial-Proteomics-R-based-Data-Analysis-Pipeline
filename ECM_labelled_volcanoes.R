
library(tidyverse)
library(data.table)
library(dplyr)
library(limma)
library(gplots)
library(ggfortify)
library(ggrepel)
library(reshape2)
library(rlang)
library(ggpubr)




make_ECM_volcano <- function(limma_results, P_val_cutoff = 0.05, fold_change_cutoff = 1.5, comparison, title, increased_in, decreased_in, xpos) {
  limma_results <- limma_results %>%
    mutate(Category = case_when(
      Matrisome_Category == "ECM Glycoproteins" ~ "ECM Glycoproteins",
      Matrisome_Category == "Collagens" ~ "Collagens",
      Matrisome_Category == "Proteoglycans" ~ "Proteoglycans",
      Matrisome_Category == "ECM-affiliated Proteins" ~ "ECM-affiliated Proteins",
      Matrisome_Category == "ECM Regulators" ~ "ECM Regulators",
      Matrisome_Category == "Secreted Factors" ~ "Secreted Factors",
      TRUE ~ "Non-ECM Proteins"
    ))
  adj_p_col <- sym(paste0(comparison, " adj p-value"))
  log2fc_col <- sym(paste0(comparison, " log2FC"))
  limma_results %>%
    ggplot(aes(y = -log10(!!adj_p_col), x = !!log2fc_col, colour = Category)) +
    geom_hline(yintercept=-log10(P_val_cutoff), col="red2", linetype = 2, size = .7) +
    geom_vline(xintercept = -log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    geom_vline(xintercept = log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    
    # small grey points for non-ecm proteins
    geom_point(data = limma_results %>% filter(Category == "Non-ECM Proteins"),
               aes(colour = Category), size = 1, alpha = 0.3) +
    
    # Larger points for BM proteins
    geom_point(data = limma_results %>% filter(Category != "Non-ECM Proteins"),
               aes(colour = Category),
               size = 2) +
    
    # black outline for BM and ECM proteins
    geom_point(data = limma_results %>%
                 filter(Category != "Non-ECM Proteins"),
               aes(colour = Category),
               size = 2, shape = 1, stroke = .5, colour = "black") +
    
    # colours for each category, and colour key
    scale_color_manual(values = c("ECM Glycoproteins" = "cyan",
                                  "Collagens" = "maroon1",
                                  "Proteoglycans" = "blue",
                                  "ECM-affiliated Proteins" = "green",
                                  "ECM Regulators" = "orangered2",
                                  "Secreted Factors" = "darkviolet",
                                  "Non-ECM Proteins" = "grey50"),
                       name = "Matrisome Category",
                       breaks = c("Collagens", "ECM Glycoproteins", "Proteoglycans",
                                  "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors",
                                  "Non-ECM Proteins")) + 
    scale_fill_manual(values = c("ECM Glycoproteins" = "cyan",
                                 "Collagens" = "maroon1",
                                 "Proteoglycans" = "blue",
                                 "ECM-affiliated Proteins" = "green",
                                 "ECM Regulators" = "orangered2",
                                 "Secreted Factors" = "darkviolet",
                                 "Non-ECM Proteins" = "grey50"),
                      guide = "none") + 
    theme_bw() +
    labs(title = title,
         x = "Log2( Fold Change )",
         y = "-Log10( P-adj )") +
    annotate("text", x = xpos, y = -.02, label = increased_in, size = 4, colour = "black") +
    annotate("text", x = -xpos, y = -.02, label = decreased_in, size = 4, colour = "black") +
    geom_label_repel(data = limma_results %>%
                       filter(Category != "Non-ECM Proteins"
                              & !!adj_p_col < P_val_cutoff
                              & !!log2fc_col > log2(fold_change_cutoff)
                              | Category != "Non-ECM Proteins"
                              & !!adj_p_col < P_val_cutoff
                              & !!log2fc_col < -log2(fold_change_cutoff)),
                     aes(label = Gene, fill = Category), colour = "black", size = 2.5, alpha = .8,
                     max.overlaps = 20, nudge_y = 0.2, nudge_x = -0.1,
                     show.legend = FALSE)
}




make_ECM_volcano_extralabs <- function(limma_results, P_val_cutoff = 0.05, fold_change_cutoff = 1.5, comparison, title, increased_in, decreased_in, xpos) {
  limma_results <- limma_results %>%
    mutate(Category = case_when(
      Matrisome_Category == "ECM Glycoproteins" ~ "ECM Glycoproteins",
      Matrisome_Category == "Collagens" ~ "Collagens",
      Matrisome_Category == "Proteoglycans" ~ "Proteoglycans",
      Matrisome_Category == "ECM-affiliated Proteins" ~ "ECM-affiliated Proteins",
      Matrisome_Category == "ECM Regulators" ~ "ECM Regulators",
      Matrisome_Category == "Secreted Factors" ~ "Secreted Factors",
      TRUE ~ "Non-ECM Proteins"
    ))
  adj_p_col <- sym(paste0(comparison, " adj p-value"))
  log2fc_col <- sym(paste0(comparison, " log2FC"))
  limma_results %>%
    ggplot(aes(y = -log10(!!adj_p_col), x = !!log2fc_col, colour = Category)) +
    geom_hline(yintercept=-log10(P_val_cutoff), col="red2", linetype = 2, size = .7) +
    geom_vline(xintercept = -log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    geom_vline(xintercept = log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    
    # small grey points for non-ecm proteins
    geom_point(data = limma_results %>% filter(Category == "Non-ECM Proteins"),
               aes(colour = Category), size = 1, alpha = 0.3) +
    
    # Larger points for BM proteins
    geom_point(data = limma_results %>% filter(Category != "Non-ECM Proteins"),
               aes(colour = Category),
               size = 2) +
    
    # black outline for BM and ECM proteins
    geom_point(data = limma_results %>%
                 filter(Category != "Non-ECM Proteins"),
               aes(colour = Category),
               size = 2, shape = 1, stroke = .5, colour = "black") +
    
    # colours for each category, and colour key
    scale_color_manual(values = c("ECM Glycoproteins" = "cyan",
                                  "Collagens" = "maroon1",
                                  "Proteoglycans" = "blue",
                                  "ECM-affiliated Proteins" = "green",
                                  "ECM Regulators" = "orangered2",
                                  "Secreted Factors" = "darkviolet",
                                  "Non-ECM Proteins" = "grey50"),
                       name = "Matrisome Category",
                       breaks = c("Collagens", "ECM Glycoproteins", "Proteoglycans",
                                  "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors",
                                  "Non-ECM Proteins")) + 
    scale_fill_manual(values = c("ECM Glycoproteins" = "cyan",
                                 "Collagens" = "maroon1",
                                 "Proteoglycans" = "blue",
                                 "ECM-affiliated Proteins" = "green",
                                 "ECM Regulators" = "orangered2",
                                 "Secreted Factors" = "darkviolet",
                                 "Non-ECM Proteins" = "grey50"),
                      guide = "none") + 
    theme_bw() +
    labs(title = title,
         x = "Log2( Fold Change )",
         y = "-Log10( P-adj )") +
    annotate("text", x = xpos, y = -.02, label = increased_in, size = 4, colour = "black") +
    annotate("text", x = -xpos, y = -.02, label = decreased_in, size = 4, colour = "black") +
    geom_label_repel(data = limma_results %>%
                       filter(!!adj_p_col < P_val_cutoff
                              & !!log2fc_col > log2(fold_change_cutoff)
                              | !!adj_p_col < P_val_cutoff
                              & !!log2fc_col < -log2(fold_change_cutoff)),
                     aes(label = Gene, fill = Category), colour = "black", size = 2.5, alpha = .8,
                     max.overlaps = 20, nudge_y = 0.2, nudge_x = -0.1,
                     show.legend = FALSE) +
    geom_label_repel(data = limma_results %>%
                       filter(Category != "Non-ECM Proteins"
                              & !!adj_p_col > P_val_cutoff
                              & !!log2fc_col > log2(fold_change_cutoff)
                              | Category != "Non-ECM Proteins"
                              & !!adj_p_col > P_val_cutoff
                              & !!log2fc_col < -log2(fold_change_cutoff)),
                     aes(label = Gene, fill = Category), colour = "black", size = 2.5, alpha = .8,
                     max.overlaps = 20, nudge_y = 0.2, nudge_x = -0.1,
                     show.legend = FALSE)
}




make_ECM_volcano_outline_core <- function(limma_results, P_val_cutoff = 0.05, fold_change_cutoff = 1.5, comparison, title, increased_in, decreased_in, xpos) {
  limma_results <- limma_results %>%
    mutate(Category = case_when(
      Matrisome_Category == "ECM Glycoproteins" ~ "ECM Glycoproteins",
      Matrisome_Category == "Collagens" ~ "Collagens",
      Matrisome_Category == "Proteoglycans" ~ "Proteoglycans",
      Matrisome_Category == "ECM-affiliated Proteins" ~ "ECM-affiliated Proteins",
      Matrisome_Category == "ECM Regulators" ~ "ECM Regulators",
      Matrisome_Category == "Secreted Factors" ~ "Secreted Factors",
      TRUE ~ "Non-ECM Proteins"
    ))
  
  # Adding a new division category for 'core matrisome'
  limma_results <- limma_results %>%
    mutate(Division = ifelse(Matrisome_Division == "Core matrisome", "Core matrisome", "Other"))
  
  adj_p_col <- sym(paste0(comparison, " adj p-value"))
  log2fc_col <- sym(paste0(comparison, " log2FC"))
  
  limma_results %>%
    ggplot(aes(y = -log10(!!adj_p_col), x = !!log2fc_col, colour = Category)) +
    geom_hline(yintercept=-log10(P_val_cutoff), col="red2", linetype = 2, size = .7) +
    geom_vline(xintercept = -log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    geom_vline(xintercept = log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    
    # small grey points for non-ecm proteins
    geom_point(data = limma_results %>% filter(Category == "Non-ECM Proteins"),
               aes(colour = Category), size = 1, alpha = 0.3) +
    
    # Larger points for BM proteins
    geom_point(data = limma_results %>% filter(Category != "Non-ECM Proteins"),
               aes(colour = Category),
               size = 2) +

    # black outline for 'core matrisome' proteins
    geom_point(data = limma_results %>%
                 filter(Division == "Core matrisome"),
               aes(colour = Category),
               size = 2, shape = 1, stroke = 1.3, colour = "black") +
    
    # colours for each category, and colour key
    scale_color_manual(values = c("ECM Glycoproteins" = "cyan",
                                  "Collagens" = "maroon1",
                                  "Proteoglycans" = "blue",
                                  "ECM-affiliated Proteins" = "green",
                                  "ECM Regulators" = "orangered2",
                                  "Secreted Factors" = "darkviolet",
                                  "Non-ECM Proteins" = "grey50"),
                       name = "Matrisome Category",
                       breaks = c("Collagens", "ECM Glycoproteins", "Proteoglycans",
                                  "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors",
                                  "Non-ECM Proteins")) + 

    scale_fill_manual(values = c("ECM Glycoproteins" = "cyan",
                                 "Collagens" = "maroon1",
                                 "Proteoglycans" = "blue",
                                 "ECM-affiliated Proteins" = "green",
                                 "ECM Regulators" = "orangered2",
                                 "Secreted Factors" = "darkviolet",
                                 "Non-ECM Proteins" = "grey50"),
                      guide = "none") + 
    theme_bw() +
    labs(title = title,
         x = "Log2( Fold Change )",
         y = "-Log10( P-adj )") +
    annotate("text", x = xpos, y = -.02, label = increased_in, size = 4, colour = "black") +
    annotate("text", x = -xpos, y = -.02, label = decreased_in, size = 4, colour = "black") +
    geom_label_repel(data = limma_results %>%
                       filter(Category != "Non-ECM Proteins"
                              & !!adj_p_col < P_val_cutoff
                              & !!log2fc_col > log2(fold_change_cutoff)
                              | Category != "Non-ECM Proteins"
                              & !!adj_p_col < P_val_cutoff
                              & !!log2fc_col < -log2(fold_change_cutoff)),
                     aes(label = Gene, fill = Category), colour = "black", size = 2.5, alpha = .8,
                     max.overlaps = 20, nudge_y = 0.2, nudge_x = -0.1,
                     show.legend = FALSE)
}



###




# import full limma results data output from analysis script (should have log FC, pvals etc, plus protein ECM categories)
xxx <- read.csv("aim3_allprots_protein_dataset_DE_results.csv", check.names = FALSE)


# plots

# just ECM categories coloured (can change colours in function above ^ ) ---
make_ECM_volcano(sample_DE, comparison = "Condition_A-Condition_B", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                 title = "Differential Expression of Proteins in Condition A and B",
                 increased_in = "Increased in condition B",
                 decreased_in = "Increased in condition A",
                 xpos = 2)

# all prots labeled p<0.05 and ECM prots labelled FC>1.5 (non-sig too)
make_ECM_volcano_extralabs(sample_DE, comparison = "Condition_A-Condition_B", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                           title = "Differential Expression of Proteins in Condition A and B",
                           increased_in = "Increased in condition B",
                           decreased_in = "Increased in condition A",
                           xpos = 2)

# ECM categories coloured, and 'core matrisome' proteins have a bold black outline ---
make_ECM_volcano_outline_core(sample_DE, comparison = "Condition_A-Condition_B", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                              title = "Differential Expression of Proteins in Condition A and B",
                              increased_in = "Increased in condition B",
                              decreased_in = "Increased in condition A",
                              xpos = 2)





# make a PDF

pdf("sample_ECM_volcanoes.pdf", width = 10, height = 6)

dev.off()









