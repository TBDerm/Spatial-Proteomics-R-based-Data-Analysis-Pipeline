
library(tidyverse)
library(data.table)
library(dplyr)
library(gplots)
library(ggfortify)
library(ggrepel)
library(reshape2)
library(rlang)
library(ggpubr)



# make_volcano_BM <- function(limma_results, P_val_cutoff = 0.05, fold_change_cutoff = 1.5, comparison, title) {
#   limma_results <- limma_results %>%
#     mutate(Category = case_when(
#       BM_Status == "Core" ~ "Core BM Protein", 
#       BM_Status == "Confirmed" ~ "Confirmed BM Protein",
#       TRUE ~ "Other"
#     ))
#   adj_p_col <- sym(paste0(comparison, " adj p-value"))
#   log2fc_col <- sym(paste0(comparison, " log2FC"))
#   limma_results %>%
#     ggplot(aes(y = -log10(!!adj_p_col), x = !!log2fc_col, colour = Category)) +
#     geom_hline(yintercept=-log10(P_val_cutoff), col="red2", linetype = 2, size = .7) +
#     geom_vline(xintercept = -log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
#     geom_vline(xintercept = log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
#     
#     # small grey points for "Other" proteins
#     geom_point(data = limma_results %>% filter(Category == "Other"),
#                aes(colour = Category), size = 1, alpha = 0.3) +
#    
#     # Larger points for BM proteins
#     geom_point(data = limma_results %>% filter(Category != "Other"),
#                aes(colour = Category),
#                size = 2) +
#     
#     # black outline for BM and ECM proteins
#     geom_point(data = limma_results %>%
#                  filter(Category != "Other"),
#                aes(colour = Category),
#                size = 2, shape = 1, stroke = .5, colour = "black") +
#     
#     # colours for each category, and colour key
#     scale_colour_manual(values = c(
#       "Core BM Protein" = "chartreuse",
#       "Confirmed BM Protein" = "cyan",
#       "Other" = "grey50"
#     ), name = "Protein Category") + 
#     
#     theme_bw() +
#     labs(title = title,
#          x = "Log2( Fold Change )",
#          y = "-Log10( P-adj )") +
#     geom_label_repel(data = limma_results %>%
#                        filter(Category != "Other" 
#                               & !!adj_p_col < P_val_cutoff 
#                               & !!log2fc_col > log2(fold_change_cutoff)
#                               | Category != "Other"
#                               & !!adj_p_col < P_val_cutoff
#                               & !!log2fc_col < -log2(fold_change_cutoff)),
#                      aes(label = Gene), colour = "black", size = 2.5,
#                      max.overlaps = 20,
#                      nudge_y = 0.2,
#                      nudge_x = -0.1)
# }



make_volcano_ECM_BM_noLabs <- function(limma_results, P_val_cutoff = 0.05, fold_change_cutoff = 1.5, comparison, title, increased_in, decreased_in, xpos) {
  limma_results <- limma_results %>%
    mutate(Category = case_when(
      BM_Status == "Core" ~ "Core BM Protein",
      BM_Status == "Confirmed" ~ "Confirmed BM Protein",
      !is.na(Matrisome_Category) ~ "ECM Protein",
      TRUE ~ "Other"
    ))
  adj_p_col <- sym(paste0(comparison, " adj p-value"))
  log2fc_col <- sym(paste0(comparison, " log2FC"))
  limma_results %>%
    ggplot(aes(y = -log10(!!adj_p_col), x = !!log2fc_col, colour = Category)) +
    geom_hline(yintercept=-log10(P_val_cutoff), col="red2", linetype = 2, size = .7) +
    geom_vline(xintercept = -log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    geom_vline(xintercept = log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    
    # small grey points for "Other" proteins
    geom_point(data = limma_results %>% filter(Category == "Other"),
               aes(colour = Category), size = 1, alpha = 0.3) +
    
    # Larger points for BM and ECM proteins
    geom_point(data = limma_results %>% filter(Category != "Other"),
               aes(colour = Category),
               size = 2) +
    
    # black outline for BM and ECM proteins
    geom_point(data = limma_results %>%
                 filter(Category != "Other"),
               aes(colour = Category),
               size = 2, shape = 1, stroke = .5, colour = "black") +
    
    # colours for each category, and colour key
    scale_colour_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), name = "Protein Category",
    breaks = c("Core BM Protein",
               "Confirmed BM Protein",
               "ECM Protein",
               "Other")) + 
    scale_fill_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), guide = "none") + 
    theme_bw() +
    labs(title = title,
         x = "Log2( Fold Change )",
         y = "-Log10( P-adj )") +
    annotate("text", x = xpos, y = -.02, label = increased_in, size = 4, colour = "black") +
    annotate("text", x = -xpos, y = -.02, label = decreased_in, size = 4, colour = "black")
}


make_volcano_ECM_BM <- function(limma_results, P_val_cutoff = 0.05, fold_change_cutoff = 1.5, comparison, title, increased_in, decreased_in, xpos) {
  limma_results <- limma_results %>%
    mutate(Category = case_when(
      BM_Status == "Core" ~ "Core BM Protein",
      BM_Status == "Confirmed" ~ "Confirmed BM Protein",
      !is.na(Matrisome_Category) ~ "ECM Protein",
      TRUE ~ "Other"
    ))
  adj_p_col <- sym(paste0(comparison, " adj p-value"))
  log2fc_col <- sym(paste0(comparison, " log2FC"))
  limma_results %>%
    ggplot(aes(y = -log10(!!adj_p_col), x = !!log2fc_col, colour = Category)) +
    geom_hline(yintercept=-log10(P_val_cutoff), col="red2", linetype = 2, size = .7) +
    geom_vline(xintercept = -log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    geom_vline(xintercept = log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    
    # small grey points for "Other" proteins
    geom_point(data = limma_results %>% filter(Category == "Other"),
               aes(colour = Category), size = 1, alpha = 0.3) +
    
    # Larger points for BM and ECM proteins
    geom_point(data = limma_results %>% filter(Category != "Other"),
               aes(colour = Category),
               size = 2) +
    
    # black outline for BM and ECM proteins
    geom_point(data = limma_results %>%
                 filter(Category != "Other"),
               aes(colour = Category),
               size = 2, shape = 1, stroke = .5, colour = "black") +
    
    # colours for each category, and colour key
    scale_colour_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), name = "Protein Category",
    breaks = c("Core BM Protein",
               "Confirmed BM Protein",
               "ECM Protein",
               "Other")) + 
    scale_fill_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), guide = "none") + 
    theme_bw() +
    labs(title = title,
         x = "Log2( Fold Change )",
         y = "-Log10( P-adj )") +
    annotate("text", x = xpos, y = -.02, label = increased_in, size = 4, colour = "black") +
    annotate("text", x = -xpos, y = -.02, label = decreased_in, size = 4, colour = "black") +
    geom_label_repel(data = limma_results %>%
                       filter(Category != "Other"
                              & !!adj_p_col < P_val_cutoff 
                              & !!log2fc_col > log2(fold_change_cutoff)
                              | Category != "Other"
                              & !!adj_p_col < P_val_cutoff
                              & !!log2fc_col < -log2(fold_change_cutoff)),
                     aes(label = Gene, fill = Category), colour = "black", size = 2.5, alpha = .8,
                     max.overlaps = 20, nudge_y = 0.2, nudge_x = -0.1,
                     show.legend = FALSE)
}


make_volcano_ECM_BM_allLabs <- function(limma_results, P_val_cutoff = 0.05, fold_change_cutoff = 1.5, comparison, title, increased_in, decreased_in, xpos) {
  limma_results <- limma_results %>%
    mutate(Category = case_when(
      BM_Status == "Core" ~ "Core BM Protein",
      BM_Status == "Confirmed" ~ "Confirmed BM Protein",
      !is.na(Matrisome_Category) ~ "ECM Protein",
      TRUE ~ "Other"
    ))
  adj_p_col <- sym(paste0(comparison, " adj p-value"))
  log2fc_col <- sym(paste0(comparison, " log2FC"))
  limma_results %>%
    ggplot(aes(y = -log10(!!adj_p_col), x = !!log2fc_col, colour = Category)) +
    geom_hline(yintercept=-log10(P_val_cutoff), col="red2", linetype = 2, size = .7) +
    geom_vline(xintercept = -log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    geom_vline(xintercept = log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    
    # small grey points for "Other" proteins
    geom_point(data = limma_results %>% filter(Category == "Other"),
               aes(colour = Category), size = 1, alpha = 0.3) +
    
    # Larger points for BM and ECM proteins
    geom_point(data = limma_results %>% filter(Category != "Other"),
               aes(colour = Category),
               size = 2) +
    
    # black outline for BM and ECM proteins
    geom_point(data = limma_results %>%
                 filter(Category != "Other"),
               aes(colour = Category),
               size = 2, shape = 1, stroke = .5, colour = "black") +
    
    # colours for each category, and colour key
    scale_colour_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), name = "Protein Category",
    breaks = c("Core BM Protein",
               "Confirmed BM Protein",
               "ECM Protein",
               "Other")) + 
    scale_fill_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), guide = "none") + 
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
                     show.legend = FALSE)
}


# all labels p<0.05 and ECM/BM labelled all FC>1.5 (not sig)
make_volcano_ECM_BM_allLabs_FC <- function(limma_results, P_val_cutoff = 0.05, fold_change_cutoff = 1.5, comparison, title, increased_in, decreased_in, xpos) {
  limma_results <- limma_results %>%
    mutate(Category = case_when(
      BM_Status == "Core" ~ "Core BM Protein",
      BM_Status == "Confirmed" ~ "Confirmed BM Protein",
      !is.na(Matrisome_Category) ~ "ECM Protein",
      TRUE ~ "Other"
    ))
  adj_p_col <- sym(paste0(comparison, " adj p-value"))
  log2fc_col <- sym(paste0(comparison, " log2FC"))
  limma_results %>%
    ggplot(aes(y = -log10(!!adj_p_col), x = !!log2fc_col, colour = Category)) +
    geom_hline(yintercept=-log10(P_val_cutoff), col="red2", linetype = 2, size = .7) +
    geom_vline(xintercept = -log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    geom_vline(xintercept = log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    
    # small grey points for "Other" proteins
    geom_point(data = limma_results %>% filter(Category == "Other"),
               aes(colour = Category), size = 1, alpha = 0.3) +
    
    # Larger points for BM and ECM proteins
    geom_point(data = limma_results %>% filter(Category != "Other"),
               aes(colour = Category),
               size = 2) +
    
    # black outline for BM and ECM proteins
    geom_point(data = limma_results %>%
                 filter(Category != "Other"),
               aes(colour = Category),
               size = 2, shape = 1, stroke = .5, colour = "black") +
    
    # colours for each category, and colour key
    scale_colour_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), name = "Protein Category",
    breaks = c("Core BM Protein",
               "Confirmed BM Protein",
               "ECM Protein",
               "Other")) + 
    scale_fill_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), guide = "none") + 
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
                       filter(Category != "Other"
                              & !!adj_p_col > P_val_cutoff 
                              & !!log2fc_col > log2(fold_change_cutoff)
                              | Category != "Other"
                              & !!adj_p_col > P_val_cutoff
                              & !!log2fc_col < -log2(fold_change_cutoff)),
                     aes(label = Gene, fill = Category), colour = "black", size = 2.5, alpha = .8,
                     max.overlaps = 20, nudge_y = 0.2, nudge_x = -0.1,
                     show.legend = FALSE)
}



# all labels p<0.05 and just BM labelled all FC>1.5 (not sig)
make_volcano_ECM_BM_allLabs_FC_BM <- function(limma_results, P_val_cutoff = 0.05, fold_change_cutoff = 1.5, comparison, title, increased_in, decreased_in, xpos) {
  limma_results <- limma_results %>%
    mutate(Category = case_when(
      BM_Status == "Core" ~ "Core BM Protein",
      BM_Status == "Confirmed" ~ "Confirmed BM Protein",
      !is.na(Matrisome_Category) ~ "ECM Protein",
      TRUE ~ "Other"
    ))
  adj_p_col <- sym(paste0(comparison, " adj p-value"))
  log2fc_col <- sym(paste0(comparison, " log2FC"))
  limma_results %>%
    ggplot(aes(y = -log10(!!adj_p_col), x = !!log2fc_col, colour = Category)) +
    geom_hline(yintercept=-log10(P_val_cutoff), col="red2", linetype = 2, size = .7) +
    geom_vline(xintercept = -log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    geom_vline(xintercept = log2(fold_change_cutoff), col="red4", linetype = 3, size = .7) +
    
    # small grey points for "Other" proteins
    geom_point(data = limma_results %>% filter(Category == "Other"),
               aes(colour = Category), size = 1, alpha = 0.3) +
    
    # Larger points for BM and ECM proteins
    geom_point(data = limma_results %>% filter(Category != "Other"),
               aes(colour = Category),
               size = 2) +
    
    # black outline for BM and ECM proteins
    geom_point(data = limma_results %>%
                 filter(Category != "Other"),
               aes(colour = Category),
               size = 2, shape = 1, stroke = .5, colour = "black") +
    
    # colours for each category, and colour key
    scale_colour_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), name = "Protein Category",
    breaks = c("Core BM Protein",
               "Confirmed BM Protein",
               "ECM Protein",
               "Other")) + 
    scale_fill_manual(values = c(
      "Core BM Protein" = "chartreuse3",
      "Confirmed BM Protein" = "dodgerblue",
      "ECM Protein" = "gold",
      "Other" = "grey50"
    ), guide = "none") + 
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
                       filter(!is.na(BM_Status)
                              & !!adj_p_col > P_val_cutoff 
                              & !!log2fc_col > log2(fold_change_cutoff)
                              | !is.na(BM_Status)
                              & !!adj_p_col > P_val_cutoff
                              & !!log2fc_col < -log2(fold_change_cutoff)
                              ),
                     aes(label = Gene, fill = Category), colour = "black", size = 2.5, alpha = .8,
                     max.overlaps = 20, nudge_y = 0.2, nudge_x = -0.1,
                     show.legend = FALSE)
}





########################




# import full limma results data (should have log FC, pvals etc, plus protein ECM/BM categories)


aim2_DE <- read.csv("sample_protein_dataset_DE_results.csv", check.names = FALSE)

aim3_DE <- read.csv("aim3_allprots_protein_dataset_DE_results.csv", check.names = FALSE)

sample_DE <- read.csv("sample_protein_dataset_DE_results.csv", check.names = FALSE)


# G_frag <- read.csv("aim3_N3_DIANN_noMods_allProts_protein_dataset_valid_imputed_limma_results.csv", check.names = FALSE)
# T_frag <- read.csv("aim3_N3_DIANN_wMods_Tess_Tuned_allProts_protein_dataset_valid_imputed_limma_results.csv", check.names = FALSE)
# G_tuned_mods_DIANN <- read.csv("aim3_N3_DIANN_noMods_allProts_protein_dataset_valid_imputed_limma_results.csv", check.names = FALSE)
# T_tuned_mods_DIANN <- read.csv("aim3_N3_DIANN_wMods_Tess_Tuned_allProts_protein_dataset_valid_imputed_limma_results.csv", check.names = FALSE)
# no_tuning_mods_DIANN <- read.csv("aim3_N3_DIANN_wMods_noTuning_allProts_protein_dataset_valid_imputed_limma_results.csv", check.names = FALSE)
# no_mods_DIANN <- read.csv("aim3_N3_DIANN_noMods_allProts_protein_dataset_valid_imputed_limma_results.csv", check.names = FALSE)





# plots -----------------------------------------------------------------------------------

make_volcano_ECM_BM(aim2_DE, comparison = "Epidermis_Intact-Epidermis_Removed", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                    title = "Differential Expression of Proteins with Epidermis Removal",
                    increased_in = "Increased with epidermis removal",
                    decreased_in = "Decreased with epidermis removal",
                    xpos = 2)

make_volcano_ECM_BM_noLabs(full_imputed_dataset_limma_results, comparison = "Epidermis_Intact-Epidermis_Removed", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                           title = "Differential Expression of Proteins with Epidermis Removal",
                           increased_in = "Increased with epidermis removal",
                           decreased_in = "Decreased with epidermis removal",
                           xpos = 2)

make_volcano_ECM_BM_allLabs(full_imputed_dataset_limma_results, comparison = "Epidermis_Intact-Epidermis_Removed", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                           title = "Differential Expression of Proteins with Epidermis Removal",
                           increased_in = "Increased with epidermis removal",
                           decreased_in = "Decreased with epidermis removal",
                           xpos = 2)

###


make_volcano_ECM_BM_noLabs(sample_DE, comparison = "Condition_A-Condition_B", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                    title = "Differential Expression of Proteins in Condition A and B",
                    increased_in = "Increased in condition B",
                    decreased_in = "Increased in condition A",
                    xpos = 2)

make_volcano_ECM_BM(sample_DE, comparison = "Condition_A-Condition_B", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                    title = "Differential Expression of Proteins in Condition A and B",
                    increased_in = "Increased in condition B",
                    decreased_in = "Increased in condition A",
                    xpos = 2)

make_volcano_ECM_BM_allLabs(sample_DE, comparison = "Condition_A-Condition_B", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                            title = "Differential Expression of Proteins in Condition A and B",
                            increased_in = "Increased in condition B",
                            decreased_in = "Increased in condition A",
                            xpos = 2)

make_volcano_ECM_BM_allLabs_FC(sample_DE, comparison = "Condition_A-Condition_B", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                               title = "Differential Expression of Proteins in Condition A and B",
                               increased_in = "Increased in condition B",
                               decreased_in = "Increased in condition A",
                               xpos = 2)

make_volcano_ECM_BM_allLabs_FC_BM(aim3_DE, comparison = "Buttock-Forearm", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                                  title = "Differential Expression of Proteins in Photoaged DEJ",
                                  increased_in = "Increased in photoaged",
                                  decreased_in = "Decreased in photoaged",
                                  xpos = 1)






pdf("aim3_touse_volcanoes.pdf", width = 10, height = 6)

dev.off()


# 
# pdf("big_volcano.pdf", width = 9, height = 5)
# 
# dev.off()



pdf("aim3_0.05_volcanoes.pdf", width = 10, height = 6)

dev.off()











pdf("all_speclibs_volcano_plots_2.pdf", width = 17, height = 5)

ggarrange(make_volcano_ECM_BM2(G_frag, comparison = "Buttock-Forearm", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                               title = "Gulsev's Fragpipe speclib", 
                               increased_in = "", decreased_in = "", xpos = 1), 
          make_volcano_ECM_BM2(G_frag, comparison = "Buttock-Forearm", P_val_cutoff = 0.1, fold_change_cutoff = 1.5, 
                               title = "Gulsev's Fragpipe speclib", 
                               increased_in = "", decreased_in = "", xpos = 1), 
          ncol = 2, nrow = 1)

ggarrange(make_volcano_ECM_BM(T_frag, comparison = "Buttock-Forearm", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                              title = "Tess' Fragpipe speclib", 
                              increased_in = "", decreased_in = "", xpos = 1), 
          make_volcano_ECM_BM(T_frag, comparison = "Buttock-Forearm", P_val_cutoff = 0.1, fold_change_cutoff = 1.5, 
                              title = "Tess' Fragpipe speclib", 
                              increased_in = "", decreased_in = "", xpos = 1), 
          ncol = 2, nrow = 1)

ggarrange(make_volcano_ECM_BM2(G_tuned_mods_DIANN, comparison = "Buttock-Forearm", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                               title = "DIANN library with mods (Gulsev's tuned)", 
                               increased_in = "", decreased_in = "", xpos = 1),
          make_volcano_ECM_BM2(G_tuned_mods_DIANN, comparison = "Buttock-Forearm", P_val_cutoff = 0.1, fold_change_cutoff = 1.5, 
                               title = "DIANN library with mods (Gulsev's tuned)", 
                               increased_in = "", decreased_in = "", xpos = 1),
          ncol = 2, nrow = 1)

ggarrange(make_volcano_ECM_BM(T_tuned_mods_DIANN, comparison = "Buttock-Forearm", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                              title = "DIANN library with mods (Tess' tuned)", 
                              increased_in = "", decreased_in = "", xpos = 1),
          make_volcano_ECM_BM(T_tuned_mods_DIANN, comparison = "Buttock-Forearm", P_val_cutoff = 0.1, fold_change_cutoff = 1.5, 
                              title = "DIANN library with mods (Tess' tuned)", 
                              increased_in = "", decreased_in = "", xpos = 1),
          ncol = 2, nrow = 1)

ggarrange(make_volcano_ECM_BM2(no_tuning_mods_DIANN, comparison = "Buttock-Forearm", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                               title = "DIANN library with mods (no tuning)", 
                               increased_in = "", decreased_in = "", xpos = 1), 
          make_volcano_ECM_BM2(no_tuning_mods_DIANN, comparison = "Buttock-Forearm", P_val_cutoff = 0.1, fold_change_cutoff = 1.5, 
                               title = "DIANN library with mods (no tuning)", 
                               increased_in = "", decreased_in = "", xpos = 1), 
          ncol = 2, nrow = 1)

ggarrange(make_volcano_ECM_BM2(no_mods_DIANN, comparison = "Buttock-Forearm", P_val_cutoff = 0.05, fold_change_cutoff = 1.5, 
                               title = "DIANN library, no mods", 
                               increased_in = "", decreased_in = "", xpos = 1), 
          make_volcano_ECM_BM2(no_mods_DIANN, comparison = "Buttock-Forearm", P_val_cutoff = 0.1, fold_change_cutoff = 1.5, 
                               title = "DIANN library, no mods", 
                               increased_in = "", decreased_in = "", xpos = 1), 
          ncol = 2, nrow = 1)

dev.off()







