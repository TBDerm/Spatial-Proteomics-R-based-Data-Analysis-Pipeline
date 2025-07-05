

library(tidyverse)
library(data.table)
library(dplyr)
library(gplots)
library(ggfortify)
library(ggrepel)
library(reshape2)
library(rlang)
library(ggpubr)



a <- both_results %>%
  ggplot(aes(x = `Buttock-Forearm log2FC`, y = -log10(`Buttock-Forearm adj p-value`))) +
  geom_hline(yintercept = -log10(0.05), col = "red2", linetype = 2, size = 0.7) + 
  geom_vline(xintercept = -log2(1.5), col = "red4", linetype = 3, size = 0.7) +
  geom_vline(xintercept = log2(1.5), col = "red4", linetype = 3, size = 0.7) +
  geom_point(size = 1, colour = "grey", alpha = 0.5) +
  geom_point(data = both_results %>% 
               filter(`Buttock-Forearm adj p-value` < 0.05 & abs(`Buttock-Forearm log2FC`) > log2(1.5)),
             size = 1.5, colour = "red", alpha = 0.5) +
  # xlim(-15,15) +
  theme_bw() +
  labs(title = paste("Volcano Plot: R script results"),
       x = "Log2( Fold Change )",
       y = "-Log10( P-adj )") +
  geom_label_repel(data = both_results %>%  
                     filter(`Buttock-Forearm adj p-value` < 0.05 & abs(`Buttock-Forearm log2FC`) > log2(1.5)),
                   aes(label = Gene), colour = "black", size = 3,
                   max.overlaps = 20,
                   nudge_y = 0.2,
                   nudge_x = -0.1)





b <- both_results %>%
  ggplot(aes(x = `foldchange.log2_deqms_contrast: buttock vs forearm # condition_variable: group`, 
             y = -log10(`pvalue_deqms_contrast: buttock vs forearm # condition_variable: group`))) +
  geom_hline(yintercept = -log10(0.05), col = "red2", linetype = 2, size = 0.7) + 
  geom_vline(xintercept = -log2(1.5), col = "red4", linetype = 3, size = 0.7) +
  geom_vline(xintercept = log2(1.5), col = "red4", linetype = 3, size = 0.7) +
  geom_point(size = 1, colour = "grey", alpha = 0.5) +
  geom_point(data = both_results %>%
               filter(`pvalue_deqms_contrast: buttock vs forearm # condition_variable: group` < 0.05 
                      & abs(`foldchange.log2_deqms_contrast: buttock vs forearm # condition_variable: group`) > log2(1.5)),
             size = 1.5, colour = "red", alpha = 0.5) +
  # xlim(-5,5) +
  theme_bw() +
  labs(title = paste("Volcano Plot: MS-DAP deqms results"),
       x = "Log2( Fold Change )",
       y = "-Log10( P-adj )") +
  geom_label_repel(data = both_results %>%  
                     filter(`pvalue_deqms_contrast: buttock vs forearm # condition_variable: group` < 0.05 
                            & abs(`foldchange.log2_deqms_contrast: buttock vs forearm # condition_variable: group`) > log2(1.5)),
                   aes(label = Gene), colour = "black", size = 3,
                   max.overlaps = 20,
                   nudge_y = 0.2,
                   nudge_x = -0.1)
    


c <- both_results %>%
  ggplot(aes(x = `foldchange.log2_msempire_contrast: buttock vs forearm # condition_variable: group`, 
             y = -log10(`pvalue_msempire_contrast: buttock vs forearm # condition_variable: group`))) +
  geom_hline(yintercept = -log10(0.05), col = "red2", linetype = 2, size = 0.7) + 
  geom_vline(xintercept = -log2(1.5), col = "red4", linetype = 3, size = 0.7) +
  geom_vline(xintercept = log2(1.5), col = "red4", linetype = 3, size = 0.7) +
  geom_point(size = 1, colour = "grey", alpha = 0.5) +
  geom_point(data = both_results %>%
               filter(`pvalue_msempire_contrast: buttock vs forearm # condition_variable: group` < 0.05 
                      & abs(`foldchange.log2_msempire_contrast: buttock vs forearm # condition_variable: group`) > log2(1.5)),
             size = 1.5, colour = "red", alpha = 0.5) +
  # xlim(-5,5) +
  theme_bw() +
  labs(title = paste("Volcano Plot: MS-DAP msempire results"),
       x = "Log2( Fold Change )",
       y = "-Log10( P-adj )") +
  geom_label_repel(data = both_results %>%  
                     filter(`pvalue_msempire_contrast: buttock vs forearm # condition_variable: group` < 0.05 
                            & abs(`foldchange.log2_msempire_contrast: buttock vs forearm # condition_variable: group`) > log2(1.5)),
                   aes(label = Gene), colour = "black", size = 3,
                   max.overlaps = 20,
                   nudge_y = 0.2,
                   nudge_x = -0.1)



d <- both_results %>%
  ggplot(aes(x = `foldchange.log2_msqrob_contrast: buttock vs forearm # condition_variable: group`, 
             y = -log10(`pvalue_msqrob_contrast: buttock vs forearm # condition_variable: group`))) +
  geom_hline(yintercept = -log10(0.05), col = "red2", linetype = 2, size = 0.7) + 
  geom_vline(xintercept = -log2(1.5), col = "red4", linetype = 3, size = 0.7) +
  geom_vline(xintercept = log2(1.5), col = "red4", linetype = 3, size = 0.7) +
  geom_point(size = 1, colour = "grey", alpha = 0.5) +
  geom_point(data = both_results %>%
               filter(`pvalue_msqrob_contrast: buttock vs forearm # condition_variable: group` < 0.05 
                      & abs(`foldchange.log2_msqrob_contrast: buttock vs forearm # condition_variable: group`) > log2(1.5)),
             size = 1.5, colour = "red", alpha = 0.5) +
  # xlim(-5,5) +
  theme_bw() +
  labs(title = paste("Volcano Plot: MS-DAP msqrob results"),
       x = "Log2( Fold Change )",
       y = "-Log10( P-adj )") +
  geom_label_repel(data = both_results %>%  
                     filter(`pvalue_msqrob_contrast: buttock vs forearm # condition_variable: group` < 0.05 
                            & abs(`foldchange.log2_msqrob_contrast: buttock vs forearm # condition_variable: group`) > log2(1.5)),
                   aes(label = Gene), colour = "black", size = 3,
                   max.overlaps = 20,
                   nudge_y = 0.2,
                   nudge_x = -0.1)




#   
#   
# pdf("MS-DAP_comp_volcano_plots.pdf", width = 17, height = 13)
# 
# ggarrange(a, 
#           b, 
#           c,
#           d,
#           ncol = 2, nrow = 2)
# 
# dev.off()
#   
#   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  