
library(tidyverse)
library(GGally)
library(dplyr)



# makes data into useable format - choose protein, then it gets sample values, puts them in two cols (one per condition), and that's it
data_formater <- function(data, protein) {
  data %>%
    filter(Gene == protein) %>%
    pivot_longer(cols = c(sample_cols_A, sample_cols_B),
                 names_to = "sample",
                 values_to = "value") %>%
    mutate(group = case_when(
      sample %in% sample_cols_A ~ condition_A,
      sample %in% sample_cols_B ~ condition_B,
      TRUE ~ NA_character_
    )) %>%
    group_by(Gene, group) %>%
    arrange(Gene, group, sample) %>%
    mutate(replicate = row_number()) %>%
    ungroup() %>%
    select(Gene, sample, value, group, replicate) %>%
    pivot_wider(
      id_cols = c(Gene, replicate),
      names_from = group,
      values_from = value
    ) %>%
    arrange(Gene, replicate) %>%
    select(condition_A, condition_B)
}



dot_line_plot <- function(data, protein, sig) {
  plot_data <- data_formater(data = data, protein = protein)
  plot_data %>%
    ggparcoord(alphaLines = 0,
               showPoints = TRUE,
               boxplot = TRUE,
               scale = "globalminmax") +
    geom_point(colour = "black", size = 3) +
    geom_line(colour = "grey22", size = 1, linetype = 2, alpha = .5) +
    theme_bw() +
    ylim((min(unlist(plot_data))), 1.015*(max(unlist(plot_data)))) +
    labs(title = protein, size = 28) +
    theme(axis.text.y = element_text(size = 26), 
          axis.text.x = element_text(size = 22),
          axis.title = element_blank(),
          plot.title = element_text(size = 28, hjust = .5)) +
    annotate("text", x = 1.5, y = 1.014*(max(unlist(plot_data))), label = sig, size = 9, colour = "red2") +
    annotate("segment", x = 1, xend = 2, y = 1.004*(max(unlist(plot_data))), yend = 1.004*(max(unlist(plot_data))), colour = "red2", size = 1.6)
}


dot_line_plot_autosig <- function(data, protein) {
  p_val_colname <- paste0(condition_A, "-", condition_B, " adj p-value")
  pval <- data[data$Gene == protein, p_val_colname][1]
  SigLabel = case_when(
    pval < 0.0001 ~ "****",
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ NA_character_
  )
  plot_data <- data_formater(data = data, protein = protein)
  p <- plot_data %>%
    ggparcoord(alphaLines = 0,
               showPoints = TRUE,
               boxplot = TRUE,
               scale = "globalminmax") +
    geom_point(colour = "black", size = 3) +
    geom_line(colour = "grey22", size = 1, linetype = 2, alpha = .5) +
    theme_bw() +
    labs(title = protein, size = 28) +
    theme(axis.text.y = element_text(size = 26), 
          axis.text.x = element_text(size = 25),
          axis.title = element_blank(),
          plot.title = element_text(size = 30, hjust = .5))
  if(!is.na(SigLabel)) {
    p <- p +
      annotate("text", x = 1.5, y = 1.005*(max(unlist(plot_data))), label = SigLabel, size = 15, colour = "red2") +
      annotate("segment", x = 1, xend = 2, y = 1.004*(max(unlist(plot_data))), yend = 1.004*(max(unlist(plot_data))), colour = "red2", size = 1.6) +
      ylim((min(unlist(plot_data))), 1.007*(max(unlist(plot_data))))
  }
  return(p)
}





# read data
aim3_DE <- read.csv("aim3_allprots_protein_dataset_DE_results.csv", check.names = FALSE)

# specify conditions, as they are in data
condition_A <- "Buttock"
condition_B <- "Forearm"

# sample colnames in data
sample_cols_A <- c("0841_Buttock", "0842_Buttock", "0863_Buttock")
sample_cols_B <- c("0841_Forearm", "0842_Forearm", "0863_Forearm")




# plot
dot_line_plot(aim3_DE, "COL14A1", sig = "p-adj < 0.1")

dot_line_plot_autosig(aim3_DE, "APCS")
dot_line_plot_autosig(aim3_DE, "COL6A1")





pdf("aim3_dot_line_plots.pdf", width = 7, height = 10)

dot_line_plot_autosig(aim3_DE, "APCS")
dot_line_plot_autosig(aim3_DE, "VTN")
dot_line_plot(aim3_DE, "COL14A1", sig = "p-adj < 0.1")

dev.off()









#------ DEJ undulation ---------------------

DEJ_scores <- read.csv("DEJ_scores.csv")

DEJ_scores %>%
  ggparcoord(alphaLines = 0,
             showPoints = TRUE,
             boxplot = TRUE,
             scale = "globalminmax") +
  geom_point(colour = "black", size = 3) +
  theme_bw() +
  labs(title = "DEJ undulation", size = 20) +
  theme(axis.text.y = element_text(size = 22), 
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 25, hjust = .5)) +
  annotate("text", x = 1.5, y = 2.27, label = "***", size = 15, colour = "red2") +
  annotate("segment", x = 1, xend = 2, y = 2.25, yend = 2.25, colour = "red2", size = 1.6)







# ------------------- dot-line graphs with mean and SD --------------------------



dot_plot <- function(data, protein, sig, sigsize = 14, ylev=1.014) {
  plot_data <- data_formater(data = data, protein = protein)
  means <- colMeans(plot_data)
  sds <- apply(plot_data, 2, sd)
  plot_data %>%
    ggparcoord(alphaLines = 0,
               showPoints = TRUE,
               boxplot = F,
               scale = "globalminmax") +     
    geom_point(aes(x = 1, y = means[1]), color = "grey50", size = 10) + 
    geom_errorbar(aes(x = 1, ymin = means[1] - sds[1], ymax = means[1] + sds[1]),
                  color = "grey50", width = 0.15, size = 1) + 
    geom_point(aes(x = 2, y = means[2]), color = "grey50", size = 10) + 
    geom_errorbar(aes(x = 2, ymin = means[2] - sds[2], ymax = means[2] + sds[2]),
                  color = "grey50", width = 0.15, size = 1) +
    geom_point(colour = "black", size = 5) +
    geom_line(colour = "grey22", size = 1, linetype = 2, alpha = .5) +
    theme_bw() +
    ylim((min(unlist(plot_data))), 1.015*(max(unlist(plot_data)))) +
    labs(title = protein, size = 38) +
    theme(axis.text.y = element_text(size = 30), 
          axis.text.x = element_text(size = 34, colour = "black"),
          axis.title = element_blank(),
          plot.title = element_text(size = 38, hjust = .5, face = "bold")) +
    annotate("text", x = 1.5, y = ylev*(max(unlist(plot_data))), label = sig, size = sigsize, colour = "red2") +
    annotate("segment", x = 1, xend = 2, y = 1.004*(max(unlist(plot_data))), yend = 1.004*(max(unlist(plot_data))), colour = "red2", size = 1.6)
}




dot_plot_autosig <- function(data, protein) {
  p_val_colname <- paste0(condition_A, "-", condition_B, " adj p-value")
  pval <- data[data$Gene == protein, p_val_colname][1]
  SigLabel = case_when(
    pval < 0.0001 ~ "****",
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ NA_character_
  )
  plot_data <- data_formater(data = data, protein = protein)
  means <- colMeans(plot_data)
  sds <- apply(plot_data, 2, sd)
  p <- plot_data %>%
    ggparcoord(alphaLines = 0,
               showPoints = TRUE,
               boxplot = F,
               scale = "globalminmax") +     
    geom_point(aes(x = 1, y = means[1]), color = "grey50", size = 10) + 
    geom_errorbar(aes(x = 1, ymin = means[1] - sds[1], ymax = means[1] + sds[1]),
                  color = "grey50", width = 0.15, size = 1) + 
    geom_point(aes(x = 2, y = means[2]), color = "grey50", size = 10) + 
    geom_errorbar(aes(x = 2, ymin = means[2] - sds[2], ymax = means[2] + sds[2]),
                  color = "grey50", width = 0.15, size = 1) +
    geom_point(colour = "black", size = 5) +
    geom_line(colour = "grey22", size = 1, linetype = 2, alpha = .5) +
    theme_bw() +
    labs(title = protein, size = 38) +
    theme(axis.text.y = element_text(size = 30), 
          axis.text.x = element_text(size = 34, colour = "black"),
          axis.title = element_blank(),
          plot.title = element_text(size = 38, hjust = .5, face = "bold")) +
  if(!is.na(SigLabel)) {
    p <- p +
      annotate("text", x = 1.5, y = 1.005*(max(unlist(plot_data))), label = SigLabel, size = 20, colour = "red2") +
      annotate("segment", x = 1, xend = 2, y = 1.004*(max(unlist(plot_data))), yend = 1.004*(max(unlist(plot_data))), colour = "red2", size = 1.6) +
      ylim(0.997*(min(unlist(plot_data))), 1.007*(max(unlist(plot_data))))
  }
  return(p)
}




dot_plot(aim3_DE, "COL14A1", sig = "p-adj = 0.0508")
dot_plot(aim3_DE, "APCS", sig = "**", sigsize = 28, ylev=1.01)
dot_plot(aim3_DE, "VTN", sig = "*", sigsize = 28, ylev=1.01)




dot_plot_autosig(aim3_DE, "APCS")
dot_plot_autosig(aim3_DE, "COL6A1")




pdf("aim3_dot_plots_wSD_biggertext.pdf", width = 7, height = 10)

dot_plot(aim3_DE, "COL14A1", sig = "p-adj = 0.0508")
dot_plot(aim3_DE, "APCS", sig = "**", sigsize = 28, ylev=1.01)
dot_plot(aim3_DE, "VTN", sig = "*", sigsize = 28, ylev=1.01)

dev.off()








