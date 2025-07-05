
library(tidyverse)
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(scales)
library(grid)



# function:
make_PLF_graphs_PDF <- function(PDF_outname, df, condition1, condition2) {
  
  df_list <- split(df, df$GeneAC)  # split the big results table into a list of dfs (one per protein)
  graphs <- list()
  
  # Convert condition1 and condition2 to symbols
  condition1_sym <- sym(condition1)
  condition2_sym <- sym(condition2)
  # Create dynamic column names for Diff and p-values
  diff_colname <- paste0("Diff ", condition2, " vs ", condition1)
  p_val_colname <- paste0("p: ", condition2, " vs ", condition1)

  pdf(PDF_outname, width = 16, height = 10)
  
  for (gene in names(df_list)) {
    current_df <- df_list[[gene]]
    
    df_long_bars_condition1 <- current_df %>%
      select(Domain_Start, !!condition1_sym) %>%
      pivot_longer(cols = c(!!condition1_sym),
                   names_to = "Condition", values_to = "Value")
    
    df_long_bars_condition2 <- current_df %>%
      select(Domain_Start, !!condition2_sym) %>%
      pivot_longer(cols = c(!!condition2_sym),
                   names_to = "Condition", values_to = "Value")
    
    df_long_line <- current_df %>%
      select(Domain_Start, diff_colname, p_val_colname) %>%
      pivot_longer(cols = c(diff_colname),
                   names_to = "Condition", values_to = "Value") %>%
      mutate(
        SigLabel = case_when(
          Condition == diff_colname & !!sym(p_val_colname) < 0.0001 ~ "****",
          Condition == diff_colname & !!sym(p_val_colname) < 0.001 ~ "***",
          Condition == diff_colname & !!sym(p_val_colname) < 0.01 ~ "**",
          Condition == diff_colname & !!sym(p_val_colname) < 0.05 ~ "*",
          TRUE ~ NA_character_
        )
      )
    
    bar_plot_condition1 <- ggplot(df_long_bars_condition1, aes(x = Domain_Start, y = Value, fill = Value)) +
      geom_bar(colour = "grey30", size = 0.01, stat = "identity") + # just makes the height of each bar directly the value on the Y axis
      facet_grid(Condition ~ ., switch = "y") +  # facet wrap, but with the labels on the left
      scale_fill_gradient(low = "lightskyblue1", high = "navy") +
      scale_x_continuous(breaks = unique(df_long_bars_condition1$Domain_Start), expand = c(.01, .01)) +   # expand makes it so there's less gaps either side of plot
      scale_y_continuous(n.breaks = 8, labels = scientific_format(), trans = "log10") +
      guides(fill = "none") +
      labs(
        title = paste(gene),
        x = NULL,
        y = "Intensity Value",
      ) +
      theme_classic() +
      theme(
        axis.text.x = element_text(size = 10),  # x-axis text size
        axis.text.y = element_text(size = 10),  # y-axis text size
        strip.text = element_text(face = "bold", size = 12),   # facet labels bold & size
        strip.background = element_blank(),                   # removes the grey box from the facet labels
        strip.placement = "outside",               # puts facet labels outside the plot
        plot.title = element_text(size = 18, face = "bold"),       # Title text size
        axis.title.y = element_text(size = 14),                     # Y-axis label size
        panel.grid.major.y = element_line(color = "grey85", linetype = 2, size = .1),
        panel.grid.minor.y = element_line(color = "grey85", linetype = 2, size = .1)
      )
    
    bar_plot_condition2 <- ggplot(df_long_bars_condition2, aes(x = Domain_Start, y = Value, fill = Value)) +
      geom_bar(colour = "grey30", size = 0.01, stat = "identity") + # just makes the height of each bar directly the value on the Y axis
      facet_grid(Condition ~ ., switch = "y") +  # facet wrap, but with the labels on the left
      scale_fill_gradient(low = "lightskyblue1", high = "navy") +
      scale_x_continuous(breaks = unique(df_long_bars_condition2$Domain_Start), expand = c(.01, .01)) +   # expand makes it so there's less gaps either side of plot
      scale_y_continuous(n.breaks = 8, labels = scientific_format(), trans = "log10") +
      guides(fill = "none") +
      labs(
        x = "Domain",
        y = "Intensity Value",
      ) +
      theme_classic() +
      theme(
        axis.text.x = element_text(size = 10),  # x-axis text size
        axis.text.y = element_text(size = 10),  # y-axis text size
        strip.text = element_text(face = "bold", size = 12),   # facet labels bold & size
        strip.background = element_blank(),                   # removes the grey box from the facet labels
        strip.placement = "outside",               # puts facet labels outside the plot
        axis.title.x = element_text(size = 14),                    # X-axis label size
        axis.title.y = element_text(size = 14),                     # Y-axis label size
        panel.grid.major.y = element_line(color = "grey85", linetype = 2, size = .1),
        panel.grid.minor.y = element_line(color = "grey85", linetype = 2, size = .1)
      )
    
    
    # Line plot for Diff
    line_plot <- ggplot(df_long_line, aes(x = Domain_Start, y = Value)) +
      geom_line(color = "grey20", size = .1) +
      geom_point(color = "darkblue", size = 5.2, alpha = .8) +
      facet_grid(Condition ~ ., switch = "y") +  # facet wrap, but with the labels on the left
      scale_x_continuous(breaks = unique(df_long_bars_condition2$Domain_Start), expand = c(.03, .03)) +   # expand makes it so there's less gaps either side of plot
      scale_y_continuous(n.breaks = 5, labels = scientific_format()) +
      coord_cartesian(ylim = c(
        -1.2 * max(abs(df_long_line$Value), na.rm = TRUE),
        1.2 * max(abs(df_long_line$Value), na.rm = TRUE)
      )) +
      geom_text(data = subset(df_long_line, !is.na(SigLabel)),
                aes(label = SigLabel),
                vjust = -0.2, color = "black", size = 7) +
      labs(
        y = "Difference in Intensity per Domain",
        x = NULL
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),  # y-axis text size
        strip.text = element_text(face = "bold", size = 12),   # facet labels bold & size
        strip.background = element_blank(),                   # removes the grey box from the facet labels
        strip.placement = "outside",               # puts facet labels outside the plot
        axis.title.y = element_text(size = 14),                     # Y-axis label size
        panel.grid.major.y = element_line(color = "grey85", linetype = 2, size = .1),
        panel.grid.minor.y = element_line(color = "grey85", linetype = 2, size = .1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    
    
    # Combine the bar plots and line plot
    composite <- grid.arrange(bar_plot_condition1,
                              line_plot,
                              bar_plot_condition2,
                              ncol = 1,
                              heights = c(1, .8, 1)  # line_plot is 80% height of the others
                              )
    
    graphs[[gene]] <- composite
    
    grid.draw(composite)
    
  }
  
  

  dev.off()
  
  return(graphs)
  
}



#############################################################################################################


# read in PLF raw data (local output tsv, or raw data tsv downloaded off MPLF website)

plf_results_og <- fread("Aim2_EDTA_PLF_allprots_PLFresults og.tsv")
plf_results_new <- fread("aim2_allprots_PLFresults new.tsv")




make_PLF_graphs_PDF("PLF graphs EDTA og.pdf", plf_results_og, "Epidermis+", "Epidermis-")
make_PLF_graphs_PDF("PLF graphs EDTA new.pdf", plf_results_new, "Epidermis+", "Epidermis-")






og <- read.csv("aim2_plf_all_PLF_segment_coverage og.csv", check.names = F)
new <- read.csv("aim2_allprots_PLFresults_segment_coverage_summary_table new.csv", check.names = F)
new_WW <- read.csv("aim2_allprots_PLFresults_segment_coverage_summary_table new WW.csv", check.names = F)

all <- og %>%
  full_join(new, by = "GeneAC") %>%
  full_join(new_WW, by = "GeneAC")

write.csv(all, "all_versions_coverages.csv", row.names = F)










