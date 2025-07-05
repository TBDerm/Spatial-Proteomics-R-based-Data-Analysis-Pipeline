
library(tidyverse)
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)




# read in PLF raw data (local output tsv, or raw data tsv downloaded off MPLF website)

plf_results <- fread("Aim2_EDTA_PLF_BMprots_PLFresults.tsv")

one_protein <- fread("one_protein.csv")


##### to make 2 bars, and a line for diff:

df_long_bars_condition1 <- one_protein %>%
  select(Domain_Start, Buttock) %>%
  pivot_longer(cols = c(Buttock),
               names_to = "Condition", values_to = "Value")

df_long_bars_condition2 <- one_protein %>%
  select(Domain_Start, Forearm) %>%
  pivot_longer(cols = c(Forearm),
               names_to = "Condition", values_to = "Value")




# Bar plots
bar_plot_condition1 <- ggplot(df_long_bars_condition1, aes(x = Domain_Start, y = Value, fill = Value)) +
  geom_bar(colour = "grey30", size = 0.01, stat = "identity") + # just makes the height of each bar directly the value on the Y axis
  facet_grid(Condition ~ ., switch = "y") +  # facet wrap, but with the labels on the left
  scale_fill_gradient(low = "lightskyblue1", high = "navy") +
  scale_x_continuous(breaks = unique(df_long_bars_condition1$Domain_Start), expand = c(.01, .01)) +   # expand makes it so there's less gaps either side of plot
  scale_y_continuous(n.breaks = 8, labels = scientific_format(), trans = "log10") +
  guides(fill = "none") +
  labs(
    title = "*** beep protein AC ***",
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
  scale_fill_gradient(low = "pink", high = "magenta4") +
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



bar_plot_condition1
bar_plot_condition2







df_long_line <- one_protein %>%
  select(Domain_Start, `Diff Forearm vs Buttock`, `p: Forearm vs Buttock`) %>%
  pivot_longer(cols = c(`Diff Forearm vs Buttock`),
               names_to = "Condition", values_to = "Value") %>%
  mutate(
    SigLabel = case_when(
      Condition == "Diff Forearm vs Buttock" & `p: Forearm vs Buttock` < 0.0001 ~ "****",
      Condition == "Diff Forearm vs Buttock" & `p: Forearm vs Buttock` < 0.001 ~ "***",
      Condition == "Diff Forearm vs Buttock" & `p: Forearm vs Buttock` < 0.01 ~ "**",
      Condition == "Diff Forearm vs Buttock" & `p: Forearm vs Buttock` < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  )


# Line plot for "Diff Forearm vs Buttock"
line_plot <- ggplot(df_long_line, aes(x = Domain_Start, y = Value)) +
  geom_line(color = "grey20", size = .1) +
  geom_point(color = "magenta3", size = 6) +
  facet_grid(Condition ~ ., switch = "y") +  # facet wrap, but with the labels on the left
  scale_x_continuous(breaks = unique(df_long_bars_condition2$Domain_Start), expand = c(.03, .03)) +   # expand makes it so there's less gaps either side of plot
  scale_y_continuous(n.breaks = 5, labels = scientific_format()) +
  ylim(-(max(df_long_line$Value)+0.1*(max(df_long_line$Value))), (max(df_long_line$Value)+0.1*(max(df_long_line$Value)))) +
  geom_text(data = subset(df_long_line, !is.na(SigLabel)),
    aes(label = SigLabel),
    vjust = -0.3, color = "black", size = 7) +
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

line_plot



# Combine the bar plots and line plot
grid.arrange(bar_plot_condition1,
             line_plot,
             bar_plot_condition2, 
             ncol = 1,
             heights = c(1, .8, 1)  # line_plot is 80% height of the others
)












