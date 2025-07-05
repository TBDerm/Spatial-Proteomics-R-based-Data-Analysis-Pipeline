
# important: venns work with a max of 4 categories :(



library(tidyverse)
library(ggvenn)
library(gplots)
library(dplyr)
library(data.table)
library(writexl)



proteins_allprots <- fread("valid_proteins_allprots_libs_venn.csv")
proteins_ECMprots <- fread("valid_proteins_ECMprots_libs_venn.csv")

peptides_allprots <- fread("valid_peptides_allprots_libs_venn.csv")
peptides_ECMprots <- fread("valid_peptides_ECMprots_libs_venn.csv")




# function: frag speclibs only
frag_only <- function(x) {
  x %>%
    select(FRAGPIPE_G, FRAGPIPE_T)
}

# function: diann speclibs only
diann_only <- function(x) {
  x %>%
    select(DIANN_noMods, DIANN_mods_noTune, DIANN_mods_Gtuned, DIANN_mods_Ttuned)
}

# function: diann G tuned v T tuned v no tuned
tunes_comp <- function(x) {
  x %>%
    select(DIANN_mods_noTune, DIANN_mods_Gtuned, DIANN_mods_Ttuned)
}

# function: diann v frag (best of)
diann_v_frag <- function(x) {
  x %>%
    select(DIANN_noMods, DIANN_mods_noTune, FRAGPIPE_G, FRAGPIPE_T)
}




#######################



# function to create data for venn diagram
make_venn_data <- function(x) {
  gene_list <- lapply(x, function(y) unique(y[y != ""]))
  setNames(as.list(gene_list), colnames(x))
}

# plot venn function
plot_venn <- function(data, venn_colours, title) {
  data %>%
    ggvenn(fill_color = venn_colours, fill_alpha = .5,
           stroke_size = .5, set_name_size = 4, text_size = 5) +
    labs(title = title) +
    theme(plot.title = element_text(size = 16, hjust = 0.5))
}

# function to make df of genes in each venn section
create_venn_dataframe <- function(protein_lists) {
  unique_genes <- unique(unlist(protein_lists))
  binary_matrix <- as.data.frame(sapply(protein_lists, function(set) {
    unique_genes %in% set  # TRUE/FALSE presence check
  }))
  binary_matrix <- as.data.frame(lapply(binary_matrix, as.integer))
  colnames(binary_matrix) <- names(protein_lists)
  rownames(binary_matrix) <- unique_genes
  v.table <- venn(binary_matrix)
  # Extract intersections
  intersections <- attr(v.table, "intersections")
  # Find the maximum length of the lists in the intersections
  max_length <- max(sapply(intersections, length))
  # Process each list: sort it alphabetically and pad with "" to match the max length
  intersections_processed <- lapply(intersections, function(x) {
    # Sort the list alphabetically
    x_sorted <- sort(x, na.last = TRUE)  # Sort and put NA values at the end
    # Extend the list with blank strings "" to match the maximum length
    length(x_sorted) <- max_length
    return(x_sorted)
  })
  # Convert the processed list to a dataframe
  df_intersections <- as.data.frame(intersections_processed, stringsAsFactors = FALSE)
  # Replace NAs with blank strings
  df_intersections[is.na(df_intersections)] <- ""
  return(df_intersections)
}




######################################################################


# noMods    noTune      Gtune     Ttune      Gfrag       Tfrag
c("green3", "skyblue2", "violet", "yellow", "darkblue", "orange")

plot_venn(make_venn_data(proteins_allprots), c("green3", "skyblue2", "violet", "darkblue"), "All proteins")


# plots:
plot_venn(make_venn_data(frag_only(proteins_allprots)), c("darkblue", "orange"), "All proteins, FRAGPIPE speclibs")
plot_venn(make_venn_data(diann_only(proteins_allprots)), c("green3", "skyblue2", "violet", "yellow"), "All proteins, DIANN speclibs")
plot_venn(make_venn_data(tunes_comp(proteins_allprots)), c("skyblue2", "violet", "yellow"), "All proteins, DIANN tuning comparison")
plot_venn(make_venn_data(diann_v_frag(proteins_allprots)), c("green3", "skyblue2", "darkblue", "orange"), "All proteins, best of DIANN vs FRAGPIPE")

plot_venn(make_venn_data(frag_only(proteins_ECMprots)), c("darkblue", "orange"), "ECM proteins, FRAGPIPE speclibs")
plot_venn(make_venn_data(diann_only(proteins_ECMprots)), c("green3", "skyblue2", "violet", "yellow"), "ECM proteins, DIANN speclibs")
plot_venn(make_venn_data(tunes_comp(proteins_ECMprots)), c("skyblue2", "violet", "yellow"), "ECM proteins, DIANN tuning comparison")
plot_venn(make_venn_data(diann_v_frag(proteins_ECMprots)), c("green3", "skyblue2", "darkblue", "orange"), "ECM proteins, best of DIANN vs FRAGPIPE")

plot_venn(make_venn_data(frag_only(proteins_BMprots)), c("darkblue", "orange"), "BM proteins, FRAGPIPE speclibs")
plot_venn(make_venn_data(diann_only(proteins_BMprots)), c("green3", "skyblue2", "violet", "yellow"), "BM proteins, DIANN speclibs")
plot_venn(make_venn_data(tunes_comp(proteins_BMprots)), c("skyblue2", "violet", "yellow"), "BM proteins, DIANN tuning comparison")
plot_venn(make_venn_data(diann_v_frag(proteins_BMprots)), c("green3", "skyblue2", "darkblue", "orange"), "BM proteins, best of DIANN vs FRAGPIPE")



# venn intersection data:
allprots_frag_comp <- create_venn_dataframe(frag_only(proteins_allprots))
allprots_diann_comp <- create_venn_dataframe(diann_only(proteins_allprots))
allprots_tuning_comp <- create_venn_dataframe(tunes_comp(proteins_allprots))
allprots_diann_v_frag <- create_venn_dataframe(diann_v_frag(proteins_allprots))

ECMprots_frag_comp <- create_venn_dataframe(frag_only(proteins_ECMprots))
ECMprots_diann_comp <- create_venn_dataframe(diann_only(proteins_ECMprots))
ECMprots_tuning_comp <- create_venn_dataframe(tunes_comp(proteins_ECMprots))
ECMprots_diann_v_frag <- create_venn_dataframe(diann_v_frag(proteins_ECMprots))
                                               
BMprots_frag_comp <- create_venn_dataframe(frag_only(proteins_BMprots))
BMprots_diann_comp <- create_venn_dataframe(diann_only(proteins_BMprots))
BMprots_tuning_comp <- create_venn_dataframe(tunes_comp(proteins_BMprots))
BMprots_diann_v_frag <- create_venn_dataframe(diann_v_frag(proteins_BMprots))


all_proteins <- create_venn_dataframe(proteins_allprots)
ECM_proteins <- create_venn_dataframe(proteins_ECMprots)


########################################################################################




pdf("all_venns_peptides_G.pdf", width = 10, height = 8)

# plots here
plot_venn(make_venn_data(proteins_allprots), c("green3", "skyblue2", "violet", "darkblue"), "All proteins")
plot_venn(make_venn_data(proteins_ECMprots), c("green3", "skyblue2", "violet", "darkblue"), "ECM proteins")
plot_venn(make_venn_data(peptides_allprots), c("green3", "skyblue2", "violet", "darkblue"), "All peptides")
plot_venn(make_venn_data(peptides_ECMprots), c("green3", "skyblue2", "violet", "darkblue"), "ECM peptides")


dev.off()







df_list <- list(all_proteins = all_proteins,
                ECM_proteins = ECM_proteins)


# Write to Excel file
write_xlsx(df_list, path = "all_venns_proteins_intersection_data.xlsx")







