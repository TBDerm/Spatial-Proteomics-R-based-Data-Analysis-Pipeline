
# this script will help you choose which type of normalisation and imputation to use in the analysis script

library(tidyverse)
library(data.table)
library(vsn) # for VSN normalisation
library(limma) # for limma normalisation
library(cowplot) # to display grids of graphs
library(ggpubr) # to display grids of graphs
library(mice) # for imputation




##################################################################################################################
# need to edit this part to suit your data
##################################################################################################################

# reading intensity data (DIA-NN output or similar)
protein_dataset <- fread('Aim3_photoage_N3_protein_dataset_all+ECM+BM_noted.csv') # proteins
peptide_dataset <- fread('Aim3_photoage_N3_peptide_dataset_all+ECM+BM_noted.csv') # peptides

# edit sample names and experimental groups (must be in correct order, must match) (each sample name must be unique)
# number of sample and group names given must match the number of sample columns in data
# (these will be used for graphs so ensure they're spelt right etc)
sample_names <- c("0841_Buttock", "0842_Buttock", "0863_Buttock",
                  "0841_Forearm", "0842_Forearm", "0863_Forearm")
# format should be rep("group name", x) - where x is the number of samples (columns) in that experimental group
group_names <- c(rep("Buttock", 3),
                 rep("Forearm", 3))

# enter the number range for the sample columns in protein and peptide data
# e.g. if columns 5, 6, 7, 8, 9, and 10 are my samples (6 in total), enter 5:10 here.
prot_sample_cols <- 5:10  # in protein data
pep_sample_cols <- 11:16  # in peptide data

# enter the number for the column containing uniprot AC numbers (in DIANN data, this is the Protein.Group column)
prot_uniprot_AC_col <- 1  # in protein data
pep_uniprot_AC_col <- 1  # in peptide data

# enter the number for the column containing Genes
prot_gene_col <- 3  # in protein data
pep_gene_col <- 4  # in peptide data

# enter the number for the column containing peptide sequence (stripped.sequence in DIANN)
pep_sequence_col <- 7 # in peptide data


##################################################################################################################
# just run the below section to tidy data etc
##################################################################################################################

metadata <- data.frame(Sample = sample_names, Group = group_names)

# ensuring sample columns are in num format, and renaming important columns to suit code
protein_dataset <- protein_dataset %>% 
  mutate(across(all_of(prot_sample_cols), as.numeric)) %>%
  rename(UniProt_AC = all_of(prot_uniprot_AC_col)) %>%
  rename(Gene = all_of(prot_gene_col))
peptide_dataset <- peptide_dataset %>% 
  mutate(across(all_of(pep_sample_cols), as.numeric)) %>%
  rename(UniProt_AC = all_of(pep_uniprot_AC_col)) %>%
  rename(Gene = all_of(pep_gene_col)) %>%
  rename(Peptide_Sequence = all_of(pep_sequence_col))

# storing the original Uniprot_AC and Genes columns, with multiple codes, in copy columns called 'multiple_ACs' and 'multiple_genes' (if not already present)
if (!"multiple_ACs" %in% colnames(protein_dataset)) {
  protein_dataset <- protein_dataset %>% mutate("multiple_ACs" = UniProt_AC)
}
if (!"multiple_genes" %in% colnames(protein_dataset)) {
  protein_dataset <- protein_dataset %>% mutate("multiple_genes" = Gene)
}
if (!"multiple_ACs" %in% colnames(peptide_dataset)) {
  peptide_dataset <- peptide_dataset %>% mutate("multiple_ACs" = UniProt_AC)
}
if (!"multiple_genes" %in% colnames(peptide_dataset)) {
  peptide_dataset <- peptide_dataset %>% mutate("multiple_genes" = Gene)
}

# remove rows that do not have a Gene (these are non-human contaminants)
# and taking only the first Uniprot_AC from the Uniprot_AC column (generally when more than one is given, the first is the main one)
# same for genes - taking only the first
protein_dataset <- protein_dataset[!(protein_dataset$Gene == "" | is.na(protein_dataset$Gene)), ] %>%
  separate(UniProt_AC, into = c("UniProt_AC", "other"), sep = ";") %>%
  separate(Gene, into = c("Gene", "otherG"), sep = ";") %>%
  mutate(Gene = toupper(Gene)) %>%
  select(-other, -otherG)
peptide_dataset <- peptide_dataset[!(peptide_dataset$Gene == "" | is.na(peptide_dataset$Gene)), ] %>%
  separate(UniProt_AC, into = c("UniProt_AC", "other"), sep = ";") %>%
  separate(Gene, into = c("Gene", "otherG"), sep = ";") %>%
  mutate(Gene = toupper(Gene)) %>%
  select(-other, -otherG)

# at this stage should filter by Q value, to remove any proteins >0.05 (protein and peptide datasets)
# note this is not needed for DIA-NN data (is already done)
### ** idk what the column would be called - change this
# protein_dataset <- protein_dataset %>% filter(Global.Q < 0.05)
# peptide_dataset <- peptide_dataset %>% filter(Global.Q < 0.05)

# rename sample columns
names(protein_dataset)[prot_sample_cols] <- sample_names
names(peptide_dataset)[pep_sample_cols] <- sample_names

# any 0 value in the protein / peptide data is replaced by NA
protein_dataset[protein_dataset == 0] <- NA
peptide_dataset[peptide_dataset == 0] <- NA

# appending any duplicate uniprot_AC rows in the data with -1, -2, -3 etc (original ACs are preserved in the 'multiple AC' column)
protein_dataset$UniProt_AC <- ave(protein_dataset$UniProt_AC, protein_dataset$UniProt_AC, FUN = function(x) make.unique(as.character(x), sep = "-"))

# Count unique peptide sequences per protein, and add this to protein and peptide datasets
peptide_counts <- peptide_dataset %>% select(UniProt_AC, Peptide_Sequence) %>% distinct() %>%
  dplyr::count(UniProt_AC, name = "Peptide_Count")
protein_dataset <- protein_dataset %>%
  left_join(peptide_counts, by = "UniProt_AC")
peptide_dataset <- peptide_dataset %>%
  left_join(peptide_counts, by = "UniProt_AC")

# exclude any proteins identified by only one peptide, from protein and peptide datasets
protein_dataset_reduced <- protein_dataset %>% 
  filter(Peptide_Count > 1) %>% 
  as.data.frame()
peptide_dataset_reduced <- peptide_dataset %>% 
  filter(Peptide_Count > 1) %>% 
  as.data.frame()

# getting just the sample columns (intensity values)
sample_intensities <- protein_dataset_reduced[, prot_sample_cols] 
rownames(sample_intensities) <- protein_dataset_reduced$UniProt_AC # getting UniProt_AC as rownames
sample_intensities[sample_intensities == 0] <- NA # any 0 value in the whole dataset is replaced with NA

# removing proteins with >40% missing values, per group (thus keeping only proteins which are present in at least 60% of samples in every group)
unique_groups <- unique(metadata$Group)
prot_names <- rownames(sample_intensities)
# Calculate missing value proportions per group for each gene
groupwise_missing_proportion <- sapply(unique_groups, function(group) {
  sampless <- metadata$Sample[metadata$Group == group]
  rowMeans(is.na(sample_intensities[, sampless, drop = FALSE]))  # Proportion of NA values
})
# Identify proteins which are present in >=60% of samples in every group (i.e. excluding prots where ANY group has >=40% missing values)
prots_to_keep <- rownames(sample_intensities)[apply(groupwise_missing_proportion, 1, max) <= 0.4]
# Subset the dataframe to keep only valid genes
sample_intensities_filtered <- sample_intensities[prots_to_keep, ]
# Ensure row names are preserved
rownames(sample_intensities_filtered) <- prots_to_keep

# Log2-transform the reduced protein data
log2_intensity_matrix <- log2(as.matrix(sample_intensities_filtered)) # creates matrix of log2-transformed values

# function to make data into 'long' format:
data_longer <- function(x) x %>%
  mutate(UniProt_AC = rownames(x)) %>%
  pivot_longer(cols = sample_names,
               names_to = "Sample",
               values_to = "Intensity") %>%
  left_join(metadata, by = "Sample") %>%
  filter(!is.na(Intensity))

# functions to plot density curves & histograms:
plot_density <- function(data, t) data %>%
  ggplot(aes(x = Intensity, colour = Sample)) +
  geom_density(size = 1) +
  labs(x = "Log2 Intensity", y = "Frequency Density", title = t) +
  theme_bw()
plot_histo <- function(data, title) data %>%
  ggplot(aes(x = Intensity)) +
  geom_histogram(binwidth = .5, alpha = .5, color = "blue", fill = "blue") +
  labs(x = "Log2 Intensity", y = "Count", title = title) +
  theme_bw()


# normalisation function
data_normaliser <- function(data, method) {
  # Median Absolute Deviation Normalisation
  if (method == "median_abs_dev") {
    return(as.data.frame(normalizeMedianAbsValues(data)))
  }
  # VSN Normalisation
  if (method == "vsn") {
    return(as.data.frame(justvsn(as.matrix(data), minDataPointsPerStratum = 3)))
  }
  # Median Centering (makes all samples have the same median)
  if (method == "median_centred") {
    return(as.data.frame(sweep(data, 2, apply(data, 2, median, na.rm = TRUE), "/")))
  }
  # Cyclic Loess Normalization (limma)
  if (method == "cyclic_loess") {
    return(as.data.frame(normalizeCyclicLoess(as.matrix(data))))
  }
  # Robust Linear Scaling (RLS)
  if (method == "rls") {
    return(as.data.frame(sweep(data, 2, apply(data, 2, median, na.rm = TRUE), "-") /
                           apply(data, 2, mad, na.rm = TRUE)))
  }
  # Median-Based Scaling
  if (method == "median_scaling") {
    return(as.data.frame(sweep(data, 2, apply(data, 2, median, na.rm = TRUE), "-") /
                           apply(data, 2, mad, na.rm = TRUE)))
  }
  # median centring then median absolute deviation scaling
  if (method == "median_centring_MAD_scaling") {
    MC <- as.data.frame(sweep(data, 2, apply(data, 2, median, na.rm = TRUE), "-"))
    return(as.data.frame(sweep(MC, 2, apply(MC, 2, mad, na.rm = TRUE), "/")))
  }
  # median centring then SD (aka z-score) scaling
  if (method == "median_centring_SD_scaling") {
    MC <- as.data.frame(sweep(data, 2, apply(data, 2, median, na.rm = TRUE), "-"))
    return(as.data.frame(sweep(MC, 2, apply(MC, 2, sd, na.rm = TRUE), "/")))
  }
}

# imputation function
data_imputer <- function(data, method, k = 3) {
  imputed_data_list <- list()
  for (group in unique_groups) {
    sampless <- metadata$Sample[metadata$Group == group]
    group_data <- data[, sampless, drop = FALSE]
    # Impute the group data based on the method specified in function
    if (method %in% c("pmm", "cart", "lasso.norm")) {
      # Use MICE for PMM, CART, or Lasso
      imputed_group_data <- complete(mice(group_data, method = method))
      imputed_group_data <- imputed_group_data %>% select(all_of(sampless))
    } else if (method == "knn") {
      # Use kNN from VIM package
      library(VIM)
      imputed_group_data <- kNN(group_data, k = k)
      imputed_group_data <- imputed_group_data %>% select(all_of(sampless))
    }
    imputed_data_list[[group]] <- imputed_group_data
  }
  imputed_data <- do.call(cbind, imputed_data_list)
  colnames(imputed_data) <- colnames(data)
  return(imputed_data)
}



##################################################################################################################
# normalisation options
##################################################################################################################

# data has been filtered and log2-transformed (now only contains proteins identified by >1 peptide, and proteins present in at least 60% of samples in every group)

# run the below lines to create normalised dataset options

# original non-normalised data (just log2-transformed)
original <- as.data.frame(log2_intensity_matrix)
# median absolute deviation normalisation
median_abs_dev <- data_normaliser(log2_intensity_matrix, method = "median_abs_dev")
# VSN
vsn_norm <- data_normaliser(log2_intensity_matrix, method = "vsn")
# median centring (makes all samples have the same median)
median_centred <- data_normaliser(log2_intensity_matrix, method = "median_centred")


### new:
# Cyclic Loess Normalization (limma)
cyclic_loess_norm <- data_normaliser(log2_intensity_matrix, method = "cyclic_loess")
# Robust Linear Scaling (RLS)
rls_norm <- data_normaliser(log2_intensity_matrix, method = "rls")
# Median-Based Scaling
median_scaled <- data_normaliser(log2_intensity_matrix, method = "median_scaling")
# median centring then median absolute deviation scaling
MC_MAD_scaled <- data_normaliser(log2_intensity_matrix, method = "median_centring_MAD_scaling")
# median centring then SD (aka z-score) scaling
sd_scaled <- data_normaliser(log2_intensity_matrix, method = "median_centring_SD_scaling")





# graphs to compare normalised datasets: -------------------
# run each section to view a grid of the graphs

plot_grid(
  plot_density(data_longer(original), "original log2-transformed"),
  plot_density(data_longer(median_abs_dev), "median absolute deviation"),
  plot_density(data_longer(vsn_norm), "variance stabilising"),
  plot_density(data_longer(median_centred), "median centring"),
  plot_density(data_longer(cyclic_loess_norm), "cyclic_loess_norm"),
  plot_density(data_longer(rls_norm), "rls_norm"),
  plot_density(data_longer(median_scaled), "median_scaled"),
  plot_density(data_longer(MC_MAD_scaled), "MC_MAD_scaled"),
  plot_density(data_longer(sd_scaled), "sd_scaled")
)

ggarrange(
  meanSdPlot(as.matrix(original))$gg,
  meanSdPlot(as.matrix(median_abs_dev))$gg,
  meanSdPlot(as.matrix(vsn_norm))$gg,
  meanSdPlot(as.matrix(median_centred))$gg,
  meanSdPlot(as.matrix(cyclic_loess_norm))$gg,
  meanSdPlot(as.matrix(rls_norm))$gg,
  meanSdPlot(as.matrix(median_scaled))$gg,
  meanSdPlot(as.matrix(MC_MAD_scaled))$gg,
  meanSdPlot(as.matrix(sd_scaled))$gg,
  labels = c("original log2-transformed", "median absolute deviation", "variance stabilising", "median centring",
             "cyclic loess", "robust linear scaling (rls)", "median_scaled", "MC_MAD_scaled", "sd_scaled"), 
  ncol = 3, nrow = 3
)


# choose the normalisation where density curves are similar across samples, 
# and meanSDplot red line is level (not trending upwards) with uniform distribution of dots


# put chosen normalised data here (options: median_abs_dev, vsn_norm, median_centred, limma_norm)
normalised_intensity_df <- cyclic_loess_norm





##################################################################################################################
# imputation options 
##################################################################################################################

# (imputation done separately per group)
# run the below lines to create imputed dataset options

# predictive mean matching:
imputed_pmm <- data_imputer(normalised_intensity_df, method = "pmm")
# classification and regression trees:
imputed_cart <- data_imputer(normalised_intensity_df, method = "cart")
# lasso linear regression:
imputed_lasso <- data_imputer(normalised_intensity_df, method = "lasso.norm")
# K nearest neighbours:   (can specify k, i.e. the number of rows compared - but 3 should be good bc is all protein abundance from same MS run)
imputed_knn <- data_imputer(normalised_intensity_df, method = "knn")



# plots to compare -------------------------------------------------------------------------

plot_grid(
  plot_histo(data_longer(normalised_intensity_df), "original"),  # original non-imputed dataset (i.e. does not count NA values)
  plot_histo(data_longer(imputed_pmm), "pmm"),
  plot_histo(data_longer(imputed_cart), "cart"),
  plot_histo(data_longer(imputed_lasso), "lasso"),
  plot_histo(data_longer(imputed_knn), "knn")
)

plot_grid(
  plot_density(data_longer(normalised_intensity_df), "original"),  # original non-imputed dataset (i.e. does not count NA values)
  plot_density(data_longer(imputed_pmm), "pmm"),
  plot_density(data_longer(imputed_cart), "cart"),
  plot_density(data_longer(imputed_lasso), "lasso"),
  plot_density(data_longer(imputed_knn), "knn")
)



# choose the imputation method which produces a histogram and curves most similar to the original data




