# This script takes in protein and peptide MS intensity datasets
# refines data, creates some lovely visualisations, and carries out PCA and limma analysis
# important: venns wont work for >4 experimental groups :( soz 


# load required packages (in order of use)
library(tidyverse)
library(data.table)
library(naniar) # for vis miss plot
library(dplyr)
library(scales)
library(limma)
library(gplots)
library(matrixStats)
library(vsn)
library(mice)
library(ggvenn) 
library(ggfortify)
library(ggrepel)
library(reshape2)
library(rlang)
library(gridExtra)




##################################################################################################################
# need to edit this part to suit your data
##################################################################################################################

# reading intensity data (DIA-NN output or similar)
protein_dataset <- fread('Aim2_EDTA_protein_dataset_all+ECM+BM_noted.csv') # proteins
peptide_dataset <- fread('Aim2_EDTA_peptide_dataset_all+ECM+BM_noted.csv') # peptides

# edit sample names and experimental groups (must be in correct order, must match) (each sample name must be unique)
# number of sample and group names given must match the number of sample columns in data
# (these will be used for graphs so ensure they're spelt right etc)
sample_names <- c("Epi+_01", "Epi+_02",
                  "Epi-_03", "Epi-_04")
# format should be rep("group name", x) - where x is the number of samples (columns) in that experimental group
group_names <- c(rep("Epidermis+", 2),
                 rep("Epdiermis-", 2))

# enter the number range for the sample columns in protein and peptide data
# e.g. if columns 5, 6, 7, 8, 9, and 10 are my samples (6 in total), enter 5:10 here.
prot_sample_cols <- 5:8  # in protein data
pep_sample_cols <- 11:14  # in peptide data

# enter the number for the column containing uniprot AC numbers (in DIANN data, this is the Protein.Group column)
prot_uniprot_AC_col <- 1  # in protein data
pep_uniprot_AC_col <- 1  # in peptide data

# enter the number for the column containing Genes
prot_gene_col <- 3  # in protein data
pep_gene_col <- 4  # in peptide data

# enter the number for the column containing peptide sequence (stripped.sequence in DIANN)
pep_sequence_col <- 7 # in peptide data

# specify the type of normalisation required (options: "none", "median_abs_dev", "vsn", "median_centred", "cyclic_loess", 
#                                                      "rls", "median_scaling", "median_centring_MAD_scaling", "median_centring_SD_scaling")
normalisation_type <- "cyclic_loess"

# specify the type of imputation required (options: "pmm", "cart", "lasso.norm", "knn")
imputation_type <- "knn"

# the script will make a new folder within current directory to store output data - provide a name for the folder
folder_name <- "sample_R_results_EDTA_cyclicloess_knn"

# output files will be named with this at the beginning - provide a short experiment ID
experiment_identifier <- "sample"

# choose a colour scheme for the plots - one colour per experimental group
# the number of colours set must match the number of experimental groups!
# (this can be changed per plot at the end if wanted)
colour_scheme <- c("darkorchid1", "seagreen3")
# colors() ## un-comment this line and run to see a list of all colours available
# aim 2: c("magenta3", "orange2")
# aim 3: c("limegreen", "blue")

# colours for heatmaps: low, high   (can be changed at bottom too)
heat_colours <- c("lemonchiffon2", "midnightblue")


##################################################################################################################
# the rest can just be run
# export datasets and graphs at the bottom
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

# Taking the uniprot AC and gene columns from the protein dataset, to make a key for use later
gene_AC_key <- protein_dataset %>% select(UniProt_AC, Gene)

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

# combining the reduced and filtered sample values (so now only valid proteins) with the full dataset to export
full_sample_intensities_filtered <- sample_intensities_filtered %>%
  mutate(UniProt_AC = rownames(sample_intensities_filtered)) %>%
  left_join((protein_dataset %>% select(-all_of(prot_sample_cols))), by = "UniProt_AC") %>%
  select(UniProt_AC, everything())
# peptide dataset, but only for valid proteins (proteins identified by >1 peptide, and with no more than 40% missing values per group)
valid_peptide_dataset <- peptide_dataset_reduced %>%
  filter(UniProt_AC %in% prots_to_keep)


# normalisation function
data_normaliser <- function(data, method) {
  # NO Normalisation
  if (method == "none") {
    return(as.data.frame(data))
  }
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

# imputation function (imputing per group):
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

# function to make data into 'long' format:
data_longer <- function(x) x %>%
  mutate(UniProt_AC = rownames(x)) %>%
  pivot_longer(cols = sample_names,
               names_to = "Sample",
               values_to = "Intensity") %>%
  left_join(metadata, by = "Sample") %>%
  filter(!is.na(Intensity))


# Log2-transform and normalise the reduced and filtered protein data (using chosen normalisation method)
log2_intensity_matrix <- log2(as.matrix(sample_intensities_filtered)) # creates matrix of log2-transformed values
normalised_intensity_df <- as.data.frame(data_normaliser(log2_intensity_matrix, method = normalisation_type))

# combining it with the full dataset to export
full_norm_prot_dataset <- normalised_intensity_df %>%
  mutate(UniProt_AC = rownames(normalised_intensity_df)) %>%
  left_join((protein_dataset %>% select(-all_of(prot_sample_cols))), by = "UniProt_AC") %>%
  select(UniProt_AC, everything())

# imputation, using chosen type:
imputed_df <- data_imputer(normalised_intensity_df, method = imputation_type)
rownames(imputed_df) <- rownames(normalised_intensity_df)

# combining it with the full dataset to export
full_imputed_dataset <- as.data.frame(imputed_df) %>%
  mutate(UniProt_AC = rownames(imputed_df)) %>%
  left_join((protein_dataset %>% select(-all_of(prot_sample_cols))), by = "UniProt_AC") %>%
  select(UniProt_AC, everything())

# long imputed dataset, to make plots (with genes too)
long_imputed_df <- data_longer(imputed_df) %>%
  left_join(protein_dataset_reduced, by = "UniProt_AC") %>%
  select(UniProt_AC, Gene, Sample, Intensity, Group)
long_imputed_df$Group <- factor(long_imputed_df$Group)

# peptide count per protein (valid proteins only), by experimental group
long_valid_peptides <- valid_peptide_dataset %>%
  select(-Peptide_Count) %>%
  pivot_longer(cols = sample_names,
               names_to = "Sample",
               values_to = "Intensity") %>%
  left_join(metadata, by = "Sample") %>%
  filter(!is.na(Intensity))
peptide_counts_per_group <- long_valid_peptides %>%
  select(UniProt_AC, Peptide_Sequence, Group) %>%
  distinct() %>%
  dplyr::count(UniProt_AC, Group, name = "Peptide_Count_byGroup")
peptide_counts_per_group_wGenes <- left_join(peptide_counts_per_group, gene_AC_key, by = "UniProt_AC")




# general plots ---------------------------------------------------------------------------------------------------------

# most plots use the imputed dataset (long_imputed_df) (values are log2-transformed, normalised, valid proteins only, imputed per group)

# missing values plot of all proteins identified by >1 peptide
vis_miss_allprots <- vis_miss(sample_intensities, cluster = T)
# missing values plot of just valid proteins (before imputation ofc)
vis_miss_validprots <- vis_miss(normalised_intensity_df, cluster = T)

# standard heatmap
heatmap_allprots <- function(colours = heat_colours) {
  long_imputed_df %>%
    ggplot(aes(Sample, UniProt_AC, fill = Intensity)) +
    geom_tile() +
    labs(title = "Heatmap of Log2 Intensity Values",
         x = "Samples",
         y = "Protein") +
    scale_fill_gradientn(colours = colorRampPalette(colours)(100), 
                         name = "Log2( Intensity )", na.value = "white") +
    scale_y_discrete(labels = NULL) + 
    theme_bw() +
    facet_wrap(~ Group, scales = "free_x")
}

# creating sub-dataset for top 30 valid proteins (by highest median intensity)
df_m <- imputed_df
df_m$UniProt_AC <- rownames(df_m)
df_m_num <- as.matrix(select(df_m, where(is.numeric)))
df_m$median_value <- rowMedians(df_m_num, na.rm = TRUE)
top_30_prots <- df_m %>%
  arrange(desc(median_value)) %>%
  slice_head(n = 30) %>%
  select(UniProt_AC, everything(), -median_value)
top_30_prots <- top_30_prots %>% left_join(gene_AC_key, by = "UniProt_AC") %>% select(-UniProt_AC)
rownames(top_30_prots) <- top_30_prots$Gene
top_30_prots <- select(top_30_prots, -Gene)

# standard heatmap for top 30:
heatmap_top30_prots <- function(colours = heat_colours) {
  data_longer(top_30_prots) %>%
    ggplot(aes(Sample, UniProt_AC, fill = Intensity)) +
    geom_tile() +
    labs(title = "Heatmap of Top 30 Most Abundant Proteins",
         x = "Samples",
         y = "Protein") +
    scale_fill_gradientn(colours = colorRampPalette(colours)(100), 
                         name = "Log2( Intensity )", na.value = "white") +
    theme_bw() +
    facet_wrap(~ Group, scales = "free_x")
}

# z-scored heatmap of all valid proteins
par(cex.main = 0.8) # Set title size before plotting
matrix_prot_intensity <- as.matrix(df_m_num)
zscore_heatmap_allprots <- function(colours = heat_colours) {
  heatmap.2(matrix_prot_intensity,
            scale = "column",              
            Colv = T,           
            Rowv = NA,
            col = colorRampPalette(colours)(100),
            trace = "none",
            key.title = NA,
            srtCol = 390,
            sepcolor = "black",
            colsep = (0:50),
            sepwidth = .001,
            labRow = NA,        
            adjCol = c(.9,.3),
            key.ylab = NA,
            cexCol = 1.2,    # Smaller column (x-axis) labels
            main = "Z-Scored Heatmap of Intensity Values",
  )
}

# z-scored heatmap of top 30
matrix_top30_prot_intensity <- as.matrix(top_30_prots)
zscore_heatmap_top30_prots <- function(colours = heat_colours) {
  heatmap.2(matrix_top30_prot_intensity,
            scale = "column", 
            Colv = T,             
            Rowv = NA,
            dendrogram = "column",
            col = colorRampPalette(colours)(100),
            trace = "none",
            key.title = NA,
            srtCol = 390,
            sepcolor = "black",
            colsep = (0:50),
            sepwidth = .001,
            adjCol = c(.9, .3),
            key.ylab = NA,
            cexCol = 1.2, 
            main = "Z-Scored Heatmap of Top 30 Most Abundant Proteins")
}

# violin plot
violin <- function(colours = colour_scheme) {
  long_imputed_df %>%
    ggplot(aes(y = Intensity, x = Sample, colour = Group)) +
    geom_violin(size = 1, aes(fill = Group), alpha = .3) +
    geom_jitter(alpha = .2) +
    theme_bw() +
    labs(title = "Distribution of Intensity Values",        # ** can change title
         x = "Samples",
         y = "Log2 Intensity") + 
    scale_colour_manual(values = colours) +         # outline colour      # * may change colours
    scale_fill_manual(values = colours) +           # fill colour
    guides(alpha = "none", colour = "none", fill = "none") +
    facet_wrap(~ Group, scales = "free_x")
}

# boxplot
boxplot <- function(colours = colour_scheme) {
  long_imputed_df %>%
    ggplot(aes(y = Intensity, x = Sample, colour = Group)) +
    geom_boxplot(size = 1, aes(fill = Group), alpha = .4) +
    # geom_jitter(alpha = .2) +         # * can add this back in if want the dots all plotted too
    theme_bw() +
    labs(title = "Distribution of Intensity Values",    
         x = "Samples",
         y = "Log2 Intensity") + 
    scale_colour_manual(values = colours) + 
    scale_fill_manual(values = colours) +  
    guides(alpha = "none", colour = "none", fill = "none") +
    facet_wrap(~ Group, scales = "free_x")
}

# Frequency distribution histogram of intensities, per group
histogram <- function(colours = colour_scheme) {
  long_imputed_df %>%
    ggplot(aes(x = Intensity, fill = Group, colour = Group)) +
    geom_histogram(binwidth = .5, alpha = .4, size = 1) +
    facet_wrap("Group") +
    labs(x = "Log2 Intensity", y = "Count", title = "Frequency Distribution Histogram of Intensities") +
    scale_fill_manual(values = colours) +
    scale_colour_manual(values = colours) +
    guides(alpha = "none", colour = "none", fill = "none") +
    theme_bw()
}

# frequency density curves per sample
density_curves_samples <- long_imputed_df %>%
  ggplot(aes(x = Intensity, colour = Sample)) +
  geom_density(size = 1.1) +
  labs(x = "Log2 Intensity", y = "Frequency Density", title = "Frequency Density Curves of Intensities") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.75)) 

# frequency density curves per group
density_curves_groups <- function(colours = colour_scheme) {
  long_imputed_df %>%
    ggplot(aes(x = Intensity, colour = Group)) +
    geom_density(size = 1.1) +
    labs(x = "Log2 Intensity", y = "Frequency Density", title = "Frequency Density Curves of Intensities") +
    scale_colour_manual(values = colours) +
    theme_bw() +
    theme(legend.position = c(0.9, 0.8)) 
}

# histogram for number of peptides per protein (valid proteins only) (across all samples) (cut off at 50)
histogram_pep_per_prot <- function(colour = "grey27") {
  full_norm_prot_dataset %>%
    ggplot(aes(x = Peptide_Count)) +
    geom_histogram(binwidth = 2, alpha = .4, size = 1, colour = colour, fill = colour) +
    labs(x = "Peptide count per protein", y = "Protein count", title = "Histogram of Peptide Count per Protein") +
    guides(alpha = "none", colour = "none", fill = "none") +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1000, by = 10)) +
    theme_bw()
}

# histogram of peps per prot as above, but per experimental group
histogram_pep_per_prot_perG <- function(colours = colour_scheme) {
  peptide_counts_per_group %>%
    ggplot(aes(x = Peptide_Count_byGroup, fill = Group, colour = Group)) +
    geom_histogram(binwidth = 2, alpha = .4, size = 1) +
    facet_wrap("Group") +
    labs(x = "Peptide count per protein", y = "Protein count", title = "Histogram of Peptide Count per Protein") +
    scale_fill_manual(values = colours) +
    scale_colour_manual(values = colours) +
    guides(alpha = "none", colour = "none", fill = "none") +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 4)) +
    scale_y_continuous(breaks = seq(0, 1000, by = 10)) +
    theme_bw()
}





# venn diagrams & data ---------------------------------------------------------------------------------

# function to identify which proteins are in each group, and create a dataset for this
compute_venn_sections <- function(sets) {
  groupsss <- names(sets)  # Get group names
  # List all prots across the groups
  all_elements <- unique(unlist(sets))  # Get all unique prots
  # Identify which groups each prot is present in
  elements_groups <- map(all_elements, function(prot) {
    groups_containing_prot <- groupsss[sapply(sets, function(group_prots) prot %in% group_prots)]
    return(groups_containing_prot)
  })
  # Create a data frame for storing genes and their corresponding group combinations
  venn_data <- tibble(
    prot = all_elements,
    groups = map_chr(elements_groups, ~paste(.x, collapse = ", "))  # Combine group names for each gene
  )
  # Create all combinations of groups dynamically (for any number of groups)
  group_combinations <- unlist(lapply(1:length(groupsss), function(i) combn(groupsss, i, simplify = FALSE)), recursive = FALSE)
  # Label each combination dynamically
  venn_data <- venn_data %>%
    mutate(
      combination_label = map_chr(groups, function(group_str) {
        group_list <- strsplit(group_str, ", ")[[1]]
        # Match the combination in the group_combinations list
        comb_idx <- which(sapply(group_combinations, function(comb) all(sort(group_list) == sort(comb))))
        comb_label <- paste("Present in", paste(group_combinations[[comb_idx]], collapse = " and "))
        return(comb_label)
      })
    ) %>%
    group_by(combination_label) %>%
    summarise(prots = list(prot))  # List of prots for each combination of groups
  return(venn_data)
}

# creating protein data for venn diagram - using protein_dataset (essentially the original input data)
venn_prots <- protein_dataset %>%
  pivot_longer(cols = sample_names,
               names_to = "Sample",
               values_to = "Intensity") %>%
  left_join(metadata, by = "Sample") %>%
  select(Gene, UniProt_AC, Sample, Group, Intensity) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Group) %>%
  summarise(genes = list(unique(Gene))) %>%
  deframe()

# creating peptide data for venn diagram - using peptide_dataset (essentially the original input data)
venn_peps <- peptide_dataset %>%
  pivot_longer(cols = sample_names,
               names_to = "Sample",
               values_to = "Intensity") %>%
  left_join(metadata, by = "Sample") %>%
  filter(!is.na(Intensity)) %>%
  group_by(Group) %>%
  summarise(sequence = list(unique(Peptide_Sequence))) %>%
  deframe()

# using the above function on protein data, to generate lists of proteins present in each group
venn_data <- compute_venn_sections(venn_prots)
# then making it into good format to export
venn_data_alphabetical <- venn_data %>%
  mutate(prots = lapply(prots, function(x) sort(unlist(x))))  # Sort prots alphabetically
# Find the maximum number of genes in any group (for column allocation)
max_prots <- max(sapply(venn_data_alphabetical$prots, length))
# Expand the list column into multiple columns (no names)
venn_data_csv <- venn_data_alphabetical %>%
  mutate(prots = lapply(prots, function(x) c(x, rep("", max_prots - length(x))))) %>%  # Pad with blanks
  unnest_wider(prots, names_sep = "") %>%
  t() %>%
  as.data.frame()
colnames(venn_data_csv) <- NULL

# protein venn
protein_venn <- function(colours = colour_scheme) {
  venn_prots %>%
    ggvenn(fill_color = colours, fill_alpha = .5,
           stroke_size = .5, set_name_size = 5, text_size = 5,
           # show_elements = TRUE, label_sep = "\n"
    ) +
    labs(title = "Proteins Identified in Each Experimental Group") +
    theme(plot.title = element_text(size = 16, hjust = 0.5)) 
}

# peptide venn
peptide_venn <- function(colours = colour_scheme) {
  venn_peps %>%
    ggvenn(fill_color = colours, fill_alpha = .3,
           stroke_size = .5, set_name_size = 5, text_size = 5) +
    labs(title = "Peptides Identified in Each Experimental Group") +
    theme(plot.title = element_text(size = 16, hjust = 0.5))
}




# PCA ---------------------------------------------------------------------------------------------------------------------

imputed_df_labelled <- as.data.frame(t(imputed_df)) %>%
  mutate(Sample = rownames(t(imputed_df))) %>%
  left_join(metadata, by = "Sample")
PCA_res <- prcomp(t(imputed_df), scale.=T)

# Plot PCA graph function:
plot_pca <- function(PCx = 1, PCy = 2, colours = colour_scheme, frame = F, ellipse = F) {
  p <- autoplot(PCA_res, 
                x = PCx, y = PCy, 
                data = imputed_df_labelled, colour = "Group", 
                size = 4, 
                scale = 0, 
                frame = frame) + 
    scale_colour_manual(values = colours) +
    scale_fill_manual(values = colours) +
    labs(title = "PCA of Intensity Values") +
    guides(shape = "none", size = "none", colour = "none", fill = "none") +
    geom_text_repel(aes(label = imputed_df_labelled$Sample, colour = Group), size = 4, box.padding = 0.5, point.padding = 0.9) +
    theme_bw()
  if (ellipse) {
    p <- p + stat_ellipse(aes(group = Group, colour = Group), type = "t", level = 0.95, linetype = "dashed", size = .5, alpha = .5)
  }
  return(p)
}

plot_pca_noTitle <- function(PCx = 1, PCy = 2, colours = colour_scheme, frame = F, ellipse = F) {
  p <- autoplot(PCA_res, 
                x = PCx, y = PCy, 
                data = imputed_df_labelled, colour = "Group", 
                size = 4, 
                scale = 0, 
                frame = frame) + 
    scale_colour_manual(values = colours) +
    scale_fill_manual(values = colours) +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2)) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 2)) +
    guides(shape = "none", size = "none", colour = "none", fill = "none") +
    geom_text_repel(aes(label = imputed_df_labelled$Sample), colour = "grey20", size = 3, box.padding = 0.5, point.padding = 0.9) +
    theme_bw()
  if (ellipse) {
    p <- p + stat_ellipse(aes(group = Group, colour = Group), type = "t", level = 0.95, linetype = "dashed", size = .5, alpha = .5)
  }
  return(p)
}

pca_plots_layout <- function(pca_options) {
  grid.arrange(
    do.call(plot_pca_noTitle, modifyList(pca_options, list(PCx = 1, PCy = 2))), nullGrob(), nullGrob(),
    do.call(plot_pca_noTitle, modifyList(pca_options, list(PCx = 1, PCy = 3))),
    do.call(plot_pca_noTitle, modifyList(pca_options, list(PCx = 2, PCy = 3))), nullGrob(),
    do.call(plot_pca_noTitle, modifyList(pca_options, list(PCx = 1, PCy = 4))),
    do.call(plot_pca_noTitle, modifyList(pca_options, list(PCx = 2, PCy = 4))),
    do.call(plot_pca_noTitle, modifyList(pca_options, list(PCx = 3, PCy = 4))),
    ncol = 3
  )
}




# Limma analysis ---------------------------------------------------------------------------------------------------------------------

groups <- factor(metadata$Group)
# Create the design matrix for your linear model
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)  # Rename columns to match group names
# Fit the linear model
fit <- lmFit(imputed_df, design)
# Generate all pairwise comparisons for every 2-pair of groups
contrast_matrix <- makeContrasts(contrasts = combn(levels(groups), 2, function(x) paste(x, collapse = "-")), 
                                 levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix) # generates pairwise comparisons
fit2 <- eBayes(fit2) # applies empirical bayes moderation (for better stats)

# Extract results and format them
comparison_names <- colnames(contrast_matrix)
# creating results tables for each pairwise comparison:
limma_results_list <- lapply(comparison_names, function(comp_name) {
  # Extract results for the specific comparison
  res <- topTable(fit2, coef = comp_name, adjust = "BH", number = Inf)
  # Add rownames (uniprot_AC) as a new column
  res$UniProt_AC <- rownames(res)
  # Compute -log10(p-values)
  res$NegLog10_P <- -log10(res$P.Value)
  res$NegLog10_AdjP <- -log10(res$adj.P.Val)
  # reverse FC values (so that it reflects the change with the treated condition)
  #res$logFC <- -res$logFC
  # Rename columns to reduce confusion when doing lots of pairs of DE analysis
  colnames(res)[colnames(res) %in% c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "NegLog10_P", "NegLog10_AdjP")] <-
    paste0(comp_name, " ", c("log2FC", "AveExpr", "t", "p-value", "adj p-value", "B", "-log10 p-value", "-log10 adj p-value"))
  return(res)
})

limma_results_list <- lapply(limma_results_list, function(df) {
  left_join(df, gene_AC_key, by = "UniProt_AC")
})

# Assign names to the list corresponding to the comparison names
names(limma_results_list) <- comparison_names

# Merge all results into one large dataframe using full join, and get gene names too
limma_final_results <- Reduce(function(x, y) full_join(x, y, by = c("UniProt_AC", "Gene")), limma_results_list) %>%
  select(UniProt_AC, Gene, everything())

# adding mean intensity value per protein, per experimental group:
# Ensure column names (samples) in imputed_df match metadata$Sample
imputed_df_t <- as.data.frame(t(imputed_df))  # Transpose so samples are rows, genes are columns
imputed_df_t$Sample <- rownames(imputed_df_t)  # Add sample names
# Merge metadata to assign groups to each sample
imputed_df_t <- inner_join(imputed_df_t, metadata, by = c("Sample" = "Sample"))
# Compute mean expression per group
group_means <- imputed_df_t %>%
  select(-Sample) %>%   # Remove Sample column
  group_by(Group) %>%   # Group by experimental condition
  summarise(across(everything(), mean, na.rm = TRUE))  # Compute mean per gene
# Reshape to long format
group_means_long <- melt(group_means, id.vars = "Group")
colnames(group_means_long) <- c("Group", "UniProt_AC", "Mean_Intensity")
# Convert back to wide format
group_means_wide <- dcast(group_means_long, UniProt_AC ~ Group, value.var = "Mean_Intensity")
# Rename columns to follow "group-name mean intensity" format
colnames(group_means_wide)[-1] <- paste0(colnames(group_means_wide)[-1], " mean intensity")
# Merge with limmma_final_results
limma_final_results <- full_join(limma_final_results, group_means_wide, by = "UniProt_AC")

# finally, merge with full imputed protein dataset to get all the info for every protein, alongside limma results
full_imputed_dataset_limma_results <- limma_final_results %>% select(-Gene) %>%
  left_join(full_imputed_dataset, by = "UniProt_AC") %>%
  select(UniProt_AC, Gene, everything())





# Volcano plots ---------------------------------------------------------------------------------------------------------------------

create_volcano_plots <- function(P_val_cutoff = 0.05, fold_change_cutoff = 1.5, sig_colour = "red") {
  volcano_plots_noLabs_list <- list()
  volcano_plots_list <- list()
  
  for (comp_name in comparison_names) {
    # Get the result for the current comparison
    finalresult <- limma_results_list[[comp_name]]
    
    # Dynamically create column names
    logFC_col <- paste0(comp_name, " log2FC")
    adjP_col <- paste0(comp_name, " adj p-value")
    
    # Base volcano plot (without labels)
    volcano_plot_base <- finalresult %>%
      ggplot(aes(x = !!sym(logFC_col), y = -log10(!!sym(adjP_col)))) +
      geom_hline(yintercept = -log10(P_val_cutoff), col = "red2", linetype = 2, size = 0.7) + 
      geom_vline(xintercept = -log2(fold_change_cutoff), col = "red4", linetype = 3, size = 0.7) +
      geom_vline(xintercept = log2(fold_change_cutoff), col = "red4", linetype = 3, size = 0.7) +
      geom_point(size = 1, colour = "grey", alpha = 0.5) +
      geom_point(data = finalresult %>% 
                   filter(get(adjP_col) < P_val_cutoff & abs(get(logFC_col)) > log2(fold_change_cutoff)),
                 size = 1.5, colour = sig_colour, alpha = 0.5) +
      theme_bw() +
      labs(title = paste("Volcano Plot:", comp_name),
           x = "Log2( Fold Change )",
           y = "-Log10( P-adj )")
    
    # Volcano plot with labels
    volcano_plot_with_labels <- volcano_plot_base +
      geom_label_repel(data = finalresult %>%  
                         filter(get(adjP_col) < P_val_cutoff & abs(get(logFC_col)) > log2(fold_change_cutoff)),
                       aes(label = Gene), colour = "black", size = 3,
                       max.overlaps = 20,
                       nudge_y = 0.2,
                       nudge_x = -0.1)
    
    # Store the plots in respective lists
    volcano_plots_noLabs_list[[comp_name]] <- volcano_plot_base
    volcano_plots_list[[comp_name]] <- volcano_plot_with_labels
  }
  
  # Return both lists as a named list
  return(list(volcanos_no_labels = volcano_plots_noLabs_list, volcanos_with_labels = volcano_plots_list))
}






##################################################################################################################
# exporting datasets - run whole thing to export all
##################################################################################################################

# run this to create an output folder (in this directory) - datasets/graphs will be saved here
dir.create(folder_name)

# choose which datasets to export (need to run all lines per section): ------------------------------------------

# full protein and peptide datasets, as imported (before filtering/normalisation/imputation)
# with proper sample names and blank genes removed, and column for peptide count
# and columns for just the first given uniprot_AC and gene
write.csv(protein_dataset, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset.csv")), row.names=FALSE)
write.csv(peptide_dataset, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_dataset.csv")), row.names=FALSE)

# protein and peptide datasets, filtered for valid proteins only (proteins identified by >1 peptide, and with no more than 40% missing values per group)
# proteins in these datasets & onwards contain only valid protein data
write.csv(full_sample_intensities_filtered, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_valid_prots_only.csv")), row.names=FALSE)
write.csv(valid_peptide_dataset, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_dataset_valid_prots_only.csv")), row.names=FALSE)

# peptide counts per valid protein, by experimental group (all the above protein datasets contain the peptide count per protein in general, across groups)
write.csv(peptide_counts_per_group_wGenes, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_counts_per_prot_byGroup.csv")), row.names=FALSE)

# protein dataset as above, but with log2-transformation and normalisation (if applied)
write.csv(full_norm_prot_dataset, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_log2_normalised.csv")), row.names=FALSE)

# imputed dataset (imputed per group, valid proteins only)
write.csv(full_imputed_dataset, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_imputed.csv")), row.names=FALSE)

# imputed dataset, with all differential expression results (every pairwise comparison) - fold change values reflect condition A -> condition B (depends on order given at start)
write.csv(full_imputed_dataset_limma_results, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_DE_results.csv")), row.names=FALSE)

# data to accompany venn diagram - lists of proteins present in each group (from original input data)
write.csv(venn_data_csv, file = file.path(folder_name, paste0(experiment_identifier, "_protein_venn_data.csv")), row.names=FALSE)



##################################################################################################################
# viewing and exporting plots - run whole thing to export all in PDF
##################################################################################################################
# run these lines to view the plots in Rstudio (can save them as an image from here)
# alternatively, can export PDF containing all plots (see bottom)

# plots will use the colour scheme set at the top by default. For that just run: plot()
# alternatively, can set an additional colour scheme here, and run: plot(colour_scheme2)
# remember, the number of colours set must match the number of experimental groups!
colour_scheme2 <- c("purple", "orange")
# to set colours per plot, run: plot(c("red", "yellow"))

# heatmaps colours were set at top, but can create additional colour scheme here too:
heat_colours2 <- c("blue", "red2")

# all graphs are plotted using the final imputed data of valid proteins only (unless stated otherwise)


# missing values plots (no colour option, always grey)
# all proteins identified by >1 peptide (essentially the original input dataset)
vis_miss_allprots
# just valid proteins (before imputation ofc)
vis_miss_validprots

# venn diagrams for the number of proteins/peptides identified in each experimental group
# made with original input protein/peptide datasets
# proteins: (pair with venn data to get gene lists per group)
protein_venn()
# peptides:
peptide_venn()


# standard heatmap of all valid proteins
heatmap_allprots()

# standard heatmap of top 30 valid proteins (by highest median intensity)
heatmap_top30_prots()

# z-scored heatmap of all valid proteins (commented out bc takes ages to load with large datasets)
# zscore_heatmap_allprots()

# z-scored heatmap of top 30 valid proteins (by highest median intensity)
# zscore_heatmap_top30_prots()

# violin plots (distribution of intensity values, per sample, separated into groups)
violin()

# boxplots (distribution of intensity values, per sample, separated into groups)
boxplot()

# frequency distribution histogram of intensity values, per group
histogram()

# density curves for samples (no colour option - every sample is a different colour of the rainbow :))
density_curves_samples

# density curves for experimental groups
density_curves_groups()

# histogram of the number of peptides per protein (just one colour for this one) - across all samples, and x-axis is cut off at 50
histogram_pep_per_prot()

# as above, but per experimental group
histogram_pep_per_prot_perG()


# options for the PCA plot - can add frames around points (T/F) and 95% confidence ellipses (T/F)
# can replace colour_scheme with colour_scheme2, or other colours of choice
pca_options <- list(frame = T, ellipse = F, colours = colour_scheme)

# Main PCA plot (PC1 v PC2)
do.call(plot_pca, pca_options)
# grid of all PCA plots, combinations up to PC4
pca_plots_layout(pca_options)


# volcano plots - specify p-value and fold change cut-off values here, then run this line (e.g. FC = 1.5 is a 50% increase/decrease)
# can also change the colour of significant points (non-sig points are grey)
all_volcanos <- create_volcano_plots(P_val_cutoff = 0.05, fold_change_cutoff = 1.5, sig_colour = "red")
# view volcano plots, no point labels:  (shows plot made per pairwise comparison)
do.call(grid.arrange, c(all_volcanos$volcanos_no_labels, ncol = 2))
# with point labels: (note if there are too many points to label, some will be missed)
do.call(grid.arrange, c(all_volcanos$volcanos_with_labels, ncol = 2))



# exporting all plots in a PDF: -------------------------------------------------------------------------------------
# again, can change colours if wanted
# must run all lines (though can remove a plot if u want)
# will use the above set options for PCA and volcano plots

pdf(file = file.path(folder_name, paste0(experiment_identifier, "_all_plots.pdf")), width = 10, height = 8)  # Adjust the size as needed

vis_miss_allprots
vis_miss_validprots
protein_venn()
peptide_venn()
heatmap_allprots()
heatmap_top30_prots()
# zscore_heatmap_allprots()
zscore_heatmap_top30_prots()
violin()
boxplot()
histogram()
density_curves_samples
density_curves_groups()
histogram_pep_per_prot()
histogram_pep_per_prot_perG()
do.call(plot_pca, pca_options)
pca_plots_layout(pca_options)
for (volcano in all_volcanos$volcanos_no_labels) {
  print(volcano)
}
for (volcano in all_volcanos$volcanos_with_labels) {
  print(volcano)
}

dev.off()

# end --------------------------------------






