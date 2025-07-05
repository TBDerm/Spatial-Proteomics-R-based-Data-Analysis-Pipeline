# input protein & peptide datasets (any format works - just need uniprot AC, gene, peptide sequence, and sample data columns)
# script removes blank genes, renames samples, and extracts the first given uniprot AC and Gene (if multiple in the same column)
# creates datasets for all proteins, ECM proteins only, and BM proteins only
# creates datasets for all proteins with ECM / BM proteins annotated (useful for graph labelling etc)
# creates PLF-ready input data for each protein subset                      (*disabled* and input data for the coverage summariser)

library(tidyverse)
library(data.table)
library(dplyr)


##################################################################################################################
# need to edit this part to suit your data
##################################################################################################################

# # reading intensity data (DIA-NN output or similar)
# protein_dataset <- fread('Aim3_photoage_N3_protein_dataset_all+ECM+BM_noted.csv') # proteins
# peptide_dataset <- fread('Aim3_photoage_N3_peptide_dataset_all+ECM+BM_noted.csv') # peptides
# 
# # edit sample names and experimental groups (must be in correct order, must match) (each sample name must be unique)
# # number of sample and group names given must match the number of sample columns in data
# # (these will be used for graphs so ensure they're spelt right etc)
# sample_names <- c("0841_Buttock", "0842_Buttock", "0863_Buttock",
#                   "0841_Forearm", "0842_Forearm", "0863_Forearm")
# # format should be rep("group name", x) - where x is the number of samples (columns) in that experimental group
# group_names <- c(rep("Buttock", 3),
#                  rep("Forearm", 3))


####################################################

protein_dataset <- fread('test_protein_dataset_all.csv') # proteins
peptide_dataset <- fread('test_peptide_dataset_all.csv') # peptides

# edit sample names and experimental groups (must be in correct order, must match) (each sample name must be unique)
# number of sample and group names given must match the number of sample columns in data
# (these will be used for graphs so ensure they're spelt right etc)

sample_names <- c("S1", "S2", "S3", "S4", "S5", "S6",
                  "S7", "S8", "S9", "S10", "S11", "S12")
# format should be rep("group name", x) - where x is the number of samples (columns) in that experimental group

group_names <- c(rep("A", 1),
                 rep("B", 1),
                 rep("C", 1),
                 rep("D", 1),
                 rep("E", 1),
                 rep("F", 1),
                 rep("G", 1),
                 rep("H", 1),
                 rep("I", 1),
                 rep("J", 1),
                 rep("K", 1),
                 rep("L", 1))



##################################################



# enter the number range for the sample columns in protein and peptide data
# e.g. if columns 5, 6, 7, 8, 9, and 10 are my samples (6 in total), enter 5:10 here.
prot_sample_cols <- 5:16  # in protein data
pep_sample_cols <- 11:22  # in peptide data

# enter the number for the column containing uniprot AC numbers (in DIANN data, this is the Protein.Group column)
prot_uniprot_AC_col <- 1  # in protein data
pep_uniprot_AC_col <- 1  # in peptide data

# enter the number for the column containing Genes
prot_gene_col <- 3  # in protein data
pep_gene_col <- 4  # in peptide data

# enter the number for the column containing peptide sequence (stripped.sequence in DIANN)
pep_sequence_col <- 7 # in peptide data

# the script will make a new folder within current directory to store output data - provide a name for the folder
folder_name <- "HDF_test_datasets"

# output files will be named with this at the beginning - provide a short experiment ID
experiment_identifier <- "HDF_test"

##################################################################################################################
# the rest can just be run (all datasets will be exported to a folder)
# or can manually export chosen datasets in last section
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

# remove rows that do not have a Gene (these are non-human contaminants)
# and taking only the first Uniprot_AC from the Uniprot_AC column (generally when more than one is given, the first is the main one)
# same for genes - taking only the first
# (stores the original Uniprot_AC and Genes columns, with multiple codes, in copy columns called 'multiple_ACs' and 'multiple_genes')
protein_dataset <- protein_dataset[!(protein_dataset$Gene == "" | is.na(protein_dataset$Gene)), ] %>%
  mutate("multiple_ACs" = UniProt_AC) %>%
  mutate("multiple_genes" = Gene) %>%
  separate(UniProt_AC, into = c("UniProt_AC", "other"), sep = ";") %>%
  separate(Gene, into = c("Gene", "otherG"), sep = ";") %>%
  select(-other, -otherG)
peptide_dataset <- peptide_dataset[!(peptide_dataset$Gene == "" | is.na(peptide_dataset$Gene)), ] %>%
  mutate("multiple_ACs" = UniProt_AC) %>%
  mutate("multiple_genes" = Gene) %>%
  separate(UniProt_AC, into = c("UniProt_AC", "other"), sep = ";") %>%
  separate(Gene, into = c("Gene", "otherG"), sep = ";") %>%
  select(-other, -otherG)

# rename sample columns
names(protein_dataset)[prot_sample_cols] <- sample_names
names(peptide_dataset)[pep_sample_cols] <- sample_names

# ECM protein-only datasets
ECM_matrisome <- read.csv("ECM_matrisome_expanded.csv") %>% select(-Gene)
protein_dataset_ECM <- inner_join(ECM_matrisome, protein_dataset, by = "UniProt_AC")
peptide_dataset_ECM <- inner_join(ECM_matrisome, peptide_dataset, by = "UniProt_AC")

# BM protein-only datasets
skin_BM_proteins <- read.csv("skin_BM_proteins.csv") %>% select(-Gene, -Uniprot_ID)
protein_dataset_BM <- inner_join(skin_BM_proteins, protein_dataset, by = "UniProt_AC")
peptide_dataset_BM <- inner_join(skin_BM_proteins, peptide_dataset, by = "UniProt_AC")

# all proteins, but ECM proteins are annotated
protein_dataset_ECM_noted <- left_join(protein_dataset, ECM_matrisome, by = "UniProt_AC")
peptide_dataset_ECM_noted <- left_join(peptide_dataset, ECM_matrisome, by = "UniProt_AC")

# all proteins, but BM proteins are annotated
protein_dataset_BM_noted <- left_join(protein_dataset, skin_BM_proteins, by = "UniProt_AC")
peptide_dataset_BM_noted <- left_join(peptide_dataset, skin_BM_proteins, by = "UniProt_AC")

# all proteins, but ECM AND BM proteins are annotated
ECM_BM_both_proteomes <- left_join(ECM_matrisome, skin_BM_proteins, by = "UniProt_AC")
protein_dataset_all_ECM_BM_noted <- left_join(protein_dataset, ECM_BM_both_proteomes, by = "UniProt_AC")
peptide_dataset_all_ECM_BM_noted <- left_join(peptide_dataset, ECM_BM_both_proteomes, by = "UniProt_AC")

# ECM proteins-only, but BM proteins are annotated
protein_dataset_ECM_BM_noted <- left_join(protein_dataset_ECM, skin_BM_proteins, by = "UniProt_AC")
peptide_dataset_ECM_BM_noted <- left_join(peptide_dataset_ECM, skin_BM_proteins, by = "UniProt_AC")

# PLF-ready data
PLF_input_allprots <- as.data.frame(peptide_dataset) %>%
  filter(Proteotypic == 1) %>%
  select(all_of("Protein.Names"), all_of("Peptide_Sequence"), all_of(c(sample_names))) %>%
  dplyr::rename(Protein = Protein.Names, Peptide = Peptide_Sequence) %>%
  pivot_longer(cols = sample_names,
               names_to = "Sample",
               values_to = "spectra") %>%
  na.omit(spectra)

PLF_input_ECMprots <- as.data.frame(peptide_dataset_ECM) %>%
  filter(Proteotypic == 1) %>%
  select(all_of("Protein.Names"), all_of("Peptide_Sequence"), all_of(c(sample_names))) %>%
  dplyr::rename(Protein = Protein.Names, Peptide = Peptide_Sequence) %>%
  pivot_longer(cols = sample_names,
               names_to = "Sample",
               values_to = "spectra") %>%
  na.omit(spectra)

PLF_input_BMprots <- as.data.frame(peptide_dataset_BM) %>%
  filter(Proteotypic == 1) %>%
  select(all_of("Protein.Names"), all_of("Peptide_Sequence"), all_of(c(sample_names))) %>%
  dplyr::rename(Protein = Protein.Names, Peptide = Peptide_Sequence) %>%
  pivot_longer(cols = sample_names,
               names_to = "Sample",
               values_to = "spectra") %>%
  na.omit(spectra)

PLF_experiment_feed <- metadata %>%
  group_by(Group) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = Group, values_from = Sample) %>%
  select(-row)

# # input data for the protein coverage summarizer, across all samples/groups
# allprots_coverage_input <- peptide_dataset %>%
#   select(UniProt_AC, Peptide_Sequence)
# ECM_prots_coverage_input <- peptide_dataset_ECM %>%
#   select(UniProt_AC, Peptide_Sequence)
# BM_prots_coverage_input <- peptide_dataset_BM %>%
#   select(UniProt_AC, Peptide_Sequence)
# 
# # input data for protein coverage summarizer, per experimental group:
# unique_groups <- unique(metadata$Group)
# # all proteins
# allprots_coverage_input_list <- list()
# for (group in unique_groups) {
#   sample_names <- metadata$Sample[metadata$Group == group]
#   coverage_input <- peptide_dataset %>%         
#     filter(rowSums(is.na(select(., all_of(sample_names)))) != length(sample_names)) %>%
#     select(UniProt_AC, Peptide_Sequence)
#   allprots_coverage_input_list[[group]] <- coverage_input
# }
# # ECM proteins
# ECMprots_coverage_input_list <- list()
# for (group in unique_groups) {
#   sample_names <- metadata$Sample[metadata$Group == group]
#   coverage_input <- peptide_dataset_ECM %>%
#     filter(rowSums(is.na(select(., all_of(sample_names)))) != length(sample_names)) %>%
#     select(UniProt_AC, Peptide_Sequence)
#   ECMprots_coverage_input_list[[group]] <- coverage_input
# }
# # BM proteins
# BMprots_coverage_input_list <- list()
# for (group in unique_groups) {
#   sample_names <- metadata$Sample[metadata$Group == group]
#   coverage_input <- peptide_dataset_BM %>%   
#     filter(rowSums(is.na(select(., all_of(sample_names)))) != length(sample_names)) %>%
#     select(UniProt_AC, Peptide_Sequence)
#   BMprots_coverage_input_list[[group]] <- coverage_input
# }




##################################################################################################################
# exporting datasets
##################################################################################################################

# run this to create an output folder (in this directory) - datasets will be saved here
dir.create(folder_name)

# choose which datasets to export (run the line)

# full protein and peptide datasets, just with proper sample names and blank genes removed
# and just the first given uniprot_AC and gene are used as row identifiers
write.csv(protein_dataset, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_all.csv")), row.names=FALSE)
write.csv(peptide_dataset, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_dataset_all.csv")), row.names=FALSE)

# as above, ECM proteins only
write.csv(protein_dataset_ECM, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_ECM.csv")), row.names=FALSE)
write.csv(peptide_dataset_ECM, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_dataset_ECM.csv")), row.names=FALSE)

# as above, BM proteins only
write.csv(protein_dataset_BM, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_BM.csv")), row.names=FALSE)
write.csv(peptide_dataset_BM, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_dataset_BM.csv")), row.names=FALSE)

# full protein and peptide datasets, with ECM proteins annotated
write.csv(protein_dataset_ECM_noted, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_all+ECM_noted.csv")), row.names=FALSE)
write.csv(peptide_dataset_ECM_noted, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_dataset_all+ECM_noted.csv")), row.names=FALSE)

# full protein and peptide datasets, with BM proteins annotated
write.csv(protein_dataset_BM_noted, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_all+BM_noted.csv")), row.names=FALSE)
write.csv(peptide_dataset_BM_noted, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_dataset_all+BM_noted.csv")), row.names=FALSE)

# full protein and peptide datasets, with ECM AND BM proteins annotated
write.csv(protein_dataset_all_ECM_BM_noted, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_all+ECM+BM_noted.csv")), row.names=FALSE)
write.csv(peptide_dataset_all_ECM_BM_noted, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_dataset_all+ECM+BM_noted.csv")), row.names=FALSE)

# ECM-only protein and peptide datasets, with BM proteins annotated
write.csv(protein_dataset_ECM_BM_noted, file = file.path(folder_name, paste0(experiment_identifier, "_protein_dataset_ECM+BM_noted.csv")), row.names=FALSE)
write.csv(peptide_dataset_ECM_BM_noted, file = file.path(folder_name, paste0(experiment_identifier, "_peptide_dataset_ECM+BM_noted.csv")), row.names=FALSE)

# PLF input data
write.table(PLF_experiment_feed, file = file.path(folder_name, paste0(experiment_identifier, "_PLF_experiment_feed.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(PLF_input_allprots, file = file.path(folder_name, paste0(experiment_identifier, "_PLF_input_allprots.csv")), row.names=FALSE)
write.csv(PLF_input_BMprots, file = file.path(folder_name, paste0(experiment_identifier, "_PLF_input_BMprots_only.csv")), row.names=FALSE)
write.csv(PLF_input_ECMprots, file = file.path(folder_name, paste0(experiment_identifier, "_PLF_input_ECMprots_only.csv")), row.names=FALSE)

# # coverage summarizer input - across all samples and groups:
# write.table(allprots_coverage_input, file = file.path(folder_name, paste0(experiment_identifier, "_coverage_input_allprots.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(ECM_prots_coverage_input, file = file.path(folder_name, paste0(experiment_identifier, "_coverage_input_ECMprots.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(BM_prots_coverage_input, file = file.path(folder_name, paste0(experiment_identifier, "_coverage_input_BMprots.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
# 
# # coverage summarizer input data - per experimental group:
# # all proteins:
# for (group in names(allprots_coverage_input_list)) {
#   file_name <- file.path(folder_name, paste0(experiment_identifier, "_coverage_input_allprots_", group, ".txt"))
#   write.table(allprots_coverage_input_list[[group]], file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
# }
# # ECM proteins:
# for (group in names(ECMprots_coverage_input_list)) {
#   file_name <- file.path(folder_name, paste0(experiment_identifier, "_coverage_input_ECMprots_", group, ".txt"))
#   write.table(ECMprots_coverage_input_list[[group]], file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
# }
# # BM proteins:
# for (group in names(BMprots_coverage_input_list)) {
#   file_name <- file.path(folder_name, paste0(experiment_identifier, "_coverage_input_BMprots_", group, ".txt"))
#   write.table(BMprots_coverage_input_list[[group]], file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
# }



# end -----------------------------------------------------------------------------------------------------------



