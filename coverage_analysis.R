# want to get rid of rows with 0 in coverage % column
# make table or something with the common proteins from two datasets, and a column for each dataset's coverage %
# make double-bar chart - bar per protein, so probs do w just ECM prots

# to tidy and merge multiple protein coverage summariser output datasets


library(tidyverse)
library(data.table)
library(dplyr)


# import data (straight from the coverage summariser program)
# include the speclib details in the variable name (will be used to name columns in final dataset)

G_frag_ECM <- fread("pap-ret_Fragpipe_G_coverage_input_ECMprots_coverage.txt")
G_frag_all <- fread("pap-ret_Fragpipe_G_coverage_input_allprots_coverage.txt")

diann__noMods_ECM <- fread("pap-ret_diann_noMods_coverage_input_ECMprots_coverage.txt")
diann__noMods_all <- fread("pap-ret_diann_noMods_coverage_input_allprots_coverage.txt")

diann__wMods_noTuning_ECM <- fread("pap-ret_diann_wMods_noTuning_coverage_input_ECMprots_coverage.txt")
diann__wMods_noTuning_all <- fread("pap-ret_diann_wMods_noTuning_coverage_input_allprots_coverage.txt")

diann__wMods_Gtuned_ECM <- fread("pap-ret_diann_wMods_Gtuned_coverage_input_ECMprots_coverage.txt")
diann__wMods_Gtuned_all <- fread("pap-ret_diann_wMods_Gtuned_coverage_input_ECMprots_coverage.txt")



# function to get the gene per row, 
# get gene, percent coverage, and peptide count cols
# and rename these cols so they're identifiable in big data
tidier <- function(dfx) {
  dataset_name <- deparse(substitute(dfx))
  dfx %>%
  separate(`Protein Description`, into = c("Protein Description", "Gene"), sep = "GN=") %>%
  separate(Gene, into = c("Gene", "other"), sep = " PE") %>%
  select(Gene, `Percent Coverage`, `Unique Peptide Count`) %>%
  filter(`Percent Coverage` != 0) %>%
  rename(!!paste0(dataset_name, "_percent_coverage") := `Percent Coverage`) %>%
  rename(!!paste0(dataset_name, "_unique_peptide_count") := `Unique Peptide Count`)
}


# use tidier() function on all the imported data
G_frag_ECM <- tidier(G_frag_ECM)
G_frag_all <- tidier(G_frag_all)

diann__noMods_ECM <- tidier(diann__noMods_ECM)
diann__noMods_all <- tidier(diann__noMods_all)

diann__wMods_noTuning_ECM <- tidier(diann__wMods_noTuning_ECM)
diann__wMods_noTuning_all <- tidier(diann__wMods_noTuning_all)

diann__wMods_Gtuned_ECM <- tidier(diann__wMods_Gtuned_ECM)
diann__wMods_Gtuned_all <- tidier(diann__wMods_Gtuned_all)



# use full_join(x, y, by = "Gene") to join the datasets
# just do it again to add a third dataset

all_libs_allProts <- G_frag_all %>%
  full_join(diann__noMods_all, by = "Gene") %>%
  full_join(diann__wMods_noTuning_all, by = "Gene") %>%
  full_join(diann__wMods_Gtuned_all, by = "Gene")


all_libs_ECMProts <- G_frag_ECM %>%
  full_join(diann__noMods_ECM, by = "Gene") %>%
  full_join(diann__wMods_noTuning_ECM, by = "Gene") %>%
  full_join(diann__wMods_Gtuned_ECM, by = "Gene")


### merge w matrisome to get categories
matrisome <- fread("ECM_matrisome_expanded.csv")
all_libs_ECMProts_matrisome <- left_join(all_libs_ECMProts, matrisome, by = "Gene") %>%
  distinct(Gene, .keep_all = TRUE)


# export datasets

write.csv(all_libs_allProts, file = "all_libs_allProts_coverage.csv", row.names=FALSE)
write.csv(all_libs_ECMProts_matrisome, file = "all_libs_ECMProts_coverage.csv", row.names=FALSE)










