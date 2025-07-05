
library(tidyverse)
library(data.table)


# function:
calculate_segment_coverage <- function(df) {
  df_name <- deparse(substitute(df))
  df_list <- split(df, df$GeneAC)
  result_df <- data.frame(GeneAC = character(),
                          num_segments = integer(),
                          num_segments_present = integer(),
                          proportion_covered = numeric(),
                          stringsAsFactors = FALSE)
  for (gene in names(df_list)) {
    current_df <- df_list[[gene]]
    num_segments <- nrow(current_df)
    num_segments_present <- sum(current_df[[condition1]] != 0)
    if (num_segments_present == 0) {
      num_segments_present <- sum(current_df[[condition2]] != 0)
    }
    proportion_covered <- ifelse(num_segments > 0, num_segments_present / num_segments, NA)
    result_df <- rbind(result_df, data.frame(GeneAC = gene,
                                             num_segments = num_segments,
                                             num_segments_present = num_segments_present,
                                             proportion_covered = proportion_covered))
  }
  # Dynamically assign the result dataframe with the appropriate name
  result_name <- paste0(df_name, "_PLF_segment_coverage")
  assign(result_name, result_df, envir = .GlobalEnv)
  write.csv(result_df, file = paste0(df_name, "_PLF_segment_coverage.csv"), row.names = F)
  return(result_df)
}



#############################################################################################################


# read in PLF raw data (local output tsv, or raw data tsv downloaded off MPLF website)

aim2_plf_all <- fread("aim2_allprots_PLFresults.tsv")
aim2_plf_BM <- fread("aim2_BMprots_PLFresults.tsv")




# write your experimental conditions (as they appear in data / experiment feed)
condition1 = "Epidermis+"
condition2 = "Epidermis-"



# run segment coverage function, with raw data (will export the csv)
calculate_segment_coverage(aim2_plf_all)









# calculating segment coverage per condition:
calculate_segment_coverage_perCondition <- function(df) {
  df_name <- deparse(substitute(df))
  df_list <- split(df, df$GeneAC)
  result_df <- data.frame(GeneAC = character(),
                          num_segments = integer(),
                          num_segments_present_epiY = integer(),
                          num_segments_present_epiN = integer(),
                          proportion_covered_epiY = numeric(),
                          proportion_covered_epiN = numeric(),
                          stringsAsFactors = FALSE)
  for (gene in names(df_list)) {
    current_df <- df_list[[gene]]
    num_segments <- nrow(current_df)
    num_segments_present_epiY <- sum(current_df[[condition1]] != 0)
    proportion_covered_epiY <- ifelse(num_segments > 0, num_segments_present_epiY / num_segments, NA)
    num_segments_present_epiN <- sum(current_df[[condition2]] != 0)
    proportion_covered_epiN <- ifelse(num_segments > 0, num_segments_present_epiN / num_segments, NA)
    result_df <- rbind(result_df, data.frame(GeneAC = gene,
                                             num_segments = num_segments,
                                             num_segments_present_epiY = num_segments_present_epiY,
                                             proportion_covered_epiY = proportion_covered_epiY,
                                             num_segments_present_epiN = num_segments_present_epiN,
                                             proportion_covered_epiN = proportion_covered_epiN))
  }
  # Dynamically assign the result dataframe with the appropriate name
  result_name <- paste0(df_name, "_PLF_segment_coverage")
  assign(result_name, result_df, envir = .GlobalEnv)
  write.csv(result_df, file = paste0(df_name, "_PLF_segment_coverage_per_condition.csv"), row.names = F)
  return(result_df)
}


calculate_segment_coverage_perCondition(aim2_plf_BM)

calculate_segment_coverage_perCondition(aim2_plf_all)







################################################


# function to rename cols, so can put all in same dataset for speclib comparison
rename_cols <- function(dfx) {
  dataset_name <- deparse(substitute(dfx))
  dfx %>%
    rename(!!paste0(dataset_name, "_num_segments") := `num_segments`) %>%
    rename(!!paste0(dataset_name, "_num_segments_present") := `num_segments_present`) %>%
    rename(!!paste0(dataset_name, "_proportion_covered") := `proportion_covered`) %>%
    return(dfx)
}


DIANN_noMods_allprots <- rename_cols(DIANN_noMods_allprots)
DIANN_Mods_noTuning_allprots <- rename_cols(DIANN_Mods_noTuning_allprots)
DIANN_G_tune_allprots <- rename_cols(DIANN_G_tune_allprots)
FRAGPIPE_G_allprots <- rename_cols(FRAGPIPE_G_allprots)

DIANN_noMods_ECMprots <- rename_cols(DIANN_noMods_ECMprots)
DIANN_Mods_noTuning_ECMprots <- rename_cols(DIANN_Mods_noTuning_ECMprots)
DIANN_G_tune_ECMprots <- rename_cols(DIANN_G_tune_ECMprots)
FRAGPIPE_G_ECMprots <- rename_cols(FRAGPIPE_G_ECMprots)



DIANN_noMods_allprots
DIANN_Mods_noTuning_allprots
DIANN_G_tune_allprots
FRAGPIPE_G_allprots

DIANN_noMods_BMprots
DIANN_Mods_noTuning_BMprots
DIANN_G_tune_BMprots
FRAGPIPE_G_BMprots




all_prots_all_libs_plf_segement_coverage <- DIANN_noMods_allprots %>%
  full_join(DIANN_Mods_noTuning_allprots, by = "GeneAC") %>%
  full_join(DIANN_G_tune_allprots, by = "GeneAC") %>%
  full_join(FRAGPIPE_G_allprots, by = "GeneAC")

ECM_prots_all_libs_plf_segement_coverage <- DIANN_noMods_ECMprots %>%
  full_join(DIANN_Mods_noTuning_ECMprots, by = "GeneAC") %>%
  full_join(DIANN_G_tune_ECMprots, by = "GeneAC") %>%
  full_join(FRAGPIPE_G_ECMprots, by = "GeneAC")





# write.csv(aim3_plf_BM_Segemnt, file = "aim3_BMprots_plf_segment_coverage.csv", row.names = F)









