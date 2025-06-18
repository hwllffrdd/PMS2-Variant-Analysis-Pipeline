#!/usr/bin/env Rscript

# Enhanced compile_results.R
# This script provides comprehensive analysis of PMS2 variant results,
# including detailed explanations for LR-PCR recommendations

library(dplyr)
library(readr)
library(stringr)
library(knitr)
library(purrr)

# Define the results directory - modify as needed
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  results_dir <- args[1]
} else {
  # Default to the current working directory
  results_dir <- "results/PMS2_analysis_2025-03-12-run116"
}

# Path to final results
final_results_dir <- file.path(results_dir, "final_results")

# Get all variant files
variant_files <- list.files(
  path = final_results_dir,
  pattern = "_PMS2_variants.txt$",
  full.names = TRUE
)

cat(paste0("Found ", length(variant_files), " variant files in ", final_results_dir, "\n"))

# Function to analyze a single variant file
analyze_variant_file <- function(file_path) {
  sample_name <- basename(file_path) %>% str_replace("_PMS2_variants.txt", "")
  cat(paste0("\nProcessing sample: ", sample_name, "\n"))
  
  # Read the variant file
  variants <- try(read_delim(file_path, delim = "\t"), silent = TRUE)
  
  if (inherits(variants, "try-error")) {
    return(list(
      sample_name = sample_name,
      error = TRUE,
      message = "Error reading file",
      variants = NULL
    ))
  }
  
  # Basic summary
  summary <- list(
    sample_name = sample_name,
    error = FALSE,
    total_variants = nrow(variants),
    in_roi = sum(variants$INROI == TRUE, na.rm = TRUE),
    pass_filter = sum(variants$FILTER.x == "PASS", na.rm = TRUE),
    pathogenic = sum(variants$class %in% c("PAT", "lPAT"), na.rm = TRUE),
    potentially_pathogenic = sum(grepl("splice|nonsense|frameshift", variants$cdna, ignore.case = TRUE), na.rm = TRUE),
    need_lrpcr = sum(variants$what_to_do == "Perform LR-PCR", na.rm = TRUE)
  )
  
  # Extract pathogenic variants
  pathogenic_variants <- variants %>%
    filter(class %in% c("PAT", "lPAT") | 
             grepl("splice|nonsense|frameshift", cdna, ignore.case = TRUE)) %>%
    select(ID, cdna, prot, class, class_paralogous, paralogous_above60, AF, FILTER.x, INROI, present_pipelines, what_to_do)
  
  # Decision analysis for potential pathogenic variants
  if (nrow(pathogenic_variants) > 0) {
    pathogenic_variants <- pathogenic_variants %>%
      mutate(
        reason_for_decision = case_when(
          what_to_do == "Perform LR-PCR" ~ "Meets all criteria for LR-PCR",
          what_to_do == "Only perform LR-PCR if IHC PMS2-" ~ "Only detected in one pipeline approach",
          what_to_do == "QUALITY FILTER NOT PASSED" ~ paste0("Failed quality filter: ", FILTER.x),
          what_to_do == "Out of ROI" ~ "Variant outside region of interest",
          what_to_do == "Do not classify, paralogous variant <60" ~ paste0("Paralogous variant with AF=", AF),
          TRUE ~ what_to_do
        )
      )
  }
  
  # Variants with high potential but weren't flagged for LR-PCR
  missed_potential <- variants %>%
    filter(
      (class %in% c("PAT", "lPAT") | 
         grepl("splice|nonsense|frameshift", cdna, ignore.case = TRUE)) &
        what_to_do != "Perform LR-PCR"
    ) %>%
    select(ID, cdna, prot, class, class_paralogous, paralogous_above60, AF, FILTER.x, INROI, present_pipelines, what_to_do)
  
  # Return all the analysis results
  return(list(
    sample_name = sample_name,
    error = FALSE,
    summary = summary,
    pathogenic_variants = pathogenic_variants,
    missed_potential = missed_potential,
    variants = variants
  ))
}

# Process all files
results <- map(variant_files, analyze_variant_file)
names(results) <- map_chr(results, ~.$sample_name)

# Create summary table
summary_table <- map_dfr(results, function(r) {
  if (r$error) {
    return(data.frame(
      sample_name = r$sample_name,
      total_variants = NA,
      in_roi = NA,
      pass_filter = NA,
      pathogenic = NA,
      potentially_pathogenic = NA,
      need_lrpcr = NA,
      num_potential_missed = NA
    ))
  }
  
  data.frame(
    sample_name = r$sample_name,
    total_variants = r$summary$total_variants,
    in_roi = r$summary$in_roi,
    pass_filter = r$summary$pass_filter,
    pathogenic = r$summary$pathogenic,
    potentially_pathogenic = r$summary$potentially_pathogenic,
    need_lrpcr = r$summary$need_lrpcr,
    num_potential_missed = nrow(r$missed_potential)
  )
})

# Print overall summary
cat("\n===== OVERALL SUMMARY =====\n")
cat(paste0("Total samples processed: ", nrow(summary_table), "\n"))
cat(paste0("Total variants found: ", sum(summary_table$total_variants, na.rm = TRUE), "\n"))
cat(paste0("Variants requiring LR-PCR: ", sum(summary_table$need_lrpcr, na.rm = TRUE), "\n"))
cat(paste0("Variants with pathogenic classification: ", sum(summary_table$pathogenic, na.rm = TRUE), "\n"))
cat(paste0("Variants with potentially pathogenic consequences: ", sum(summary_table$potentially_pathogenic, na.rm = TRUE), "\n"))
cat(paste0("Potentially pathogenic variants not flagged for LR-PCR: ", sum(summary_table$num_potential_missed, na.rm = TRUE), "\n"))

# Print sample summary table
cat("\n===== SAMPLE SUMMARY =====\n")
print(summary_table)

# Create detailed report for each sample
cat("\n===== DETAILED SAMPLE REPORTS =====\n")
for (r in results) {
  if (r$error) {
    cat(paste0("\nSample: ", r$sample_name, " - ERROR: ", r$message, "\n"))
    next
  }
  
  cat(paste0("\n\nSample: ", r$sample_name, "\n"))
  cat("-------------------------------\n")
  
  # Print basic metrics
  cat(paste0("Total variants: ", r$summary$total_variants, "\n"))
  cat(paste0("Variants in ROI: ", r$summary$in_roi, "\n"))
  cat(paste0("PASS filter variants: ", r$summary$pass_filter, "\n"))
  cat(paste0("Classified pathogenic variants: ", r$summary$pathogenic, "\n"))
  cat(paste0("Variants with pathogenic consequences: ", r$summary$potentially_pathogenic, "\n"))
  cat(paste0("Variants requiring LR-PCR: ", r$summary$need_lrpcr, "\n"))
  
  # Print pathogenic variants
  if (nrow(r$pathogenic_variants) > 0) {
    cat("\n== Potentially Pathogenic Variants ==\n")
    print(r$pathogenic_variants, n = 100)
  } else {
    cat("\n== No Potentially Pathogenic Variants Found ==\n")
  }
  
  # Print missed potential variants
  if (nrow(r$missed_potential) > 0) {
    cat("\n== Potentially Pathogenic Variants NOT Flagged for LR-PCR ==\n")
    print(r$missed_potential, n = 100)
  }
  
  # Print variant counts by what_to_do
  if (nrow(r$variants) > 0) {
    cat("\n== Variant Decision Distribution ==\n")
    decision_counts <- table(r$variants$what_to_do)
    print(decision_counts)
  }
}

# Create comprehensive summary CSV
all_pathogenic_variants <- bind_rows(lapply(results, function(r) {
  if (r$error || nrow(r$pathogenic_variants) == 0) return(NULL)
  
  r$pathogenic_variants %>%
    mutate(sample_name = r$sample_name)
}))

all_missed_potential <- bind_rows(lapply(results, function(r) {
  if (r$error || nrow(r$missed_potential) == 0) return(NULL)
  
  r$missed_potential %>%
    mutate(sample_name = r$sample_name)
}))

# Write summary tables to CSV
write_csv(summary_table, file.path(results_dir, "enhanced_variant_summary.csv"))

if (nrow(all_pathogenic_variants) > 0) {
  write_csv(all_pathogenic_variants, file.path(results_dir, "pathogenic_variants.csv"))
}

if (nrow(all_missed_potential) > 0) {
  write_csv(all_missed_potential, file.path(results_dir, "missed_pathogenic_variants.csv"))
}

cat("\nAnalysis complete. Summary files written to:\n")
cat(paste0("  ", file.path(results_dir, "enhanced_variant_summary.csv"), "\n"))
if (nrow(all_pathogenic_variants) > 0) {
  cat(paste0("  ", file.path(results_dir, "pathogenic_variants.csv"), "\n"))
}
if (nrow(all_missed_potential) > 0) {
  cat(paste0("  ", file.path(results_dir, "missed_pathogenic_variants.csv"), "\n"))
}