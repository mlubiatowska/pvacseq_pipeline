library(tidyverse)

# Loop through each sample directory - function to extract all result for combined class I and class II neoag)
merging_pvacseq_combined_results <- function(list_of_paths) {
  all_samples <- list()
  
  for (sample_path in list_of_paths) {
    sample_name <- basename(sample_path)
    combined_dir <- file.path(sample_path, 'combined')
    
    # Find the filtered TSV file
    tsv_files <- list.files(combined_dir, pattern = '\\.Combined\\.filtered\\.tsv$', full.names = TRUE)
    
    if (length(tsv_files) > 0) {
      # Read the TSV file with consistent column types
      data <- read_delim(file = tsv_files[1], delim = '\t', col_types = cols(.default = col_character()))
      
      # Add sample name column
      data <- data |>
        mutate(SAMPLE = sample_name) |>
        distinct(SAMPLE, `MT Epitope Seq`, .keep_all = TRUE)
      
      # Add to list
      all_samples[[sample_name]] <- data
      
      out_df <- bind_rows(all_samples)
      
      cat(sprintf("Loaded: %s\n", sample_name))
    }
  }
  
  out_df <- out_df |> select(SAMPLE, everything())
  
  cat(sprintf("Combined %d samples into data frame\n", length(all_samples)))
  
  return(out_df)
}  

#function for plotting 
neoag_counts_plot <- function(out_df) {
  out_df <- out_df |>
    filter(!is.na(SAMPLE))
  
  burden_all <- out_df |>
    group_by(SAMPLE) |>
    summarise(Total_neoantigens = n(), .groups = 'drop')
  
  burden_classI <- out_df |>
    filter(str_detect(`HLA Allele`, "HLA")) |>
    group_by(SAMPLE) |>
    summarise(Class_I_neoantigens = n(), .groups = 'drop')
  
  burden_class2 <- out_df |>
    filter(str_detect(`HLA Allele`, "D")) |>
    group_by(SAMPLE) |>
    summarise(Class_II_neoantigens = n(), .groups = 'drop')
  
  # Join and pivot for plotting
  all_long <- burden_all |>
    left_join(burden_classI, by = "SAMPLE") |>
    left_join(burden_class2, by = "SAMPLE") |>
    add_column(Expression = "Yes") |>
    pivot_longer(
      cols = Total_neoantigens:Class_II_neoantigens,
      names_to = "Neoantigens",
      values_to = "Count"
    )
  
  # Create plot
  p <- ggplot(all_long, aes(x = Neoantigens, y = Count)) +
    geom_boxplot(fill = "steelblue") +
    geom_point(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.7), alpha = 0.7) +
    theme_bw()
  
  # Calculate and print statistics
  med_total <- median(burden_all$Total_neoantigens)
  sd_total <- sd(burden_all$Total_neoantigens)
  med_classI <- median(burden_classI$Class_I_neoantigens)
  sd_classI <- sd(burden_classI$Class_I_neoantigens)
  med_class2 <- median(burden_class2$Class_II_neoantigens)
  sd_class2 <- sd(burden_class2$Class_II_neoantigens)
  
  cat(sprintf("Overall median ± sd = %.1f ± %.1f\n", med_total, sd_total))
  cat(sprintf("Class I neoag median ± sd = %.1f ± %.1f\n", med_classI, sd_classI))
  cat(sprintf("Class II neoag median ± sd = %.1f ± %.1f\n", med_class2, sd_class2))
  
  # Return results as a list
  return(list(
    long_table = all_long,
    plot = p,
    statistics = data.frame(
      Type = c("Total", "Class I", "Class II"),
      Median = c(med_total, med_classI, med_class2),
      SD = c(sd_total, sd_classI, sd_class2)
    )
  ))
}

comparing_neoag_paired <- function(out_df_relapse, out_df_baseline) {
  graphs_rel <- neoag_counts_plot(out_df_relapse)
  
  rel_long <- graphs_rel$long_table |>
    mutate(Condition = "Relapse")
  
  graphs_base <- neoag_counts_plot(out_df_baseline)
  
  base_long <- graphs_base$long_table |>
    mutate(Condition = "Baseline") |>
    filter(SAMPLE %in% c(rel_long$SAMPLE))
  
  long_base_rel <- bind_rows(base_long, rel_long)
  
  p <- ggplot(long_base_rel, aes(x = Condition, y = Count, fill = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_brewer(palette = "Set1") +
    
    geom_line(aes(group = SAMPLE),
              colour = "grey50", alpha = 0.5) +
    
    geom_point(position = position_jitter(width = 0.1),
               alpha = 0.7) +
    
    facet_wrap(~ Neoantigens) +
    theme_bw()
  
  return(p)
}

recurrence <- function(out_df, epitope_or_HGVSp_or_both) {
  # Count unique patients per neoantigen
  if (epitope_or_HGVSp_or_both == "both") {
    result <- out_df |>
      group_by(`MT Epitope Seq`, HGVSp) |>
      summarise(
        Patient_count = n_distinct(SAMPLE),
        Total_count = n(),
        .groups = 'drop'
      ) |>
      mutate(Neoantigen = paste0(`MT Epitope Seq`, " (", HGVSp, ")")) |>
      select(Neoantigen, Patient_count, Total_count) |>
      filter(Patient_count > 1) |>
      arrange(desc(Patient_count))
    
  } else if (epitope_or_HGVSp_or_both == "epitope") {
    result <- out_df |>
      group_by(`MT Epitope Seq`) |>
      summarise(
        Patient_count = n_distinct(SAMPLE),
        Total_count = n(),
        .groups = 'drop'
      ) |>
      rename(Neoantigen = `MT Epitope Seq`) |>
      filter(Patient_count > 1) |>
      arrange(desc(Patient_count))
    
  } else if (epitope_or_HGVSp_or_both == "HGVSp") {
    result <- out_df |>
      group_by(HGVSp) |>
      summarise(
        Patient_count = n_distinct(SAMPLE),
        Total_count = n(),
        .groups = 'drop'
      ) |>
      rename(Neoantigen = HGVSp) |>
      filter(Patient_count > 1) |>
      arrange(desc(Patient_count))
    
  } else {
    stop("epitope_or_HGVSp_or_both must be 'epitope' or 'HGVSp' or 'both'")
  }
  
  # Get summary statistics
  total_patients <- n_distinct(out_df$SAMPLE)
  recurring_neoags <- nrow(result)
  
  # Calculate total unique neoantigen (recurrent + non-recurrent)
  if (epitope_or_HGVSp_or_both == "both") {
    total_unique_neoags <- out_df |>
      distinct(`MT Epitope Seq`, HGVSp) |>
      nrow()
  } else if (epitope_or_HGVSp_or_both == "epitope") {
    total_unique_neoags <- out_df |>
      distinct(`MT Epitope Seq`) |>
      nrow()
  } else if (epitope_or_HGVSp_or_both == "HGVSp") {
    total_unique_neoags <- out_df |>
      distinct(HGVSp) |>
      nrow()
  }
  
  # Create horizontal bar plot
  p <- ggplot(result, aes(x = reorder(Neoantigen, Patient_count), y = Patient_count)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(
      title = paste0("Neoantigen Recurrence across Patients"),
      subtitle = paste0("Total Patients: ", total_patients, " | Total Unique Neoags: ", total_unique_neoags, 
                        " | Recurring Neoags (>1 patient): ", recurring_neoags),
      x = "Neoantigen",
      y = "Number of Patients"
    ) +
    theme_bw() +
    theme(
      plot.subtitle = element_text(size = 10, face = "italic"),
      axis.text.y = element_text(size = 9)
    )
  
  # Return results as a list
  return(list(
    data = result,
    plot = p,
    summary = data.frame(
      Total_Patients = total_patients,
      Total_Unique_Neoags = total_unique_neoags,
      Recurring_Neoags = recurring_neoags,
      Non_Recurring_Neoags = total_unique_neoags - recurring_neoags,
      Neoags_per_Patient = round(total_unique_neoags / total_patients, 2)
    )
  ))
  
}

analyses <- function(out_df) {
  
  pvac <- out_df |>
    mutate(
      `Median MT IC50 Score` = as.numeric(`Median MT IC50 Score`),
      `Median Fold Change`   = as.numeric(`Median Fold Change`),
      `Peptide Length`       = as.integer(`Peptide Length`)
    ) |>
    
    mutate(
      Binder = case_when(
        `Median MT IC50 Score` < 50  ~ "Strong (IC50 < 50)",
        `Median MT IC50 Score` < 500 ~ "Intermediate (IC50 < 500)",
        TRUE                         ~ "Weak"
      )
    ) |>
    
    mutate(
      Fold_change_bins = case_when(
        `Median Fold Change` < 0.9 ~ "FC < 0.9",
        `Median Fold Change` < 1.1 ~ "FC 0.9–1.1",
        `Median Fold Change` < 3   ~ "FC 1.1–3",
        TRUE                       ~ "FC > 3"
      )
    ) |>
    
    mutate(
      HLA_class = case_when(
        str_starts(`HLA Allele`, "HLA-") ~ "Class I",
        str_starts(`HLA Allele`, "D")    ~ "Class II",
        TRUE                             ~ NA_character_
      )
    ) |>
      mutate(
        Fold_change_bins = factor(
          Fold_change_bins,
          levels = c(
            "FC < 0.9",
            "FC 0.9–1.1",
            "FC 1.1–3",
            "FC > 3"
          )
        )
      )
  
  IC50_dist <- ggplot(pvac, aes(x = `Median MT IC50 Score`)) +
    geom_histogram(bins = 50, fill = "steelblue") +
    #scale_x_log10() +
    theme_bw() +
    labs(x = "Median MT IC50", y = "Count") +
    theme(text = element_text(size = 20))
  
  IC50_strength <- ggplot(pvac, aes(x = Binder, fill = Binder)) +
    geom_bar(fill = "steelblue") +
    scale_fill_brewer(palette = "Blues", direction = -1) +
    theme_bw() +
    theme(text = element_text(size = 20))
  
  IC50_Fc <- ggplot(pvac, aes(x = `Median Fold Change`)) +
    geom_histogram(bins = 50,  fill = "steelblue") +
    theme_bw() +
    labs(x = "Median Fold Change") +
    theme(text = element_text(size = 20))

  IC50_FC_bins <- ggplot(pvac, aes(x = Fold_change_bins, fill = Fold_change_bins)) +
    geom_bar(fill = "steelblue") +
    theme_bw() +
    labs(x = "Median fold change in IC50 (relative to REF)") +
    theme(text = element_text(size = 20))
  
  epitope_length <- ggplot(pvac, aes(x = `Peptide Length`, fill = HLA_class)) +
    scale_fill_brewer(palette = "Set1") +
    geom_bar() +
    theme_bw() +
    scale_x_continuous(breaks = sort(unique(pvac$`Peptide Length`))) +
    labs(x = "Epitope Length", y = "Count") +
    theme(text = element_text(size = 20))
  
  return(list(
    IC50_dist = IC50_dist,
    IC50_strength = IC50_strength,
    IC50_Fc = IC50_Fc,
    IC50_FC_bins = IC50_FC_bins,
    epitope_length = epitope_length
    )
  )
}

#####################################################################
##########Applying functions
#####################################################################
setwd('/Volumes/DATA/DGE/DUDGE/MOPOPGEN/PROJECTS/CG2025P032_Myeloma_neoantigens/CGA0084_Neoantigen_prediction')

# Read all filtered TSV files from MUK9 directory and combine with sample names
#muk9_dir <- './MUK9_expressionfiltered/Relapse_expressed' #give path to folders with pvacseq output 
muk9_dir <- './MUK9_NOexpression/Relapse_noexpression'
sample_dirs <- list.dirs(muk9_dir, recursive = FALSE)

neoag_MUK9_relapse <- merging_pvacseq_combined_results(sample_dirs)

graphs_rel <- neoag_counts_plot(neoag_MUK9) 

test <- recurrence(neoag_MUK9_relapse, epitope_or_HGVSp_or_both = "both")

#####################################################
##########.  Baseline
#####################################################
muk9_dir <- './MUK9_NOexpression/Baseline_noexpression'
sample_dirs <- list.dirs(muk9_dir, recursive = FALSE)
neoag_MUK9_baseline <- merging_pvacseq_combined_results(sample_dirs)

neoag_MUK9_baseline <- neoag_MUK9_baseline |>
  filter(SAMPLE %in% c(rel_long$SAMPLE))

graphs_base <- neoag_counts_plot(neoag_MUK9_baseline) 
test <- recurrence(neoag_MUK9_baseline, epitope_or_HGVSp_or_both = "both")

comparing_neoag_paired(neoag_MUK9_relapse, neoag_MUK9_baseline)


graphs <- analyses(neoag_MUK9_baseline)
