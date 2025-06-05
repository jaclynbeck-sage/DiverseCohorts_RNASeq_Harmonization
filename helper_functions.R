# ---- include-libraries ----
library(synapser)
library(ggplot2)
library(viridis)
library(patchwork)
library(matrixStats)
library(dplyr)
library(stringr)
library(sageRNAUtils)

# ---- download-data-functions ----

download_metadata <- function(configs) {
  synLogin()

  ind <- synGet(configs$individual_metadata_synid,
                downloadLocation = "downloads",
                ifcollision = "overwrite.local")$path |>
    read.csv()
  bio <- synGet(configs$biospecimen_metadata_synid,
                downloadLocation = "downloads",
                ifcollision = "overwrite.local")$path |>
    read.csv()
  assay <- synGet(configs$assay_metadata_synid,
                  downloadLocation = "downloads",
                  ifcollision = "overwrite.local")$path |>
    read.csv()

  duplicates_remove <- configs$Rush$remove_specimenIDs # 7 samples
  duplicates_batch <- configs$Rush$remove_specimenIDs_batch # B74

  metadata <- assay |>
    merge(bio) |>
    merge(ind) |>
    mutate(specimenID = make.names(specimenID)) |>
    # Remove 7 duplicate Rush samples that are in batch B74
    subset(!(specimenID %in% duplicates_remove & sequencingBatch == duplicates_batch)) |>
    # Add some extra statistics
    mutate(
      apoe4Status = case_when(
        grepl("4", apoeGenotype) ~ "True",
        apoeGenotype == "missing or unknown" ~ "missing or unknown",
        .default = "False"
      ),
      ageDeathNumeric = suppressWarnings(as.numeric(ageDeath)),
      PMI = suppressWarnings(as.numeric(PMI)),
      ageDeath_binned = case_when(
        ageDeathNumeric < 65 ~ "Under 65",
        ageDeathNumeric >= 65 & ageDeathNumeric < 70 ~ "65 to 69",
        ageDeathNumeric >= 70 & ageDeathNumeric < 75 ~ "70 to 74",
        ageDeathNumeric >= 75 & ageDeathNumeric < 80 ~ "75 to 79",
        ageDeathNumeric >= 80 & ageDeathNumeric < 85 ~ "80 to 74",
        ageDeathNumeric >= 85 & ageDeathNumeric < 90 ~ "85 to 90",
        .default = ageDeath
      ),
      # Sequencing batches need to be re-named to be unique to each data set --
      # Mayo sequenced "Mayo", "Emory", and sample swaps in the same batches
      # Rush did the same with "Rush" and sample swaps
      # NYGC sequenced "Columbia" separately from "MSSM" -- TODO sample swaps?
      # TODO this is wrong for MSSM, sample swaps were prepped at NYGC instead of at the originating site
      across(c(rnaBatch, libraryBatch, sequencingBatch),
        ~ case_match(dataGenerationSite,
          c("Mayo", "Rush") ~ ifelse(is.na(.x), .x, paste(dataGenerationSite, .x, sep = "_")),
          "NYGC" ~ ifelse(is.na(.x), .x, paste(dataGenerationSite, dataContributionGroup, .x, sep = "_"))
        )
      )
    ) |>
    select(-ageDeathNumeric)

  # Fill missing batch information where it was left out for sample swaps
  metadata <- fill_missing_categorical(metadata, "rnaBatch")
  metadata <- fill_missing_categorical(metadata, "libraryBatch")

  # Fill missing RIN and DV200, most of which are from sample swaps
  metadata <- fill_missing_numeric(metadata, "DV200")
  metadata <- fill_missing_numeric(metadata, "RIN")
}

fill_missing_numeric <- function(metadata, col_name) {
  missing <- subset(metadata, is.na(metadata[, col_name]))

  for (row_ind in 1:nrow(missing)) {
    m_rows <- subset(metadata, individualID == missing$individualID[row_ind] &
                       tissue == missing$tissue[row_ind])

    # Note: If this sample is not from a sample swap, the value will likely
    # still be NA
    missing[row_ind, col_name] <- mean(m_rows[, col_name], na.rm = TRUE)
  }

  metadata[is.na(metadata[, col_name]), col_name] <- missing[, col_name]
  return(metadata)
}


fill_missing_categorical <- function(metadata, col_name) {
  missing <- subset(metadata, is.na(metadata[, col_name]))

  for (row_ind in 1:nrow(missing)) {
    m_rows <- subset(metadata, individualID == missing$individualID[row_ind] &
                       tissue == missing$tissue[row_ind])

    # This is highly likely to be a sample swap that was not prepped at the
    # data generation site. Find the originating site's batch information
    m_rows <- subset(m_rows, (dataGenerationSite == dataContributionGroup) |
                       (dataGenerationSite == "Mayo" & dataContributionGroup == "Emory") |
                       (dataGenerationSite == "NYGC" & dataContributionGroup %in% c("Columbia", "MSSM")))

    # Note: Sample swaps that originated from MSSM do not currently have batch
    # information so this will remain NA for those samples
    if (length(na.omit(m_rows[, col_name])) == 1) {
      missing[row_ind, col_name] <- na.omit(m_rows[, col_name])
    }
  }

  metadata[is.na(metadata[, col_name]), col_name] <- missing[, col_name]
  return(metadata)
}


download_rsem <- function(syn_id) {
  synLogin()
  synGet(syn_id, downloadLocation = "downloads")$path |>
    read.table(header = TRUE) |>
    dplyr::select(-transcript_id.s.) |>
    # Genes that end with version number followed by _PAR_Y should be removed,
    # as the counts are identical to their non-PAR_Y counterparts
    dplyr::filter(!grepl("\\.[0-9]+_PAR_Y", gene_id)) |>
    tibble::column_to_rownames("gene_id") |>
    as.matrix()
}


# ---- download-qc ----

download_fastqc <- function(dataset_config, load_saved_stats = FALSE) {
  fastq_dir <- file.path("downloads", paste0(dataset_config$name, "_fastqc"))
  fq_stats_file <- file.path("data", "QC",
                             paste0(dataset_config$name, "_fq_stats_raw.rds"))

  # Avoid repeating the very long process of downloading and reading in the
  # fastqc files if possible
  if (load_saved_stats & file.exists(fq_stats_file)) {
    message("Loading pre-existing FastQC file...")
    fq_stats <- readRDS(fq_stats_file)
  } else {
    synLogin()

    # Get a list of all fastq.zip files in the folder on Synapse and download them
    children <- synGetChildren(dataset_config$fastqc_folder_synid)

    message("Downloading FastQC data from Synapse...")

    files <- sapply(children$asList(), function(child) {
      synGet(child$id, downloadLocation = fastq_dir)$path
    })

    # Save the names of the original files so that if this is Columbia data,
    # renaming the files doesn't lose the sample names we need
    sample_names <- str_replace(basename(files), "_fastqc.zip", "")

    if (all(grepl("NYBB_", files))) {
      files <- remap_columbia_fastqc(files)
    }

    # Save time by not loading unneeded stats from the files
    modules <- c("Basic Statistics", "Per base sequence quality",
                 "Per base sequence content", "Sequence Duplication Levels")
    fq_stats <- load_fastqc_output(files, sample_names, modules = modules)

    # Save in case we need to re-run some of this processing
    saveRDS(fq_stats, file = fq_stats_file)
  }

  message("Summarizing statistics...")

  # Summary and modifications to the statistics
  basic_stats <- fq_stats$basic_statistics |>
    dplyr::mutate(
      # This pattern matches strings that end in "_1", "_2", or that end in a
      # pattern like "_S21_1". Anything before these patterns should be the
      # specimen ID for every study's files
      specimenID = str_replace(sample, "(_S[0-9]+)?_(1|2)$", ""),
      specimenID = make.names(specimenID),
    ) |>
    dplyr::rename(total_sequences = `Total Sequences`,
                  percent_gc_content = `%GC`)

  id_df <- dplyr::select(basic_stats, sample, specimenID, read)

  phred_per_base <- merge(id_df, fq_stats$per_base_sequence_quality) |>
    dplyr::rename(position = Base)

  base_content <- merge(id_df, fq_stats$per_base_sequence_content) |>
    dplyr::rename(position = Base) |>
    tidyr::pivot_longer(cols = c("G", "A", "T", "C"),
                        names_to = "base",
                        values_to = "percent") |>
    dplyr::mutate(difference = abs(percent - 25)) |>
    dplyr::group_by(sample, specimenID, read, base) |>
    dplyr::summarize(mean_base_deviation = sum(difference),
                     .groups = "drop")

  return(list(basic_statistics = basic_stats,
              phred_per_base = phred_per_base,
              base_content = base_content))
}


# Columbia samples inside the fastqc files have a different sample name than the
# file name, which breaks fastqcr::read. Here we rename all the zip files to
# match the sample name inside. The original names are retained prior to calling
# this function so we can map new names onto old ones.
remap_columbia_fastqc <- function(fastqc_files) {
  new_files <- sapply(fastqc_files, function(zip_file) {
    # Lists what's in the zip file without extracting it
    contents <- unzip(zip_file, list = TRUE)
    old_name <- paste0(str_replace(contents$Name[1], "/", ""),
                       ".zip")
    old_name <- file.path(dirname(zip_file), old_name)

    # If we've done this before and there are already renamed files in the
    # directory, they will be overwritten by this command
    file.rename(zip_file, old_name)
    return(old_name)
  })

  return(new_files)
}


# ---- validate-fastqc-function ----

validate_fastqc <- function(metadata, fastqc_data, thresholds) {
  # TODO jitter the outlier points
  plt1 <- ggplot(fastqc_data$basic_statistics,
                 aes(x = tissue, y = percent_gc_content, fill = read)) +
    geom_boxplot(width = 0.5) +
    theme_bw() +
    ggtitle("Percent GC Content")

  print(plt1)

  plt2 <- ggplot(fastqc_data$base_content, aes(x = base, y = mean_base_deviation, fill = base)) +
    geom_boxplot(width = 0.5) +
    theme_bw() +
    facet_grid(rows = vars(read), cols = vars(tissue)) +
    ggtitle("Mean deviation from expected base proportions")

  print(plt2)

  # Using the mean content deviation from 25% of each base at each position.
  # Mayo uses sum of deviation instead of mean, but they're equivalent for the
  # purposes of outlier finding.
  base_content_outliers <- fastqc_data$base_content |>
    dplyr::group_by(tissue, read, base) |>
    mutate(is_outlier = is_outlier_IQR(mean_base_deviation, tail = "upper")) |>
    subset(is_outlier == TRUE) |>
    pull(specimenID)

  # Any base falling below the Phred threshold marks the sample as failing QC
  phred_fail <- fastqc_data$phred_per_base |>
    subset(Median < thresholds$phred) |>
    pull(specimenID)

  phred_outliers <- fastqc_data$phred_per_base |>
    group_by(position) |>
    mutate(is_outlier = is_outlier_IQR(Median, tail = "lower")) |>
    subset(is_outlier == TRUE) |>
    pull(specimenID)

  metadata$base_content_warn <- metadata$specimenID %in% base_content_outliers
  metadata$phred_score_valid <- !(metadata$specimenID %in% phred_fail)
  metadata$phred_score_warn <- metadata$specimenID %in% phred_outliers

  return(metadata)
}


# ---- validate-multiqc-function ----

# In the absence of RSeQC data, we only check percentage of reads mapped and
# percentage of duplicated reads for now.
validate_multiqc <- function(metadata, multiqc_stats, thresholds) {
  reads_mapped_fail <- multiqc_stats |>
    subset(samtools_reads_mapped_percent < thresholds$reads_mapped) |>
    pull(specimenID)

  reads_mapped_outliers <- multiqc_stats |>
    mutate(is_outlier = is_outlier_IQR(samtools_reads_mapped_percent, tail = "lower")) |>
    subset(is_outlier == TRUE) |>
    pull(specimenID)

  reads_dupe_fail <- multiqc_stats |>
    subset(picard_PERCENT_DUPLICATION > thresholds$reads_duplicated) |>
    pull(specimenID)

  # Use Q3 + 3*IQR for outliers here
  reads_dupe_outliers <- multiqc_stats |>
    mutate(is_outlier = is_outlier_IQR(picard_PERCENT_DUPLICATION,
                                       tail = "upper", IQR_mult = 3)) |>
    subset(is_outlier == TRUE) |>
    pull(specimenID)

  metadata$reads_mapped_valid <- !(metadata$specimenID %in% reads_mapped_fail)
  metadata$reads_mapped_warn <- metadata$specimenID %in% reads_mapped_outliers
  metadata$reads_duplicated_valid <- !(metadata$specimenID %in% reads_dupe_fail)
  metadata$reads_duplicated_warn <- metadata$specimenID %in% reads_dupe_outliers

  multiqc_stats <- multiqc_stats |>
    mutate(
      mapped_status = case_when(
        specimenID %in% reads_mapped_fail ~ "Fail",
        specimenID %in% reads_mapped_outliers ~ "Warn",
        .default = "Pass"
      ),
      duplicated_status = case_when(
        specimenID %in% reads_dupe_fail ~ "Fail",
        specimenID %in% reads_dupe_outliers ~ "Warn",
        .default = "Pass"
      )
    )

  stat_colors <- c("Pass" = "black", "Warn" = "orange", "Fail" = "red")
  plt1 <- ggplot(multiqc_stats,
                 aes(x = tissue, y = samtools_reads_mapped_percent, fill = tissue)) +
    geom_boxplot(outliers = FALSE, width = 0.1) +
    geom_jitter(aes(color = mapped_status),
                width = 0.3,
                size = ifelse(multiqc_stats$mapped_status == "Pass", 0.5, 1)) +
    theme_bw() +
    scale_color_manual(values = stat_colors) +
    ggtitle("Percentage of reads mapped")

  plt2 <- ggplot(multiqc_stats,
                 aes(x = tissue, y = picard_PERCENT_DUPLICATION, fill = tissue)) +
    geom_boxplot(outliers = FALSE, width = 0.1) +
    geom_jitter(aes(color = duplicated_status),
                width = 0.3,
                size = ifelse(multiqc_stats$duplicated_status == "Pass", 0.5, 1)) +
    theme_bw() +
    scale_color_manual(values = stat_colors) +
    ggtitle("Percentage of reads duplicated")

  print(plt1 + plt2)

  return(metadata)
}


# ---- validate-sex-function ----

validate_sex <- function(metadata, data, thresholds) {
  mismatches <- find_sex_mismatches(metadata, data,
                                    y_expr_threshold = thresholds$sex)

  plts <- plot_sex_mismatch_results(mismatches$sex_check_df, thresholds$sex,
                                    print_plot = FALSE)
  print(plts[[1]] + plts[[2]])

  metadata$sex_valid <- !(metadata$specimenID %in% mismatches$mismatches)

  metadata
}


# ---- validate-pca-function ----

validate_pca <- function(metadata, data, gene_info, n_sds = 4) {
  # Do PCA outlier detection on a per-tissue or per-group basis.
  # Note: When doing a PCA of Rush data, the points clearly separate by batch.
  # This is most evident in the DLPFC but is also mildly visible in the other
  # tissues. Combining all batches causes batch-specific outliers to be missed,
  # so instead Rush outliers are detected on a tissue + batch basis.

  if (unique(metadata$dataGenerationSite) == "Rush") {
    metadata$pca_group <- paste(metadata$tissue, metadata$sequencingBatch, sep = " / ")
  } else {
    metadata$pca_group <- metadata$tissue
  }

  results <- find_pca_outliers_by_group(data,
                                        pca_group = "pca_group",
                                        n_sds = n_sds,
                                        metadata = metadata,
                                        gene_info = gene_info)

  plts <- lapply(names(results$group_results), function(res_name) {
    plt <- plot_pca_outliers(results$group_results[[res_name]]$pca_df,
                             results$group_results[[res_name]]$pc1_threshold,
                             results$group_results[[res_name]]$pc2_threshold,
                             print_plot = FALSE) +
      ggtitle(res_name)
  })

  print(Reduce("+", plts))

  metadata$pca_valid <- !(metadata$specimenID %in% results$outliers)

  metadata
}


# ---- validate-dv200-function ----

validate_DV200 <- function(metadata, thresholds) {
  plt1 <- ggplot(metadata, aes(x = tissue, y = RIN, fill = tissue)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(size = 0.5) +
    theme_bw()

  plt2 <- ggplot(metadata, aes(x = tissue, y = DV200, fill = tissue)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(size = 0.5) +
    theme_bw()

  plt3 <- ggplot(metadata, aes(x = RIN, y = DV200, color = tissue)) +
    geom_jitter(size = 0.5) +
    theme_bw()

  plt4 <- ggplot(metadata, aes(x = rank(DV200), y = DV200, color = tissue)) +
    geom_point(size = 0.5) +
    geom_hline(yintercept = thresholds$DV200) +
    theme_bw()

  print((plt1 + plt2) / (plt3 + plt4))

  metadata$DV200_valid <- case_when(
    !is.na(metadata$DV200) ~ metadata$DV200 > thresholds$DV200,
    # If DV200 is NA, use RIN instead
    is.na(metadata$DV200) & !is.na(metadata$RIN) ~ metadata$RIN > thresholds$RIN,
    # If both are NA, fail QC
    .default = FALSE)

  metadata
}


# ---- download-multiqc-function ----
# TODO move to library?
download_multiqc_json <- function(syn_id) {
  synLogin()
  js_file <- synGet(syn_id, downloadLocation = "downloads")

  js_txt <- readr::read_file(js_file$path)
  js_txt <- str_replace_all(js_txt, "NaN", str_escape('"NaN"'))
  js <- jsonlite::fromJSON(js_txt, simplifyVector = FALSE)

  data <- lapply(js$report_saved_raw_data, function(item) {
    do.call(rbind, item)
  })

  # Helper function to convert items from "data", which are matrices, to
  # data frames with only the specified columns. Rownames are moved to a new
  # "specimenID" column and the specified prefix is added to the other column
  # names
  reformat_stats <- function(df, prefix, cols_keep) {
    df |>
      as.data.frame() |>
      dplyr::select({{cols_keep}}) |>
      # Add prefix to column names
      dplyr::rename_with(~ paste0(prefix, "_", .x), everything()) |>
      # All columns are actually named lists, unlist them
      dplyr::mutate(across(everything(), ~unlist(.x))) |>
      # Make specimenID column
      tibble::rownames_to_column("specimenID") |>
      dplyr::mutate(specimenID = make.names(specimenID))
  }

  # TODO decide on:
  # multiqc_samtools_flagstat: ??
  # multiqc_fastqc: %GC, total_deduplicated_percentage <-- identical to fastqc data, but phred and base content not available

  picard_stats <- data$multiqc_picard_dups |>
    reformat_stats(prefix = "picard", cols_keep = PERCENT_DUPLICATION) |>
    dplyr::mutate(picard_PERCENT_DUPLICATION = picard_PERCENT_DUPLICATION * 100)

  rsem_stats <- data$multiqc_rsem |>
    reformat_stats(prefix = "rsem",
                   cols_keep = c(alignable_percent, Unique, Total)) |>
    dplyr::mutate(
      rsem_uniquely_aligned_percent = rsem_Unique / rsem_Total * 100
    ) |>
    dplyr::select(-rsem_Unique, -rsem_Total)

  samtools_stats <- data$multiqc_samtools_stats |>
    reformat_stats(
      prefix = "samtools",
      cols_keep = c(average_quality, insert_size_average, reads_mapped_percent,
                    reads_duplicated_percent, reads_MQ0_percent,
                    reads_QC_failed_percent)
    )

  # Each item in idxstats is a list with 2 numbers. The first number is the
  # number of reads mapped to the chromosome, the second number is the length
  # of the chromosome. We want the first number for X and Y
  xy_stats <- apply(data$multiqc_samtools_idxstats, 1, function(row) {
    c(chrX = row[["chrX"]][[1]], chrY = row[["chrY"]][[1]])
  }) |>
    t() |>
    reformat_stats(prefix = "samtools", cols_keep = c(chrX, chrY)) |>
    dplyr::mutate(total = samtools_chrX + samtools_chrY,
                  samtools_percent_mapped_X = samtools_chrX / total * 100,
                  samtools_percent_mapped_Y = samtools_chrY / total * 100) |>
    dplyr::select(-total, -samtools_chrX, -samtools_chrY)

  # This matrix has separate rows for read 1 and read 2, which need to be
  # combined into a single row.
  cutadapt_stats <- data$multiqc_cutadapt |>
    reformat_stats(prefix = "cutadapt", cols_keep = percent_trimmed) |>
    dplyr::mutate(read = ifelse(grepl("_1$", specimenID), "R1", "R2"),
                  specimenID = str_replace(specimenID, "_(1|2)$", "")) |>
    tidyr::pivot_wider(names_from = "read",
                       values_from = "cutadapt_percent_trimmed") |>
    dplyr::rename(cutadapt_percent_trimmed_R1 = R1,
                  cutadapt_percent_trimmed_R2 = R2) |>
    dplyr::mutate(
      cutadapt_mean_percent_trimmed = (cutadapt_percent_trimmed_R1 +
                                         cutadapt_percent_trimmed_R2) / 2
    )

  technical_stats <- purrr::reduce(
    list(picard_stats, rsem_stats, samtools_stats, xy_stats, cutadapt_stats),
    dplyr::full_join
  )

  if (any(is.na(technical_stats))) {
    warning("NA values exist in technical stats data frame.")
  }

  return(technical_stats)

  # Possible QC metrics:
  # FastQC -
  #   median PHRED
  #   base content at each position
  # RSeqC -
  #   reads explained by first strand
  #   percent mapped reads
  #   number of reads mapped to genes
  #   percent reads mapped to junctions
  #   percent reads mapped to genes
  #   insert size inner distance
  #   percent known junctions
  #   80/20 ratio of gene body coverage
}


# ---- remove-correlated-vars-function ----
# TODO move to library
remove_correlated_variables <- function(cor_mat, na_vars, R2_threshold = 0.5) {
  r2 <- cor_mat^2
  diag(r2) <- NA

  r2_melt <- r2
  r2_melt[upper.tri(r2_melt, diag = TRUE)] <- 0  # Avoids picking up both (a vs b) and (b vs a)
  r2_melt <- r2_melt |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "var1") |>
    tidyr::pivot_longer(cols = -var1,
                        names_to = "var2", values_to = "value") |>
    subset(value >= R2_threshold) |>  # pairs with R^2 > 0.5 (cor ~ 0.7) only
    dplyr::arrange(desc(value))

  to_remove <- c()

  for (R in 1:nrow(r2_melt)) {
    vars <- as.character(c(r2_melt$var1[R], r2_melt$var2[R]))

    if (any(vars %in% to_remove)) {
      # No need to re-check if we're already removing one of these variables
      next
    } else {
      # If either variable has any NA values, and they don't have the same number
      # of NAs, remove the one with the most NAs
      if (any(na_vars[1, vars] > 0) && (na_vars[1, vars[1]] != na_vars[1, vars[2]])) {
        to_remove <- c(to_remove,
                       vars[which.max(na_vars[1, vars])])
      } else {
        cur_vars <- setdiff(rownames(r2), to_remove)
        mean_r2 <- rowMeans(r2[cur_vars, cur_vars], na.rm = TRUE)

        # Remove the variable with the largest mean R^2 with the other remaining variables
        to_remove <- c(to_remove,
                       vars[which.max(mean_r2[vars])])
      }
    }
  }

  return(unique(to_remove))
}
