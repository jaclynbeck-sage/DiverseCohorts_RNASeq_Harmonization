# This script is set up so that sections of it will be embedded in the Quarto
# QC notebooks. All this code is here instead of directly in the notebooks
# because there are 4 datasets / 4 notebooks, and moving the code here allows me
# to change something in one place and have it propagate to all 4 notebooks
# instead of having to change it 4 times.

# ---- include-libraries ----

library(synapser)
library(ggplot2)
library(viridis)
library(patchwork)
library(matrixStats)
library(dplyr)
library(stringr)
library(sageRNAUtils)

configs <- config::get(file = "config.yml")

# Defaults that should be over-written in the qmd file
dataset <- NULL
upload_to_synapse <- FALSE


# ---- download-metadata ----

# Function to copy DV200 and RIN values from original sample to sample swaps.
# In the assay metadata, specimens that were sample swaps have these fields
# set to NA, but the originating specimen has this information. Here we find
# each originating specimen and copy its RIN and DV200 over to the corresponding
# sample swap specimen. If there is more than 1 matching specimen with a non-NA
# value, the non-NA values are averaged to produce an imputed value.
fill_missing_numeric <- function(metadata, col_name) {
  # All rows with NA values in the field
  missing <- subset(metadata, is.na(metadata[, col_name]))

  for (row_ind in 1:nrow(missing)) {
    # All rows from the same individual and from the same tissue
    m_rows <- subset(metadata,
                     individualID == missing$individualID[row_ind] &
                       tissue == missing$tissue[row_ind])

    # Average all non-NA values that match the specimen with the missing value.
    # In most cases, there will only be one matching non-NA value, but for a few
    # samples there might be 2 non-NA values, in cases where two centers
    # processed the RNA separately.
    # Note: If this sample is not from a sample swap, the value will still be NA
    missing[row_ind, col_name] <- mean(m_rows[, col_name], na.rm = TRUE)
  }

  # Fill metadata with replacements for the missing values
  metadata[is.na(metadata[, col_name]), col_name] <- missing[, col_name]
  return(metadata)
}


# Function to copy batch values from original sample to sample swaps. In the
# assay metadata, specimens that were sample swaps have these fields set to NA,
# but the originating specimen has this information. Here we find each
# originating specimen and copy its rnaBatch over to the corresponding sample
# swap specimen.
fill_missing_rnaBatch <- function(metadata, col_name) {
  # All rows with NA values in the field
  missing <- subset(metadata, is.na(metadata[, col_name]))

  for (row_ind in 1:nrow(missing)) {
    # All rows from the same individual and from the same tissue
    m_rows <- subset(metadata, individualID == missing$individualID[row_ind] &
                       tissue == missing$tissue[row_ind])

    # This specimen is highly likely to be a sample swap that was not prepped at
    # the data generation site. Find the originating site's batch information
    m_rows <- subset(m_rows,
                     # Rush, Mayo
                     (dataGenerationSite == dataContributionGroup) |
                       # Mayo also processed Emory
                       (dataGenerationSite == "Mayo" & dataContributionGroup == "Emory") |
                       # NYGC processed Columbia and MSSM
                       (dataGenerationSite == "NYGC" & dataContributionGroup %in% c("Columbia", "MSSM")))

    # If exactly 1 unique non-NA value was found from the originating site, use that value
    if (length(unique(na.omit(m_rows[, col_name]))) == 1) {
      missing[row_ind, col_name] <- unique(na.omit(m_rows[, col_name]))
    }
  }

  # Fill metadata with replacements for the missing values
  metadata[is.na(metadata[, col_name]), col_name] <- missing[, col_name]
  return(metadata)
}


# Download the individual, biospecimen, and assay metadata
synLogin()

ind <- synGet(configs$download$individual_metadata_synid,
              downloadLocation = "downloads",
              ifcollision = "overwrite.local")$path |>
  read.csv()
bio <- synGet(configs$download$biospecimen_metadata_synid,
              downloadLocation = "downloads",
              ifcollision = "overwrite.local")$path |>
  read.csv()
assay <- synGet(configs$download$assay_metadata_synid,
                downloadLocation = "downloads",
                ifcollision = "overwrite.local")$path |>
  read.csv()

# Rush has 7 samples in a specific batch that are duplicates and should be removed
duplicates_remove <- configs$Rush$remove_specimenIDs # 7 samples
duplicates_batch <- configs$Rush$remove_specimenIDs_batch # B74

# Combine the 3 metadatas and fix a few fields
metadata <- assay |>
  merge(bio) |>
  merge(ind) |>

  # Remove 7 duplicate Rush samples that are in batch B74
  subset(!(specimenID %in% duplicates_remove & sequencingBatch == duplicates_batch)) |>

  # Remove the single STG sample from Columbia
  subset(specimenID != configs$Columbia$stg_remove_id) |>

  # Fix or alter some fields
  mutate(
    # Make PMI numeric
    PMI = suppressWarnings(as.numeric(PMI)),

    # Fix missing batch information for MSSM sample swaps sequenced at NYGC.
    # Batch information was confirmed by NYGC.
    rnaBatch = case_when(
      is.na(rnaBatch) & dataGenerationSite == "NYGC" &
        sampleExchangeOrigin == "MSSM" ~ "B01",
      .default = rnaBatch
    ),
    libraryBatch = case_when(
      is.na(libraryBatch) & dataGenerationSite == "NYGC" &
        sampleExchangeOrigin == "MSSM" ~ "B04",
      specimenID == "213917-2" ~ "B04",
      is.na(libraryBatch) & dataGenerationSite == "NYGC" &
        individualID %in% c("29600", "29631", "R1263350", "R3914030") ~ "B01",
      is.na(libraryBatch) & dataGenerationSite == "NYGC" &
        individualID %in% c("6434", "R7674931", "R9500594") ~ "B02",
      is.na(libraryBatch) & dataGenerationSite == "NYGC" &
        individualID %in% c("28871", "29691", "R5508487", "R9652199") ~ "B03",
      .default = libraryBatch
    ),

    # Batches need to be re-named to be unique to each data set --
    # * Mayo sequenced "Mayo", "Emory", and sample swaps in the same batches
    # * Rush did the same with "Rush" and sample swaps
    # * NYGC sequenced "Columbia" separately from "MSSM"
    across(c(rnaBatch, libraryBatch, sequencingBatch),
      ~case_match(
        dataGenerationSite,
        # Pre-pend "Mayo" or "Rush" to the batch, but leave NAs alone
        c("Mayo", "Rush") ~ ifelse(is.na(.x), .x,
                                   paste(dataGenerationSite, .x, sep = "_")),
        # Pre-pend "NYGC" + "MSSM" or "Columbia" to the batch, but leave NAs alone
        "NYGC" ~ case_when(
          is.na(.x) ~ .x,
          dataContributionGroup == "Columbia" ~ paste("NYGC_Columbia", .x, sep = "_"),
          # All MSSM and sample swaps
          .default = paste("NYGC_MSSM", .x, sep = "_")
        )
      )
    ),

    # Make specimenID match column names of count matrix -- needs to be done
    # last so we don't break references to original specimenIDs above
    specimenID = make.names(specimenID)
  )

# Fill missing batch information where it was left out for sample swaps
metadata <- fill_missing_rnaBatch(metadata, "rnaBatch")

# Fill missing RIN and DV200, most of which are from sample swaps
metadata <- fill_missing_numeric(metadata, "DV200")
metadata <- fill_missing_numeric(metadata, "RIN")


# ---- download-counts ----

# Download the count matrix from Synapse. The "dataset" variable must be set
# prior to including this code and should match one of the data set names in
# config.yml.
syn_ids <- configs[[dataset]]$count_matrix_synids
counts <- lapply(syn_ids, function(syn_id) {
  synGet(syn_id, downloadLocation = "downloads")$path |>
    read.table(header = TRUE) |>

    # RSEM adds a transcript ID column that we don't need
    dplyr::select(-transcript_id.s.) |>

    # Genes that end with version number followed by _PAR_Y should be removed,
    # as the counts are identical to their non-PAR_Y counterparts
    dplyr::filter(!grepl("\\.[0-9]+_PAR_Y", gene_id)) |>

    # Set the rownames to the gene_id column and remove the column
    tibble::column_to_rownames("gene_id") |>
    as.data.frame()
})

# Check for duplicate samples
samps <- unlist(lapply(counts, colnames))
stopifnot(all(table(samps) == 1))

counts <- purrr::list_cbind(counts) |>
  as.matrix()

# Make sure metadata and counts contain the same samples in the same order
metadata <- subset(metadata, specimenID %in% colnames(counts))
counts <- counts[, metadata$specimenID]

stopifnot(length(unique(metadata$specimenID)) == nrow(metadata))
stopifnot(all(metadata$specimenID == colnames(counts)))

# Normalize the counts matrix
orig_size <- ncol(counts)
counts_log <- sageRNAUtils::simple_log2norm(counts)


# ---- download-qc-stats ----

# These files are generated by running "01_Download_QC_Files.R" first
fastqc_data <- readRDS(file.path("data", "QC",
                                 paste0(dataset, "_fastqc_stats.rds")))
multiqc_stats <- readRDS(file.path("data", "QC",
                                   paste0(dataset, "_multiqc_stats.rds")))

gene_file <- synGet(configs$download$gene_metadata_synid,
                    downloadLocation = "downloads")
gene_info <- read.csv(gene_file$path)

fastqc_data <- lapply(fastqc_data, function(df) {
  merge(dplyr::select(metadata, specimenID, tissue), df)
})

# Columbia and Rush have to alter specimen IDs before merging, taken care of elsewhere
if (dataset != "Columbia" & dataset != "Rush") {
  multiqc_stats <- merge(dplyr::select(metadata, specimenID, tissue), multiqc_stats)
}


# ---- columbia-remap-ids ----

# Columbia only -- remap specimen IDs in the mulitqc stats data frame
if (dataset == "Columbia") {
  # New vs old ID information is contained in the fastqc filenames
  id_map <- fastqc_data$basic_statistics |>
    select(specimenID, Filename) |>
    mutate(old_id = str_replace(Filename, "_(1|2).gz", ""),
           old_id = make.names(old_id)) |>
    select(-Filename) |>
    distinct()

  multiqc_stats <- multiqc_stats |>
    dplyr::rename(old_id = specimenID) |>
    merge(id_map) |>
    select(-old_id)

  multiqc_stats <- merge(dplyr::select(metadata, specimenID, tissue), multiqc_stats)
}


# ---- rush-remove-duplicates ----

# Rush only -- Remove duplicate samples
if (dataset == "Rush") {
  to_remove <- lapply(configs$Rush$remove_samples_fastqc, function(id) {
    grep(id, fastqc_data$basic_statistics$sample, value = TRUE)
  })

  fastqc_data <- lapply(fastqc_data, function(df) {
    subset(df, !(sample %in% unlist(to_remove)))
  })

  multiqc_stats <- subset(multiqc_stats, !(specimenID %in% configs$Rush$remove_samples_fastqc)) |>
    mutate(specimenID = str_replace(specimenID, "_S[0-9]+", ""))
  multiqc_stats <- merge(dplyr::select(metadata, specimenID, tissue), multiqc_stats)
}


# ---- validate-fastqc ----

thresholds <- configs$thresholds

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

n_fqc_warn_fail <- sum(metadata$base_content_warn | metadata$phred_score_warn |
                         !metadata$phred_score_valid)


# ---- print-fastqc-results ----

meta_sub <- subset(metadata, base_content_warn | phred_score_warn | !phred_score_valid) |>
  select(specimenID, tissue, base_content_warn, phred_score_warn, phred_score_valid) |>
  mutate(
    `Phred Score` = case_when(
      !phred_score_valid ~ "FAIL",
      phred_score_warn ~ "WARN",
      .default = "."
    ),
    `Base Content` = ifelse(base_content_warn == TRUE, "WARN", ".")
  ) |>
  select(specimenID, tissue, `Phred Score`, `Base Content`) |>
  dplyr::rename(Tissue = tissue, `Specimen ID` = specimenID) |>
  arrange(Tissue, `Specimen ID`)

if (nrow(meta_sub) > 0) {
  meta_sub
}

# ---- validate-multiqc ----

thresholds <- configs$thresholds

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

mqc_plot <- multiqc_stats |>
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
plt1 <- ggplot(mqc_plot,
               aes(x = tissue, y = samtools_reads_mapped_percent, fill = tissue)) +
  geom_boxplot(outliers = FALSE, width = 0.1) +
  geom_jitter(aes(color = mapped_status),
              width = 0.3,
              size = ifelse(mqc_plot$mapped_status == "Pass", 0.5, 1)) +
  theme_bw() +
  scale_color_manual(values = stat_colors) +
  ggtitle("Percentage of reads mapped")

plt2 <- ggplot(mqc_plot,
               aes(x = tissue, y = picard_PERCENT_DUPLICATION, fill = tissue)) +
  geom_boxplot(outliers = FALSE, width = 0.1) +
  geom_jitter(aes(color = duplicated_status),
              width = 0.3,
              size = ifelse(mqc_plot$duplicated_status == "Pass", 0.5, 1)) +
  theme_bw() +
  scale_color_manual(values = stat_colors) +
  ggtitle("Percentage of reads duplicated")

print(plt1 + plt2)

n_mqc_warn_fail <- sum(metadata$reads_mapped_warn |
                         metadata$reads_duplicated_warn |
                         !metadata$reads_mapped_valid |
                         !metadata$reads_duplicated_valid)


# ---- print-multiqc-results ----

meta_sub <- subset(metadata, reads_mapped_warn | reads_duplicated_warn |
                     !reads_mapped_valid | !reads_duplicated_valid) |>
  select(specimenID, tissue, reads_mapped_warn, reads_duplicated_warn,
         reads_mapped_valid, reads_duplicated_valid) |>
  mutate(
    `Reads Mapped` = case_when(
      !reads_mapped_valid ~ "FAIL",
      reads_mapped_warn ~ "WARN",
      .default = "."
    ),
    `Reads Duplicated` = case_when(
      !reads_duplicated_valid ~ "FAIL",
      reads_duplicated_warn ~ "WARN",
      .default = "."
    )
  ) |>
  select(specimenID, tissue, `Reads Mapped`, `Reads Duplicated`) |>
  dplyr::rename(Tissue = tissue, `Specimen ID` = specimenID) |>
  arrange(Tissue, `Specimen ID`)

if (nrow(meta_sub) > 0) {
  meta_sub
}


# ---- validate-sex ----

thresholds <- configs$thresholds

mismatches <- sageRNAUtils::find_sex_mismatches(
  metadata, counts_log, y_expr_threshold = thresholds$sex
)

plts <- sageRNAUtils::plot_sex_mismatch_results(
  mismatches$sex_check_df, thresholds$sex, print_plot = FALSE
)
print(plts[[1]] + plts[[2]])

metadata$sex_valid <- !(metadata$specimenID %in% mismatches$mismatches)


# ---- print-sex-mismatches ----

meta_sub <- subset(metadata, !sex_valid) |>
  select(specimenID, tissue) |>
  dplyr::rename(`Specimen ID` = specimenID,
                Tissue = tissue) |>
  arrange(Tissue, `Specimen ID`)

if (nrow(meta_sub) > 0) {
  meta_sub
}


# ---- validate-pca ----

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

results <- sageRNAUtils::find_pca_outliers_by_group(
  counts_log,
  pca_group = "pca_group",
  n_sds = 4,
  metadata = metadata,
  gene_info = gene_info
)

plts <- lapply(names(results$group_results), function(res_name) {
  plt <- sageRNAUtils::plot_pca_outliers(
    results$group_results[[res_name]]$pca_df,
    results$group_results[[res_name]]$thresholds,
    print_plot = FALSE
  ) +
    ggtitle(res_name)
})

# TODO wrap plots differently
print(Reduce("+", plts))

metadata$pca_valid <- !(metadata$specimenID %in% results$outliers)


# ---- print-pca-outliers ----

meta_sub <- subset(metadata, !pca_valid) |>
  select(specimenID, tissue) |>
  dplyr::rename(`Specimen ID` = specimenID, Tissue = tissue) |>
  arrange(Tissue, `Specimen ID`)

if (nrow(meta_sub) > 0) {
  meta_sub
}


# ---- validate-dv200 ----

thresholds <- configs$thresholds

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

print(plt1 + plt2)
print(plt3 + plt4)

metadata$DV200_valid <- case_when(
  !is.na(metadata$DV200) ~ metadata$DV200 >= thresholds$DV200,
  # If DV200 is NA, use RIN instead
  is.na(metadata$DV200) & !is.na(metadata$RIN) ~ metadata$RIN >= thresholds$RIN,
  # If both are NA, fail QC
  .default = FALSE
)


# ---- print-dv200-results ----

meta_sub <- subset(metadata, !DV200_valid) |>
  select(specimenID, tissue, DV200, RIN) |>
  dplyr::rename(`Specimen ID` = specimenID, Tissue = tissue) |>
  arrange(Tissue, DV200, RIN)

if (nrow(meta_sub) > 0) {
  meta_sub
}


# ---- save-samples ----

metadata$valid <- Reduce("&", metadata[, grepl("_valid", colnames(metadata))])
metadata$warn <- Reduce("+", metadata[, grepl("_warn", colnames(metadata))])

metadata$valid <- metadata$valid & metadata$warn < 2

n_passes <- table(metadata$tissue, metadata$valid)
colnames(n_passes) <- c("Fail", "Pass")
failures <- subset(metadata, valid == FALSE)

metadata <- subset(metadata, valid == TRUE)
counts <- counts[, metadata$specimenID]
multiqc_stats <- subset(multiqc_stats, specimenID %in% metadata$specimenID)

# Remove genes that are all 0's
zeros <- rowSums(counts) == 0

data_final <- list("metadata" = metadata, "counts" = counts[!zeros, ],
                   "multiqc_stats" = multiqc_stats)

saveRDS(data_final, file.path("data", "QC",
                              paste0(dataset, "_qc.rds")))

# Save counts to CSV and upload to Synapse
counts_filename <- file.path("data", "counts_post_qc",
                             paste0(dataset, "_counts_filtered.csv"))

write.csv(counts[!zeros, ], counts_filename, quote = FALSE)

# Set `upload_to_synapse` to TRUE or FALSE in the qmd notebook
if (upload_to_synapse) {
  synLogin()

  syn_file <- File(counts_filename, parent = configs$upload$counts_folder_synid)

  provenance <- c(configs$download$individual_metadata_synid,
                  configs$download$biospecimen_metadata_synid,
                  configs$download$assay_metadata_synid,
                  configs$download$gene_metadata_synid,
                  configs[[dataset]]$fastq_folder_synids,
                  configs[[dataset]]$multiqc_json_synids,
                  configs[[dataset]]$count_matrix_synids)

  github <- c(
    paste0(
      "https://github.com/jaclynbeck-sage/DiverseCohorts_RNASeq_Harmonization/",
      "blob/main/02_", dataset, "_QC.qmd"
      ),
    paste0(
      "https://github.com/jaclynbeck-sage/DiverseCohorts_RNASeq_Harmonization/",
      "blob/main/qc_functions.R"
    )
  )

  syn_file <- synStore(
    syn_file,
    forceVersion = FALSE,
    used = provenance,
    executed = github
  )
}



# ---- print-final-qc-results ----

n_passes |>
  as.data.frame() |>
  tidyr::pivot_wider(names_from = Var2, values_from = Freq) |>
  dplyr::rename(Tissue = Var1)


# ---- print-qc-failures-summary ----

failures |>
  group_by(tissue) |>
  summarize(`Specimen IDs` = paste(str_replace(specimenID, "^X", ""),
                                   collapse = ", "))


# ---- print-qc-failures-detail ----

failures_detail <- failures |>
  select(specimenID, tissue, isSampleExchange, sampleExchangeOrigin, cohort,
         contains("_warn"), contains("_valid")) |>
  mutate(
    `Specimen ID` = str_replace(specimenID, "^X", ""),
    `Sample Exchange?` = case_when(
      isSampleExchange & sampleExchangeOrigin != dataset ~ paste0("Yes (", sampleExchangeOrigin, ")"),
      .default = "No"
    ),
    `Base content` = ifelse(base_content_warn == TRUE, "WARN", "."),
    `Phred score` = case_when(
      !phred_score_valid ~ "FAIL",
      phred_score_warn ~ "WARN",
      .default = "."
    ),
    `Reads Mapped` = case_when(
      !reads_mapped_valid ~ "FAIL",
      reads_mapped_warn ~ "WARN",
      .default = "."
    ),
    `Reads Duplicated` = case_when(
      !reads_duplicated_valid ~ "FAIL",
      reads_duplicated_warn ~ "WARN",
      .default = "."
    ),
    `Sex check` = ifelse(sex_valid, ".", "FAIL"),
    `PCA check` = ifelse(pca_valid, ".", "FAIL"),
    `DV200 check` = ifelse(DV200_valid, ".", "FAIL")
  ) |>
  select(`Specimen ID`, tissue, cohort, `Sample Exchange?`,
         `Base content`, `Phred score`, `Reads Mapped`, `Reads Duplicated`,
         `Sex check`, `PCA check`, `DV200 check`) |>
  dplyr::rename(Tissue = tissue, Cohort = cohort) |>
  arrange(Tissue, `Specimen ID`)

failures_detail
