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
library(forcats)
library(plotly)

configs <- config::get(file = "config.yml")

synLogin()


# ---- not-run ----

# Defaults that should be over-written in the qmd file. This is only here to
# allow for testing the script directly outside of the QC notebooks.
dataset <- NULL
upload_to_synapse <- FALSE
render_for_wiki <- FALSE


# ---- optional-wiki-setup ----

if (render_for_wiki) {
  # This resolution and aspect ratio, scaled to ~45% on the wiki, looks good
  knitr::opts_chunk$set(dpi = 300)
  knitr::opts_chunk$set(fig.height = 4)
  knitr::opts_chunk$set(fig.width = 6.5)

  # Code will not be rendered to the Synapse wiki page
  knitr::opts_chunk$set(echo = FALSE)

  # Prints out a markdown header for the Synapse wiki and a reminder to download
  # the HTML file for better display, but only if render_for_wiki is set to TRUE.
  dataset_display <- stringr::str_replace(dataset, "_", " / ")
  title <- paste0("#! Diverse Cohorts RNA-Seq QC (", dataset_display, ")")
  cat(title, "\n\n")
  cat("**Author:** Jaclyn Beck (Sage Bionetworks)", "\n\n")
  cat("> **Note:** This wiki page shows the output of the R notebook but with",
      "code and potentially identifying information removed. **Download the",
      "rendered HTML file** with the `Download Options` button to view the",
      "notebook with sample-level details, code, and interactive PCA graphs.",
      "\n\n", "${toc}", "\n\n") # wiki table of contents
}


# ---- download-metadata ----

# Download the individual, biospecimen, and assay metadata
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

# Columbia only has a single individual with "race" = "Asian" and a single
# individual with "race" = "White", so we need to remove those individuals.
# Otherwise we can't use race as a covariate in a regression.
columbia_remove_ids <- subset(ind, dataContributionGroup == "Columbia" &
                                race %in% c("Asian", "White")) |>
  pull(individualID)

stopifnot(length(columbia_remove_ids) == 2)

# Combine the 3 metadatas and fix a few fields
metadata <- assay |>
  merge(bio) |>
  merge(ind) |>

  # Remove 7 duplicate Rush samples that are in batch B74, which are marked as "exclude"
  subset(!exclude) |>

  # Remove the single STG sample and single Asian sample from Columbia
  subset(specimenID != configs$Columbia$stg_remove_id &
           !(individualID %in% columbia_remove_ids)) |>

  # Fix or alter some fields
  mutate(
    # Make PMI numeric
    PMI = suppressWarnings(as.numeric(PMI)),

    # Make specimenID match column names of count matrix -- needs to be done
    # last so we don't break references to original specimenIDs above
    specimenID = make.names(specimenID)
  )


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

# Shorten tissue names for display but save the original tissue names
metadata <- metadata |>
  mutate(tissue_orig = tissue,
         tissue = as.character(tissue),
         tissue = case_match(tissue,
                             "caudate nucleus" ~ "CN",
                             "dorsolateral prefrontal cortex" ~ "DLPFC",
                             "superior temporal gyrus" ~ "STG",
                             "temporal pole" ~ "TP",
                             .default = tissue),
         tissue = factor(tissue))


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
  merge(dplyr::select(metadata, specimenID, tissue, tissue_orig), df)
})

# Columbia and Rush have to alter specimen IDs before merging, taken care of elsewhere
if (dataset != "Columbia" & dataset != "Rush") {
  multiqc_stats <- merge(dplyr::select(metadata, specimenID, tissue, tissue_orig),
                         multiqc_stats)
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

  multiqc_stats <- merge(dplyr::select(metadata, specimenID, tissue, tissue_orig),
                         multiqc_stats)
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
  multiqc_stats <- merge(dplyr::select(metadata, specimenID, tissue, tissue_orig),
                         multiqc_stats)
}


# ---- plot-demographic-info ----

make_bar_plot <- function(metadata, var_of_interest, facet_var = "tissue") {
  ord <- metadata |>
    mutate(variable = fct_infreq(factor(get(var_of_interest)))) |>
    pull(variable) |>
    levels()

  meta_tmp <- metadata |>
    group_by_at(c(var_of_interest, facet_var)) |>
    dplyr::count() |>
    mutate(variable = factor(get(var_of_interest), levels = ord))

  # Add a little extra space for the number label at the top of each bar
  max_y <- ceiling(max(meta_tmp$n) * 1.02)

  # Decrease bar label text size because a large number of facets makes the
  # graphs shorter. Also add even more space for the number label to account
  # for shorter graphs
  if (!is.null(facet_var) && length(unique(meta_tmp[[facet_var]])) > 3) {
    text_size = 3
    max_y <- ceiling(max_y * 1.2)
  } else if (length(unique(metadata[[var_of_interest]])) > 4) {
    text_size = 3
  } else {
    text_size = 4
  }

  plt <- ggplot(meta_tmp, aes(x = variable, y = n, fill = variable)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.5, size = text_size) +
    theme_bw() +
    xlab(NULL) +
    ylab("count") +
    ylim(0, max_y) +
    labs(fill = var_of_interest) +
    scale_fill_viridis(discrete = TRUE, begin = 0.2) +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"),
          title = element_text(face = "bold"))

  if (!is.null(facet_var)) {
    plt <- plt + facet_wrap(facet_var) +
      labs(title = paste(facet_var, "vs", var_of_interest))
  } else {
    plt <- plt + labs(title = var_of_interest)
  }

  # Shorten legend title for dataContributionGroup to save space
  if (var_of_interest == "dataContributionGroup") {
    plt <- plt + labs(fill = "Contrib. Group")
  }

  plt
}

make_bar_plot(metadata, "tissue", facet_var = NULL)

make_bar_plot(metadata, "ADoutcome")

make_bar_plot(metadata, "sex")

make_bar_plot(metadata, "race")

make_bar_plot(metadata, "isHispanic")

make_bar_plot(metadata, "cohort")

format_histogram_plot <- function(plt) {
  plt + geom_histogram(aes(fill = after_stat(count)), bins = 30, na.rm = TRUE) +
    theme_bw() +
    facet_wrap(~tissue) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          title = element_text(face = "bold"))
}

meta_age <- metadata |>
  mutate(ageDeath = case_when(
    is.na(ageDeath) | ageDeath == "missing or unknown" ~ NA,
    ageDeath == "90+" ~ 90,
    .default = suppressWarnings(as.numeric(ageDeath))
  ))

(ggplot(meta_age, aes(x = ageDeath)) + labs(title = "ageDeath")) |>
 format_histogram_plot()

(ggplot(metadata, aes(x = PMI)) + labs(title = "PMI")) |>
  format_histogram_plot()

(ggplot(metadata, aes(x = RIN)) + labs(title = "RIN")) |>
  format_histogram_plot()

(ggplot(metadata, aes(x = DV200)) + labs(title = "DV200")) |>
  format_histogram_plot()


# ---- validate-fastqc ----

thresholds <- configs$thresholds

base_content <- fastqc_data$base_content |>
  group_by(read, tissue, base) |>
  mutate(outlier_color = ifelse(is_outlier_IQR(mean_base_deviation),
                                "black", NA)) |>
  ungroup()

plt1 <- ggplot(base_content,
               aes(x = base, y = mean_base_deviation, fill = base)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(color = base_content$outlier_color,
              size = 0.5, width = 0.05, na.rm = TRUE) +
  theme_bw() +
  facet_grid(rows = vars(read), cols = vars(tissue)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, alpha = 0.8) +
  theme(strip.background = element_rect(fill = "#EEEEEE", color = "black"),
        strip.text = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  ggtitle("Mean deviation from expected base proportions")

print(plt1)

# Using the mean content deviation from 25% of each base at each position.
# Mayo uses sum of deviation instead of mean, but they're equivalent for the
# purposes of outlier finding. Unlike for the graph, we only want outliers that
# are on the upper tail (over-representation of a base).
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

if (nrow(meta_sub) > 0 && !render_for_wiki) {
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

stat_colors <- c("Pass" = "black", "Warn" = "darkorange", "Fail" = "red")
plt1 <- ggplot(mqc_plot,
               aes(x = tissue, y = samtools_reads_mapped_percent, fill = tissue)) +
  geom_boxplot(outliers = FALSE, width = 0.1) +
  geom_jitter(aes(color = mapped_status),
              width = 0.2,
              size = ifelse(mqc_plot$mapped_status == "Pass", 0.5, 1)) +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) +
  scale_color_manual(values = stat_colors) +
  scale_fill_viridis(discrete = TRUE, begin = 0.2, alpha = 0.7) +
  ggtitle("Percentage of reads mapped")

plt2 <- ggplot(mqc_plot,
               aes(x = tissue, y = picard_PERCENT_DUPLICATION, fill = tissue)) +
  geom_boxplot(outliers = FALSE, width = 0.1) +
  geom_jitter(aes(color = duplicated_status),
              width = 0.2,
              size = ifelse(mqc_plot$duplicated_status == "Pass", 0.5, 1)) +
  theme_bw() +
  theme(title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) +
  labs(color = "Status") +
  scale_color_manual(values = stat_colors) +
  scale_fill_viridis(discrete = TRUE, begin = 0.2, alpha = 0.7) +
  ggtitle("Percentage of reads duplicated") +
  theme(plot.margin = unit(c(0, 0, 0, 50), units = "pt"))

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

if (nrow(meta_sub) > 0 && !render_for_wiki) {
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

plts[[1]] <- plts[[1]] + scale_color_viridis(discrete = TRUE,
                                             begin = 0.2, end = 0.8)

# Make the valid points a little lighter so the mismatches stand out more
if (length(mismatches$mismatches) > 0) {
  plts[[2]] <- plts[[2]] +
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "gray"))
}

print(plts[[1]] + plts[[2]])

metadata$sex_valid <- !(metadata$specimenID %in% mismatches$mismatches)


# ---- print-sex-mismatches ----

meta_sub <- subset(metadata, !sex_valid) |>
  select(specimenID, tissue) |>
  dplyr::rename(`Specimen ID` = specimenID,
                Tissue = tissue) |>
  arrange(Tissue, `Specimen ID`)

if (nrow(meta_sub) > 0 && !render_for_wiki) {
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
  pca_df <- results$group_results[[res_name]]$pca_df
  pc_thresholds <- results$group_results[[res_name]]$thresholds

  # Color points red if "is_outlier" is TRUE
  stat_colors <- c("FALSE" = "black", "TRUE" = "red")

  plt <- plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3,
                 color = ~is_outlier, colors = stat_colors,
                 mode = "markers", type = "scatter3d",
                 size = 0.5,
                 legendgrouptitle = list(text = "Outlier", font = list(size = 16))) |>
    plotly::layout(title = res_name)

  return(plt)
})

if (render_for_wiki) {
  cat("Download HTML notebook for PCA plots.", "\n")
} else {
  htmltools::tagList(plts)
}

metadata$pca_valid <- !(metadata$specimenID %in% results$outliers)


# ---- print-pca-outliers ----

meta_sub <- subset(metadata, !pca_valid) |>
  select(specimenID, tissue) |>
  dplyr::rename(`Specimen ID` = specimenID, Tissue = tissue) |>
  arrange(Tissue, `Specimen ID`)

if (nrow(meta_sub) > 0 && !render_for_wiki) {
  meta_sub
}


# ---- validate-dv200 ----

thresholds <- configs$thresholds

plt1 <- ggplot(metadata, aes(x = tissue, y = RIN, fill = tissue)) +
  geom_boxplot(outliers = FALSE, na.rm = TRUE, width = 0.5) +
  geom_jitter(size = 0.5, na.rm = TRUE, width = 0.3) +
  xlab(NULL) +
  theme_bw() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, alpha = 0.7) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"))

plt2 <- ggplot(metadata, aes(x = tissue, y = DV200, fill = tissue)) +
  geom_boxplot(outliers = FALSE, na.rm = TRUE, width = 0.5) +
  geom_jitter(size = 0.5, na.rm = TRUE, width = 0.3) +
  geom_hline(yintercept = thresholds$DV200, linetype = "dotdash", color = "darkgray") +
  xlab(NULL) +
  theme_bw() +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 20), units = "pt"))

plt3 <- ggplot(metadata, aes(x = RIN, y = DV200, color = tissue)) +
  geom_jitter(size = 0.5, na.rm = TRUE) +
  geom_hline(yintercept = thresholds$DV200, linetype = "dotdash", color = "darkgray") +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"))

dv200_rank <- metadata |>
  group_by(tissue) |>
  mutate(rank_DV200 = rank(DV200))

plt4 <- ggplot(dv200_rank, aes(x = rank_DV200, y = DV200, color = tissue)) +
  geom_point(size = 0.5, na.rm = TRUE) +
  geom_hline(yintercept = thresholds$DV200, linetype = "dotdash", color = "darkgray") +
  theme_bw() +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(face = "bold"),
        plot.margin = unit(c(0, 0, 0, 20), units = "pt"))

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

if (nrow(meta_sub) > 0 && !render_for_wiki) {
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
fastqc_basic_stats <- subset(fastqc_data$basic_statistics,
                             specimenID %in% metadata$specimenID)

# Undo tissue renames
metadata <- mutate(metadata, tissue = tissue_orig) |> select(-tissue_orig)
multiqc_stats <- mutate(multiqc_stats, tissue = tissue_orig) |> select(-tissue_orig)
fastqc_basic_stats <- mutate(fastqc_basic_stats, tissue = tissue_orig) |> select(-tissue_orig)

# Remove genes that are all 0's
zeros <- rowSums(counts) == 0

data_final <- list("metadata" = metadata, "counts" = counts[!zeros, ],
                   "multiqc_stats" = multiqc_stats,
                   "fastqc_stats" = fastqc_basic_stats)

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
  dplyr::rename(Tissue = Var1) |>
  as.data.frame()


# ---- print-qc-failures-summary ----

if (render_for_wiki) {
  cat("Download HTML notebook for sample-level details.", "\n")
} else {
  failures |>
    group_by(tissue) |>
    summarize(`Specimen IDs` = paste(str_replace(specimenID, "^X", ""),
                                     collapse = ", "))
}


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

if (render_for_wiki) {
  cat("Download HTML notebook for sample-level details.", "\n")
} else {
  failures_detail
}
