source("helper_functions.R")

configs <- config::get(file = "config.yml")

metadata <- download_metadata(configs)
counts <- download_rsem(configs$Rush$count_matrix_synid)

metadata <- subset(metadata, specimenID %in% colnames(counts))
counts <- counts[, metadata$specimenID]

stopifnot(length(unique(metadata$specimenID)) == nrow(metadata))

fastqc_data <- readRDS(file.path("data", "QC", "Rush_fastqc_stats.rds"))
multiqc_stats <- readRDS(file.path("data", "QC", "Rush_multiqc_stats.rds"))
gene_info <- read.csv(file.path("data", "gene_metadata.csv"))

# Remove the duplicate samples from Rush
to_remove <- lapply(configs$Rush$remove_samples_fastqc, function(id) {
  grep(id, fastqc_data$basic_statistics$sample, value = TRUE)
})

fastqc_data <- lapply(fastqc_data, function(df) {
  df <- subset(df, !(sample %in% unlist(to_remove)))
  merge(dplyr::select(metadata, specimenID, tissue), df)
})

multiqc_stats <- subset(multiqc_stats, !(specimenID %in% configs$Rush$remove_samples_fastqc)) |>
  mutate(specimenID = str_replace(specimenID, "_S[0-9]+", ""))
multiqc_stats <- merge(dplyr::select(metadata, specimenID, tissue), multiqc_stats)


# Validate ---------------------------------------------------------------------

orig_size <- ncol(counts)

counts_log <- simple_log2norm(counts)

metadata <- validate_fastqc(metadata, fastqc_data, configs$thresholds)
metadata <- validate_multiqc(metadata, multiqc_stats, configs$thresholds)
metadata <- validate_sex(metadata, counts_log, configs$thresholds)
metadata <- validate_pca(metadata, counts_log, gene_info)
metadata <- validate_DV200(metadata, configs$thresholds)


# Save samples that passed QC --------------------------------------------------

metadata$valid <- Reduce("&", metadata[, grepl("_valid", colnames(metadata))])
metadata$warn <- Reduce("+", metadata[, grepl("_warn", colnames(metadata))])

metadata$valid <- metadata$valid & metadata$warn < 2

print(table(metadata$tissue, metadata$valid))

metadata <- subset(metadata, valid == TRUE)
counts <- counts[, metadata$specimenID]
multiqc_stats <- subset(multiqc_stats, specimenID %in% metadata$specimenID)

message(str_glue("{ncol(counts)} of {orig_size} samples passed QC."))

# Remove genes that are all 0's
zeros <- rowSums(counts) == 0

data_final <- list("metadata" = metadata, "counts" = counts[!zeros, ],
                   "multiqc_stats" = multiqc_stats)

saveRDS(data_final, file.path("data", "QC", "Rush_qc.rds"))
