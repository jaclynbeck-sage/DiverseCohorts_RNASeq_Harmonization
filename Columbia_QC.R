source("helper_functions.R")

configs <- config::get(file = "config.yml")

# TODO multiqc stats still has old IDs

# TODO counts matrix is missing 8 samples, 7 of which are from Biggs Institute
# and 1 that is listed as a sample exchange but from Columbia
metadata <- download_metadata(configs)
counts <- download_rsem(configs$Columbia$count_matrix_synid)

fastqc_data <- readRDS(file.path("data", "QC", "Columbia_fastqc_stats.rds"))
multiqc_stats <- readRDS(file.path("data", "QC", "Columbia_multiqc_stats.rds"))
gene_info <- read.csv(file.path("data", "gene_metadata.csv"))

metadata <- subset(metadata, specimenID %in% colnames(counts))
counts <- counts[, metadata$specimenID]

stopifnot(length(unique(metadata$specimenID)) == nrow(metadata))

fastqc_data <- lapply(fastqc_data, function(df) {
  merge(dplyr::select(metadata, specimenID, tissue), df)
})

remapping <- fastqc_data$basic_statistics |>
  mutate(old_id = make.names(str_replace(Filename, "_(1|2).gz", ""))) |>
  select(specimenID, old_id) |>
  dplyr::rename(new_id = specimenID) |>
  distinct()

multiqc_stats <- merge(multiqc_stats, remapping,
                       by.x = "specimenID", by.y = "old_id") |>
  mutate(specimenID = new_id) |>
  dplyr::select(-new_id)

multiqc_stats <- merge(dplyr::select(metadata, specimenID, tissue), multiqc_stats)


# Validate ---------------------------------------------------------------------

orig_size <- ncol(counts)

counts_log <- simple_lognorm(counts)

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

message(str_glue("{ncol(counts)} of {orig_size} samples passed QC."))

# Remove genes that are all 0's
zeros <- rowSums(counts) == 0

data_final <- list("metadata" = metadata, "counts" = counts[!zeros, ],
                   "multiqc_stats" = multiqc_stats)

saveRDS(data_final, file.path("data", "QC", "Columbia_qc.rds"))
