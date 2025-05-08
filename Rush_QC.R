library(edgeR)
source("helper_functions.R")

metadata <- download_metadata()
counts <- download_rsem("syn64176441")
fastq_stats <- download_fastq("syn66358698", load_saved_stats = TRUE)
fastq_summary <- fastq_stats$fastq_summary
gc_distribution <- fastq_stats$gc_distribution

gene_info <- read.csv(file.path("data", "gene_lengths_gc.csv"))

metadata <- subset(metadata, specimenID %in% colnames(counts))
counts <- counts[, metadata$specimenID]

# Remove the duplicate samples from Rush -- there will be 4 files instead of 2
# for these specimen IDs. The files to remove are all labeled as "<id>_S3xx"
dupes <- which(table(fastq_summary$specimenID) == 4)
to_remove <- lapply(names(dupes), function(id) {
  grep(paste0(id, "_S3[0-9]+"), fastq_summary$sample, value = TRUE)
})

fastq_summary <- subset(fastq_summary,
                        specimenID %in% metadata$specimenID &
                          !(sample %in% unlist(to_remove)))
gc_distribution <- subset(gc_distribution,
                          sample %in% fastq_summary$sample)

stopifnot(all(duplicated(metadata$specimenID) == FALSE))

orig_size <- ncol(counts)

counts_log <- simple_lognorm(counts)

metadata <- validate_sex(metadata, counts_log)
metadata <- outlier_pca(metadata, counts_log, gene_info)

metadata$lib_size <- colSums(counts)

ggplot(metadata, aes(x = tissue, y = lib_size, fill = tissue)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 0.5) +
  theme_bw()

metadata <- validate_DV200(metadata)


# Save samples that passed QC

metadata$valid <- metadata$sex_valid & metadata$pca_valid & metadata$dv200_valid
print(table(metadata$tissue, metadata$valid))

metadata <- subset(metadata, valid == TRUE)
counts <- counts[, metadata$specimenID]

message(str_glue("{ncol(counts)} of {orig_size} samples passed QC."))

data_final <- DGEList(counts, samples = metadata, remove.zeros = TRUE)
data_final <- normLibSizes(data_final, method = "TMM")

saveRDS(data_final, file.path("data", "QC", "Rush_qc.rds"))
