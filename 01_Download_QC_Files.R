# This script downloads QC-related files for each study from Synapse and
# pre-processes them so they can be used in QC. This will take some time as
# there are separate fastqc files for each sample, and downloading and parsing
# each one can be slow.
#
# Note that if the fastqc download has been done before, the raw data is already
# saved and won't need to be re-downloaded or re-parsed.
source("helper_functions.R")

QC_folder <- file.path("data", "QC")
configs <- config::get(file = "config.yml")

# FastQC output ----------------------------------------------------------------

columbia_fastqc <- download_fastqc(configs$Columbia, load_saved_stats = FALSE)
saveRDS(columbia_fastqc, file.path(QC_folder, "Columbia_fastqc_stats.rds"))

mayo_emory_fastqc <- download_fastqc(configs$Mayo_Emory, load_saved_stats = FALSE)
saveRDS(mayo_emory_fastqc, file.path(QC_folder, "Mayo_Emory_fastqc_stats.rds"))

mssm_fastqc <- download_fastqc(configs$MSSM, load_saved_stats = FALSE)
saveRDS(mssm_fastqc, file.path(QC_folder, "MSSM_fastqc_stats.rds"))

rush_fastqc <- download_fastqc(configs$Rush, load_saved_stats = FALSE)
saveRDS(rush_fastqc, file.path(QC_folder, "Rush_fastqc_stats.rds"))


# MultiQC output ---------------------------------------------------------------

columbia_multiqc <- download_multiqc_json(configs$Columbia$multiqc_json_synids)
saveRDS(columbia_multiqc, file.path(QC_folder, "Columbia_multiqc_stats.rds"))

mayo_emory_multiqc <- download_multiqc_json(configs$Mayo_Emory$multiqc_json_synids)
saveRDS(mayo_emory_multiqc, file.path(QC_folder, "Mayo_Emory_multiqc_stats.rds"))

mssm_multiqc <- download_multiqc_json(configs$MSSM$multiqc_json_synids)
saveRDS(mssm_multiqc, file.path(QC_folder, "MSSM_multiqc_stats.rds"))

rush_multiqc <- download_multiqc_json(configs$Rush$multiqc_json_synids)
saveRDS(rush_multiqc, file.path(QC_folder, "Rush_multiqc_stats.rds"))


# Gene length and GC content ---------------------------------------------------

gc_content <- sageRNAUtils::get_gc_content_gtf(configs$gtf_file,
                                               configs$fasta_file,
                                               include_introns = FALSE)

write.csv(gc_content, file.path("data", "gene_metadata.csv"),
          row.names = FALSE, quote = FALSE)

# TODO upload to Synapse
