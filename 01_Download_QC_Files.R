# This script downloads QC-related files for each study from Synapse and
# pre-processes them so they can be used in QC. This will take some time as
# there are separate fastqc files for each sample, and downloading and parsing
# each one can be slow.
#
# Note that if the fastqc download has been done before, the raw data is already
# saved and won't need to be re-downloaded or re-parsed.
source("helper_functions.R")

QC_folder <- file.path("data", "QC")

# FastQC output ----------------------------------------------------------------

columbia_fastqc <- download_fastq("Columbia", "syn66353726", load_saved_stats = TRUE)
saveRDS(columbia_fastqc, file.path(QC_folder, "Columbia_fastqc_stats.rds"))

mayo_emory_fastqc <- download_fastq("Mayo_Emory", "syn66354409", load_saved_stats = TRUE)
saveRDS(mayo_emory_fastqc, file.path(QC_folder, "Mayo_Emory_fastqc_stats.rds"))

mssm_fastqc <- download_fastq("MSSM", "syn66354135", load_saved_stats = TRUE)
saveRDS(mssm_fastqc, file.path(QC_folder, "MSSM_fastqc_stats.rds"))

rush_fastqc <- download_fastq("Rush", "syn66358698", load_saved_stats = TRUE)
saveRDS(rush_fastqc, file.path(QC_folder, "Rush_fastqc_stats.rds"))


# MultiQC output ---------------------------------------------------------------

columbia_multiqc <- download_multiqc_json("syn64558076")
saveRDS(columbia_multiqc, file.path(QC_folder, "Columbia_multiqc_stats.rds"))

mayo_emory_multiqc <- download_multiqc_json("syn64558505")
saveRDS(mayo_emory_multiqc, file.path(QC_folder, "Mayo_Emory_multiqc_stats.rds"))

mssm_multiqc <- download_multiqc_json("syn64558496")
saveRDS(mssm_multiqc, file.path(QC_folder, "MSSM_multiqc_stats.rds"))

rush_multiqc <- download_multiqc_json("syn64558758")
saveRDS(rush_multiqc, file.path(QC_folder, "Rush_multiqc_stats.rds"))
