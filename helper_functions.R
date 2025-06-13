library(synapser)
library(ggplot2)
library(viridis)
library(patchwork)
library(matrixStats)
library(dplyr)
library(stringr)
library(sageRNAUtils)

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
