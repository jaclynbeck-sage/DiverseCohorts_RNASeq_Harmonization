library(synapser)
library(ggplot2)
library(viridis)
library(patchwork)
library(matrixStats)
library(dplyr)
library(stringr)
library(sageRNAUtils)

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
      RINbin = case_when(
        RIN <= 2.5 ~ "0 to 2.5",
        RIN > 2.5 & RIN <= 5 ~ "2.6 to 5",
        RIN > 5 & RIN <= 7.5 ~ "5.1 to 7.5",
        RIN > 7.5 ~ "7.6+"
      ),
      apoe4Status = case_when(
        grepl("4", apoeGenotype) ~ "True",
        apoeGenotype == "missing or unknown" ~ "missing or unknown",
        .default = "False"
      ),
      ageDeathNumeric = suppressWarnings(as.numeric(ageDeath)),
      pmiNumeric = suppressWarnings(as.numeric(PMI)),
      ageBin = case_when(
        ageDeathNumeric < 65 ~ "Under 65",
        ageDeathNumeric >= 65 & ageDeathNumeric < 70 ~ "65 to 69",
        ageDeathNumeric >= 70 & ageDeathNumeric < 75 ~ "70 to 74",
        ageDeathNumeric >= 75 & ageDeathNumeric < 80 ~ "75 to 79",
        ageDeathNumeric >= 80 & ageDeathNumeric < 85 ~ "80 to 74",
        ageDeathNumeric >= 85 & ageDeathNumeric < 90 ~ "85 to 90",
        .default = ageDeath
      ),
      pmiBin = case_when(
        pmiNumeric < 10 ~ "Under 10",
        pmiNumeric >= 10 & pmiNumeric < 20 ~ "10 to 19",
        pmiNumeric >= 20 & pmiNumeric < 30 ~ "20 to 29",
        pmiNumeric >= 30 ~ "30+",
        .default = PMI
      ),
      # Sequencing batches need to be re-named to be unique to each data set --
      # Mayo sequenced "Mayo", "Emory", and sample swaps in the same batches
      # Rush did the same with "Rush" and sample swaps
      # NYGC sequenced "Columbia" separately from "MSSM" -- TODO sample swaps?
      sequencingBatch = case_match(dataGenerationSite,
        c("Mayo", "Rush") ~ paste(dataGenerationSite, sequencingBatch, sep = "_"),
        "NYGC" ~ paste(dataGenerationSite, dataContributionGroup, sequencingBatch, sep = "_")
      )
    )
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


# tail can be "upper", "lower", or "both"
# TODO move to library
is_outlier_IQR <- function(data, tail = "both") {
  iqr <- IQR(data, na.rm = TRUE)
  q1 <- quantile(data, 0.25, na.rm = TRUE)
  q3 <- quantile(data, 0.75, na.rm = TRUE)

  switch(
    tail,
    "both" = (data < q1 - 1.5 * iqr) | (data > q3 + 1.5 * iqr),
    "upper" = data > q3 + 1.5 * iqr,
    "lower" = data < q1 - 1.5 * iqr
  ) | is.na(data)
}

is_outlier_SD <- function(data, n_sds = 4, tail = "both") {
  mean_d <- mean(data, na.rm = TRUE)
  stdev <- n_sds * sd(data, na.rm = TRUE)

  switch(
    tail,
    "both" = (data < mean_d - stdev) | (data > mean_d + stdev),
    "upper" = data > mean_d + stdev,
    "lower" = data < mean_d - stdev
  ) | is.na(data)
}


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


# In the absence of RSeQC data, we only check percentage of reads mapped for
# now. Could consider duplicated reads, rsem uniquely aligned, and cutadapt
# trimming, all of which have clear outliers.
validate_multiqc <- function(metadata, multiqc_stats, thresholds) {
  plt1 <- ggplot(multiqc_stats,
                 aes(x = tissue, y = samtools_reads_mapped_percent, fill = tissue)) +
    geom_boxplot() +
    theme_bw() +
    ggtitle("Percentage of reads mapped")

  print(plt1)

  corr_mat <- multiqc_stats |>
    dplyr::select(where(is.numeric), -samtools_reads_QC_failed_percent) |>
    cor()

  corrplot::corrplot(corr_mat)

  reads_mapped_fail <- multiqc_stats |>
    subset(samtools_reads_mapped_percent < thresholds$reads_mapped) |>
    pull(specimenID)

  reads_mapped_outliers <- multiqc_stats |>
    mutate(is_outlier = is_outlier_IQR(samtools_reads_mapped_percent, tail = "lower")) |>
    subset(is_outlier == TRUE) |>
    pull(specimenID)

  metadata$reads_mapped_valid <- !(metadata$specimenID %in% reads_mapped_fail)
  metadata$reads_mapped_warn <- metadata$specimenID %in% reads_mapped_outliers

  return(metadata)
}


validate_sex <- function(metadata, data, thresholds) {
  # source: https://www.ncbi.nlm.nih.gov/pubmed/23829492
  y_genes <- c(RPS4Y1 = "ENSG00000129824",
               EIF1AY = "ENSG00000198692",
               DDX3Y = "ENSG00000067048",
               KDM5D = "ENSG00000012817")
  sex_genes <- c(XIST = "ENSG00000229807", y_genes)

  # Strip any version numbers off of the gene names
  rownames(data) <- str_replace(rownames(data), "\\.[0-9]+", "")

  sex_check <- metadata %>%
    dplyr::select(specimenID, sex) %>%
    merge(
      t(data[sex_genes, ]),
      by.x = "specimenID", by.y = "row.names"
    ) %>%
    dplyr::rename(XIST = ENSG00000229807) %>%
    mutate(
      mean_Y = rowMeans(dplyr::select(., all_of(y_genes))),
      est_sex = case_when(mean_Y > thresholds$sex ~ "male",
                          .default = "female"),
      sex_valid = sex == est_sex
    ) %>%
    group_by(sex) %>%
    mutate(
      sex_warn = sex_valid & is_outlier_SD(mean_Y),
      status = case_when(sex_valid & !sex_warn ~ "Pass",
                         sex_warn ~ "Warning",
                         .default = "Fail")
    ) %>%
    ungroup() %>%
    # Put FALSE last so they are plotted on top of other dots
    arrange(desc(sex_valid))

  plt1 <- ggplot(sex_check, aes(x = sex, y = XIST, fill = sex)) +
    geom_boxplot(width = 0.25, outliers = FALSE) +
    geom_jitter(size = 0.5) +
    theme_bw() +
    theme(legend.position = "bottom")

  plt2 <- ggplot(sex_check, aes(x = sex, y = mean_Y, fill = sex)) +
    geom_boxplot(width = 0.25, outliers = FALSE) +
    geom_jitter(size = 0.5) +
    theme_bw() +
    theme(legend.position = "bottom")

  plt3 <- ggplot(sex_check, aes(x = XIST, y = mean_Y, color = sex)) +
    geom_point(size = 0.5) +
    theme_bw() +
    theme(legend.position = "bottom")

  status_colors <- c("Pass" = "black", "Warning" = "orange", "Fail" = "red")
  plt4 <- ggplot(sex_check, aes(x = XIST, y = mean_Y, color = status)) +
    geom_point(size = ifelse(sex_check$sex_valid & !sex_check$sex_warn, 0.5, 1)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = status_colors) +
    geom_hline(yintercept = thresholds$sex, linetype = "dotdash",
               color = "red", alpha = 0.8)

  print(plt1 + plt2 + plt3 + plt4)

  rownames(sex_check) <- sex_check$specimenID
  sex_check <- sex_check[metadata$specimenID, ]

  metadata$sex_valid <- sex_check$sex_valid
  metadata$sex_warn <- sex_check$sex_warn

  metadata
}


outlier_pca <- function(metadata, data, gene_info, n_sds = 4) {
  metadata$pca_valid <- FALSE

  # Use protein-coding autosomal genes, per Mayo's processing pipeline
  protein_coding <- subset(gene_info, gene_biotype == "protein_coding" &
                             grepl("chr[1-9]+", chromosome_name)) |>
    pull(ensembl_gene_id)

  pca_dfs <- list()
  ellipse_dfs <- list()

  # Formula for an ellipse is (x^2 / a^2) + (y^2 / b^2) = 1, so
  # y = +/- sqrt((1 - x^2 / a^2) * b^2)
  ellipse_points <- function(axis_A, axis_B) {
    x <- seq(from = -axis_A, to = axis_A, length.out = 1000)
    y <- sqrt((1 - x^2 / axis_A^2) * axis_B^2)

    data.frame(x = c(x, rev(x)), y = c(y, rev(-y)))
  }

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

  for (grp in unique(metadata$pca_group)) {
    meta_group <- subset(metadata, pca_group == grp)
    data_group <- data[protein_coding, meta_group$specimenID]

    # Remove genes that are mostly 0's, which are -1 in log2-scale. For PCA,
    # restrict to genes expressed in >= 80% of samples
    genes_keep <- rowSums(data_group > -1) >= 0.8 * ncol(data_group)
    data_group <- data_group[genes_keep, ]

    # Remove genes with very low variance in this group, per Mayo's processing pipeline
    variance <- rowVars(data_group)
    genes_keep <- names(variance)[variance > 0.001]
    data_group <- data_group[genes_keep, ]

    pc_group <- prcomp(t(data_group), center = TRUE, scale. = TRUE)
    pc_df <- merge(pc_group$x, meta_group,
                   by.x = "row.names", by.y = "specimenID") |>
      dplyr::rename(specimenID = Row.names)

    pc1_thresh <- sd(pc_df$PC1) * n_sds
    pc2_thresh <- sd(pc_df$PC2) * n_sds

    in_ellipse <- (pc_df$PC1^2 / pc1_thresh^2) + (pc_df$PC2^2 / pc2_thresh^2)

    pc_df$pca_valid <- in_ellipse < 1

    pca_dfs[[grp]] <- pc_df |>
      dplyr::select(specimenID, PC1, PC2, tissue, sequencingBatch, pca_group, pca_valid)

    ellipse_dfs[[grp]] <- ellipse_points(pc1_thresh, pc2_thresh) |>
      mutate(
        tissue = unique(meta_group$tissue),
        pca_group = unique(meta_group$pca_group),
        # This variable is only needed for Rush
        sequencingBatch = ifelse(
          unique(meta_group$dataGenerationSite == "Rush"),
          unique(meta_group$sequencingBatch),
          1
        )
      )
  }

  pca_dfs <- do.call(rbind, pca_dfs)
  ellipse_dfs <- do.call(rbind, ellipse_dfs)

  plt1 <- ggplot(pca_dfs, aes(x = PC1, y = PC2, color = sequencingBatch)) +
    geom_point(size = 0.8) +
    geom_path(data = ellipse_dfs,
              aes(x = x, y = y),
              linetype = "dotdash",
              color = "blue") +
    theme_bw() +
    theme(legend.position = "bottom")

  plt2 <- ggplot(pca_dfs, aes(x = PC1, y = PC2, color = pca_valid)) +
    geom_point(size = ifelse(pca_dfs$pca_valid, 0.5, 1)) +
    geom_path(data = ellipse_dfs,
              aes(x = x, y = y),
              linetype = "dotdash",
              color = "blue") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red"))

  if (unique(metadata$dataGenerationSite == "Rush")) {
    n_columns <- length(unique(ellipse_dfs$pca_group))
    plt_final <- (plt1 + facet_wrap(~tissue + sequencingBatch, ncol = n_columns)) /
      (plt2 + facet_wrap(~tissue + sequencingBatch, ncol = n_columns))
    print(plt_final)
  } else {
    print((plt1 + facet_wrap(~tissue)) / (plt2 + facet_wrap(~tissue)))
  }

  pca_dfs <- subset(pca_dfs, pca_valid)
  metadata$pca_valid[metadata$specimenID %in% pca_dfs$specimenID] <- TRUE

  metadata
}


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
