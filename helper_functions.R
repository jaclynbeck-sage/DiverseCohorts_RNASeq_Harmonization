library(synapser)
library(ggplot2)
#library(ggrepel)
library(viridis)
library(patchwork)
library(matrixStats)
library(dplyr)
library(stringr)
library(biomaRt)
library(GenomicRanges)
library(BSgenome)

download_metadata <- function() {
  synLogin()

  ind <- synGet("syn51757646",
                downloadLocation = "downloads",
                ifcollision = "overwrite.local")$path |>
    read.csv()
  bio <- synGet("syn51757645",
                downloadLocation = "downloads",
                ifcollision = "overwrite.local")$path |>
    read.csv()
  assay <- synGet("syn51757643",
                  downloadLocation = "downloads",
                  ifcollision = "overwrite.local")$path |>
    read.csv()

  duplicates_remove <- c("Div_487", "Div_498", "Div_500", "Div_503", "Div_508",
                         "Div_509", "Div_545")

  metadata <- assay |>
    merge(bio) |>
    merge(ind) |>
    mutate(specimenID = make.names(specimenID)) |>
    # Remove 7 duplicate Rush samples that are in batch B74
    subset(!(specimenID %in% c(duplicates_remove) & sequencingBatch == "B74")) |>
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
    # Remove Ensembl version from Ensembl IDs
    mutate(gene_id = str_replace(gene_id, "\\.[0-9]+", "")) |>
    tibble::column_to_rownames("gene_id")
}


lognorm <- function(data) {
  data <- sweep(data, 2, colSums(data), "/") * 1e6
  log2(data + 0.5) |>
    as.matrix()
}


validate_sex <- function(metadata, data, sex_thresh = 2.0) {
  # source: https://www.ncbi.nlm.nih.gov/pubmed/23829492
  y_genes <- c(RPS4Y1 = "ENSG00000129824",
               EIF1AY = "ENSG00000198692",
               DDX3Y = "ENSG00000067048",
               KDM5D = "ENSG00000012817")
  sex_genes <- c(XIST = "ENSG00000229807", y_genes)

  sex_check <- metadata %>%
    dplyr::select(specimenID, sex) %>%
    merge(
      t(data[sex_genes, ]),
      by.x = "specimenID", by.y = "row.names"
    ) %>%
    dplyr::rename(XIST = ENSG00000229807) %>%
    mutate(
      mean_Y = rowMeans(dplyr::select(., all_of(y_genes))),
      est_sex = case_when(mean_Y > sex_thresh ~ "male",
                          .default = "female"),
      sex_valid = sex == est_sex
    ) %>%
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

  plt4 <- ggplot(sex_check, aes(x = XIST, y = mean_Y, color = sex_valid)) +
    geom_point(size = ifelse(sex_check$sex_valid, 0.5, 1)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) +
    geom_hline(yintercept = sex_thresh, linetype = "dotdash", color = "blue")
    #geom_label_repel(aes(label = ifelse(!sex_valid, specimenID, NA)),
    #                 segment.color = "darkgray",
    #                 na.rm = TRUE)

  print(plt1 + plt2 + plt3 + plt4)

  rownames(sex_check) <- sex_check$specimenID
  sex_check <- sex_check[metadata$specimenID, ]

  metadata$sex_valid <- sex_check$sex_valid

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


validate_DV200 <- function(metadata) {
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
    geom_hline(yintercept = 60) +
    theme_bw()

  print((plt1 + plt2) / (plt3 + plt4))


  # TODO arbitrary threshold, need to look at other data sources. 80 is where
  # the decrease starts.
  metadata$dv200_valid <- metadata$DV200 > 60

  metadata
}

get_gc_content_biomart <- function(gene_list) {
  mart <- biomaRt::useEnsembl(biomart = "ensembl",
                              dataset = "hsapiens_gene_ensembl",
                              version = 109)

  # Note: It is possible to query attribute "percentage_gene_gc_content",
  # however the percent returned by Biomart is for the entire gene sequence,
  # including introns. We want the GC content of exons only.
  gene_info <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id",
      "exon_chrom_start",  "exon_chrom_end", "strand", "start_position",
      "end_position", "external_gene_name", "chromosome_name", "gene_biotype",
      "transcript_biotype"
    ),
    filters = "ensembl_gene_id",
    values = gene_list,
    mart = mart
  )

  gene_info$hgnc_symbol <- gene_info$external_gene_name

  # Fix strand to be + and - to work with GenomicRanges
  gene_info$strand <- as.character(gene_info$strand) |>
    str_replace("-1", "-") |>
    str_replace("1", "+")

  gene_info <- gene_info |>
    # "start_position" is gene start. If the gene is on the reverse strand,
    # coordinates need to be reversed since the sequence returned by Biomart
    # is the reverse complement instead of the forward strand.
    mutate(start = case_match(strand,
                              "+" ~ exon_chrom_start - start_position + 1,
                              "-" ~ end_position - exon_chrom_end + 1),
           end = case_match(strand,
                            "+" ~ exon_chrom_end - start_position + 1,
                            "-" ~ end_position - exon_chrom_start + 1))

  seqs <- getSequence(id = unique(gene_info$ensembl_gene_id),
                      type = "ensembl_gene_id",
                      seqType = "gene_exon_intron",
                      mart = mart)

  seqs2 <- DNAStringSet(seqs$gene_exon_intron)
  names(seqs2) <- seqs$ensembl_gene_id

  grange <- GenomicRanges::makeGRangesFromDataFrame(
    gene_info,
    keep.extra.columns = TRUE,
    seqnames.field = "ensembl_gene_id"
  )
  grange$ensembl_gene_id <- as.character(seqnames(grange))

  message("Calculating length and GC content for each gene...")
  gene_data <- .calculate_gc_content(grange, seqs2)

  message("Calculating length and GC content for each transcript...")
  transcript_data <- .calculate_gc_content(grange, seqs2, by_transcript = TRUE)

  return(list("genes" = gene_data, "transcripts" = transcript_data))
}


# Alternate to biomart ---------------------------------------------------------

# Uses original GTF and FASTA files as input instead of querying Biomart. This
# is the preferred method since Biomart doesn't keep annotations on all genes
# in the reference and it occasionally doesn't respond to requests.

# GTF file: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz
# FASTA file: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
get_gc_content_gtf <- function(gtf_file, fasta_file) {
  message("Importing GTF and FASTA data...")
  # Import GTF and FASTA once, use for both gene-level and transcript-level calculations
  gtf <- rtracklayer::import(gtf_file,
                             format = "gtf",
                             feature.type = "exon")

  # Some renames to match Biomart
  gtf$ensembl_gene_id <- gtf$gene_id
  gtf$transcript_id <- gtf$transcript_id
  gtf$chromosome_name <- as.character(seqnames(gtf))
  gtf$hgnc_symbol <- gtf$gene_name
  gtf$gene_biotype <- gtf$gene_type
  gtf$transcript_biotype <- gtf$transcript_type

  fasta <- rtracklayer::import(fasta_file,
                               format = "fasta",
                               type = "DNA")
  # Change chromosome names from things like "chr2 2" to "chr2" to match
  # GTF file
  names(fasta) <- str_replace(names(fasta), " .*", "")

  message("Calculating length and GC content for each gene...")
  gene_data <- .calculate_gc_content(gtf, fasta)

  message("Calculating length and GC content for each transcript...")
  transcript_data <- .calculate_gc_content(gtf, fasta, by_transcript = TRUE)

  return(list("genes" = gene_data, "transcripts" = transcript_data))
}


# Helper function to calculate GC content and gene [or transcript] length.
#
# Arguments:
#   ranges_obj - A `GRanges` object containing the start and end coordinates of
#     each exon for each gene. Overlapping exon coordinates should be merged
#     with `reduce()` first. This object must contain either an
#     `ensembl_transcript_id` and/or `ensembl_gene_id` metadata column,
#     depending on whether GC content/length are being calculated for
#     transcripts or for genes. It should also have an `hgnc_symbol` column.
#     Values in `seqnames(ranges_obj)` must match the values in
#     `names(sequences)`.
#   sequences - a `DNAStringSet` object where each item is a `DNAString`
#     containing a sequence to subset exon regions from. `names(sequences)`
#     must match the values in `seqnames(ranges_obj)`.
#
# If the original input was imported from GTF and FASTA files, each `DNAString`
# in `sequences` should be the sequence of an entire chromosome, and the
# coordinates in `ranges_obj` should be relative to the start of the chromosome,
# which is the default in GTF files.
#
# If the original input was from Biomart, each `DNAString` in `sequences` should
# be the sequence of a single gene, and the coordinates in `ranges_obj` should
# be relative to the start of the gene. Note that for genes on the reverse
# strand, Biomart returns the gene sequence as the reverse complement of the
# reference strand, so that position 1 on the reference strand is position
# <gene_length> on the reverse strand, and position <gene_length> of the
# reference strand is position 1 of the reverse strand. However, the coordinates
# given by Biomart are relative to the reference sequence. This means that one
# of two adjustments need to be made *prior* to calling this function.
# Either:
#   1. `sequences` needs to be updated so that the sequence for every gene on
#      the negative/reverse strand is reverse-complemented,
# or
#   2. coordinates for exons on the reverse strand need to be adjusted to be
#      relative to the start of reverse complement sequence, e.g. if Biomart
#      returns coordinates `start = 3` and `end = 10`, the start position is
#      now the *end* coordinate relative to the reverse complement strand and is
#      adjusted as position <gene_length - 3 + 1>. The original end position is
#      now the start coordinate and is adjusted as <gene_length - 10 + 1>.
#
# No such adjustment is necessary for data from GTF and FASTA files, as the
# sequences are always the reference strand.
.calculate_gc_content <- function(ranges_obj, sequences, by_transcript = FALSE) {
  # Genes will be a GRangesList where each list item is a GRanges object for a
  # specific gene [or transcript]. The GRanges object will have one or more ranges in it
  # corresponding to the de-duplicated ranges for all exonic regions in that
  # gene/transcript
  id_field <- ifelse(by_transcript, "ensembl_transcript_id", "ensembl_gene_id")
  genes <- GenomicRanges::split(ranges_obj,
                                f = mcols(ranges_obj)[[id_field]]) |>
    GenomicRanges::reduce()

  gene_length <- sum(width(genes))

  # Creates a DNAStringSetList where each list item is a set of sequences
  # (DNAStringSet) for a specific gene, corresponding to the exonic regions in
  # that gene
  exon_seqs <- getSeq(sequences, genes)

  gc <- sapply(exon_seqs, function(gene) {
    # Count all bases across all exonic sequences for this single gene
    gc_counts <- colSums(alphabetFrequency(gene))

    # Proportion of C and G bases. Some values in the sequence might be "N" and
    # we don't count "N" in the total so we can't use the calculated gene length
    # in the denominator
    sum(gc_counts[c("C", "G")]) / sum(gc_counts[c("A", "C", "G", "T")])
  })

  if (by_transcript) {
    # Collapse ranges_obj to transcript level and saves some extra fields
    fields_keep <- c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol",
                     "gene_biotype", "transcript_biotype", "chromosome_name")
  } else {
    # No transcript ID, collapses to gene level
    fields_keep <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype",
                     "chromosome_name")
  }

  fields_keep <- intersect(fields_keep, colnames(mcols(ranges_obj)))

  # Get data frame with one row per gene [or transcript] and add the calculated
  # length and GC content. Some columns are re-named to match how they are
  # named in Biomart
  gene_final <- ranges_obj |>
    as.data.frame() |>
    dplyr::select(all_of(fields_keep)) |>
    distinct() |>
    # Make sure gene length and gc are in the same order as the IDs in this data frame
    mutate(
      gene_length = gene_length[.data[[id_field]]],
      percent_gc_content = gc[.data[[id_field]]],
      # Some HGNC symbols might have been assigned the Ensembl ID because there's
      # no symbol for that gene. This fixes those to be blank.
      hgnc_symbol = case_when(
        grepl("ENSG00", hgnc_symbol) ~ "",
        .default = hgnc_symbol
      )
    ) |>
    # Sort for neatness and easier diff in case of updates
    arrange(ensembl_gene_id)

  if (by_transcript) {
    # Change "gene_length" to "transcript_length"
    gene_final <- gene_final |>
      dplyr::rename(transcript_length = gene_length)

    out_file <- file.path("data", "transcript_lengths_gc.csv")
  } else {
    out_file <- file.path("data", "gene_lengths_gc.csv")
  }

  # Save results for later
  write.csv(gene_final, out_file, row.names = FALSE, quote = FALSE)

  return(gene_final)
}


qc_json_test <- function() {
  js_txt <- readr::read_file("downloads/multiqc_data.json") # Mayo
  js_txt <- str_replace_all(js_txt, "NaN", str_escape('"NaN"'))
  js <- jsonlite::fromJSON(js_txt, simplifyVector = FALSE)

  data <- lapply(js$report_saved_raw_data, function(item) {
    do.call(rbind, item)
  })

  # multiqc_picard_dups: PERCENT_DUPLICATION
  # multiqc_rsem: Unalignable, Alignable, Filtered, Total, alignable_percent, Unique, Multi, Uncertain
  # multiqc_samtools_stats: reads_mapped, reads_mapped_and_paired, reads_unmapped, reads_duplicated, reads_QC_failed, bases_duplicated, mismatches, error_rate
  # multiqc_samtools_flagstat: ??
  # multiqc_samtools_idxstats: 1st number of chrX and chrY
  # multiqc_fastqc: %GC, total_deduplicated_percentage
  # multiqc_cutadapt: ??

  # Each item in idxstats is a list with 2 numbers. The first number is the
  # number of reads mapped to the chromosome, the second number is the length
  # of the chromosome. We want the first number for X and Y
  xy_stats <- apply(data$multiqc_samtools_idxstats, 1, function(row) {
    c(chrX = row[["chrX"]][[1]], chrY = row[["chrY"]][[1]])
  }) |>
    t() |>
    as.data.frame() |>
    mutate(total_reads = chrX + chrY,
           pct_X = chrX / total_reads,
           pct_Y = chrY / total_reads)

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
