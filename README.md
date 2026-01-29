# Diverse Cohorts RNASeq Harmonization

Code to do quality control and harmonize RNA Seq data from Diverse Cohorts. The
current outputs of this pipeline are:
-   count matrices containing only samples that passed QC
-   CQN-normalized matrices that have been adjusted for gene length and GC content

--------------------------------------------------------------------------------

## Pipeline information

Each step in the pipeline has been numbered so the order of operations is clear:

-   **00_Install.R:** installs all packages required to run the pipeline

-   **01_Download_QC_Files.R:** Downloads FastQC and MultiQC files for each
    dataset [from Synapse](https://www.synapse.org/Synapse:syn68755487) and
    calculates gene length and GC content for each gene.

-   **02\_\<dataset\>\_QC.qmd:** Quarto notebooks, one per data set, that
    perform QC on the raw counts.

    -   Each notebook file contains a series of code chunks that are empty in
        the notebook file but reference sections of **qc_functions.R**. The code
        from these sections is rendered seamlessly into the notebook when it is
        rendered to HTML.

    -   Things are set up this way to ensure that all four data sets use the
        exact same QC code and ensures there are no copy/paste errors between
        notebooks.

-   **03\_\<dataset\>\_CQN.qmd:** Quarto notebooks, one per data set, that
    perform conditional quantile normalization (CQN) on the post-QC count
    matrices.

    -   Like step 2, these notebooks reference code from **cqn_functions.R**,
        which is rendered into the notebook.

-   **04_Find_Models.R:** *(experimental)* Tries to algorithmically determine
    the best variables to include in a model for residualizing out technical
    variation but leaving biological variation in the data. This script should
    be considered unfinished.

-   **05\_\<dataset\>\_Regression.qmd:** *(experimental)* Quarto notebooks that
    try to regress technical variation out of the data based on the models
    recommended by step 4 plus my own judgment. These notebooks should be
    considered unfinished due to unresolved exploration into untangling batch
    from race in the MSSM and Rush data sets (see [Notes and Warnings] below).

Other files:

-   **config.yml:** Contains general and dataset-specific configuration like
    Synapse IDs for necessary files, thresholds for QC, and any samples that
    need removing prior to QC.

-   **cqn_functions.R:** Contains the code to perform CQN, formatted so that
    each piece gets rendered into the CQN notebooks.

-   **helper_functions.R:** Functions used across various parts of the pipeline
    that don't fit into a single category. This code is sourced by the other R
    scripts and notebooks but is not rendered into the notebooks.

-   **qc_functions.R:** Contains the code to perform QC, formatted so that each
    piece gets rendered into the QC notebooks.

-   **scratch/:** This folder contains early experimental code that I do not
    want to delete yet. No code in this folder is sourced or used in the final
    pipeline.

--------------------------------------------------------------------------------

## Data sources

Data sets from multiple groups were sequenced by three centers:

-   **Mayo:** Mayo Clinic + Emory
-   **New York Genome Center:** Columbia + MSSM
-   **Rush:** Rush

Each center also sequenced several samples from each of the other centers
(sample swaps).

FASTQ or BAM files from each center were processed using the
[nf-core/rnaseq](https://nf-co.re/rnaseq/) pipeline. Data was aligned and
quantified using the STAR/RSEM path, and quality statistics were calculated with
MultiQC.

Data was aligned to [**Gencode release 43 (primary
assembly)**](https://www.gencodegenes.org/human/release_43.html). GTF and FASTA
files used for alignment are here:

-   GTF:
    [gencode.v43.primary_assembly.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz)
-   FASTA:
    [GRCh38.primary_assembly.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz)

--------------------------------------------------------------------------------

# Notes and Warnings

## General

RIN values are generally low for Rush samples due to the RNA prep method, but
DV200 is on par with the other sites, so we use DV200 as the quality metric for
all data in QC.

The batch design for MSSM and Rush samples is **highly entangled** with race.
There are several library or sequencing batches in each dataset that contain
samples from mostly or only one race, and also some batches that only have 1-2
samples total.

-   I have done some exploration into regressing out technical variables, but
    the batch design makes it difficult to produce a residualized matrix that
    has batch removed but does not remove potential variation due to race.

-   Any analysis done that uses regression should instead use the normalized
    matrix and include both race and batch in the model if race is a variable of
    interest.

## Rush data

Seven specimen IDs are duplicated with two samples each, due to the same
sample/specimen being sequenced twice. One of the duplicates for each ID needs
to be removed from the data prior to analysis. Removal instructions:

-   In the assay metadata, the samples to remove from each pair of duplicates
    are marked as TRUE in the "exclude" column of the assay metadata.
    
-   In the raw counts matrix, two columns will have the same name (specimen ID)
    for each of the IDs belonging to "exclude" samples above. The **first column
    from the left** in each pair is always the correct column to keep.
    
    -   Note: In R, when reading with `read.table`, the columns are renamed as
        `<ID>` and `<ID>.1`. The column that is just `<ID>` **without** the ".1"
        is the correct column to keep.
        
    -   In the post-QC matrices, **these samples have already been removed.**
    
-   In Nextflow / QC files, the specimenIDs have slightly different formatting.
    The config.yml file contains the correct samples to remove, in the format
    matching these files.
