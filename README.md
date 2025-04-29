# DiverseCohorts RNASeq Harmonization
Code to do quality control and harmonize RNA Seq data from Diverse Cohorts.

-----

# Data sources

Data sets from multiple groups were sequenced by three centers:

* **Mayo:** Mayo Clinic + Emory
* **New York Genome Center:** Columbia + MSSM
* **Rush:** Rush

Each center also sequenced several samples from each of the other centers
(sample swaps).

FASTQ or BAM files from each center were processed using the [nf-core/rnaseq](https://nf-co.re/rnaseq/)
pipeline. Data was aligned and quantified using the STAR/RSEM path, and quality
statistics were calculated with MultiQC.

Data was aligned to [**Gencode release 43 (primary assembly)**](https://www.gencodegenes.org/human/release_43.html).
GTF and FASTA files used for alignment are here:

* GTF: [gencode.v43.primary_assembly.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz)
* FASTA: [GRCh38.primary_assembly.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz)

----

# Warnings

## General

Batch names are not unique across data sets and need to be re-named prior to
analysis so they are unique. Info on batches:

* All data sequenced at Mayo (includes samples from Emory, Mayo, and sample
swaps) was run together in the same batches and not separated by data source.
Batches should be renamed to `Mayo_<batch>`.
* All data sequenced at Rush (includes samples from Rush and sample swaps) was
run together in the same batches and not separated by data source. Batches
should be renamed to `Rush_<batch>`.
* Data sequenced at NYGC (includes samples from Columbia, MSSM, and sample
swaps) was separated by data source and not run in the same batches. Batches
should be renamed to `NYGC_<source>_<batch>`.

RIN values are generally low for Rush samples, but DV200 is on par with the
other sites, so we use DV200 as the quality metric for all data in QC.

## Rush data

Specimen IDs "Div_487", "Div_498", "Div_500", "Div_503", "Div_508", "Div_509",
and "Div_545" are duplicated with two samples each. One of the duplicates for
each ID needs to be removed from the data prior to analysis. Removal
instructions:

* In the metadata, samples with these specimen IDs in batch "B74" should be
removed.
* In the counts matrix, two columns will have the same name (specimen ID) for
each of the IDs listed above. The **first column from the left** is always the
correct column to keep.
    * Note: In R, when reading with `read.table`, the columns are renamed as
    `<ID>` and `<ID>.1`. The column that is just `<ID>` **without** the ".1" is
    the correct column to keep.
* In Nextflow / QC files, the correct samples to keep are "Div_487_S128",
"Div_498_S126", "Div_500_S127", "Div_503_S118", "Div_508_S119",
"Div_509_S120", and "Div_545_S121".
    * Any samples with these specimen IDs and "S3xx" (where xx are any two
    numbers) should be discarded.

When we release this data we need to:

a) rename the unwanted duplicate specimen IDs to something else, or delete them,
and
b) make a note on the wiki and wherever else people might see as they access
this data that they should throw away the metadata for batch B74 for those
samples only.
