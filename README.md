# prep_bbc_shared

This is a Snakemake pipeline for preparing index files and other shared resources for use in the BBC at VAI. The goal is to automate and standardize the creation of index files for a number of different tools and across different species.

Table of Contents
=================

   * [How to run](#how-to-run)
   * [Incremental backups after each run](#incremental-backups-after-each-run)
   * [The 'species' file](#the-species-file)
   * [The workflow](#the-workflow)
   * [Hybrid genomes and spike-in references](#hybrid-genomes-and-spike-in-references)

# How to run

1. For BBC members, don't clone this repo, just work from the previously cloned directory.

   ```cd /varidata/research/projects/bbc/research/prep_bbc_shared_current```

2. Add new rows to the 'species.tsv', 'hybrid_genomes.tsv', and/or 'spikeins.tsv' as needed. If you simply want to compute the index files for a new species, you only need to add a row with the relevant information to 'species.tsv'. For more information about the other two files, see [Hybrid genomes and spike-in references](#hybrid-genomes-and-spike-in-references).

3. Run `snakemake -npr` and look over all the jobs that will be run. Typically, only jobs related to the new information entered in step 1 need to be run. To rerun specific species or hybrid genomes (for example if making the index with a new aligner version), you will need to delete the corresponding directories in `data/`.

4. To run the above jobs, type `qsub -q bbc bin/run_snakemake.sh`.

# Incremental backups after each run

In order to allow specific index or resource files to be permanently accessible, after every run, this workflow creates a timestamped copy of the entire workflow directory at `/varidata/research/projects/bbc/versioned_references/<timestamp>_v<version_num>`. This backup mechanism is specified by the 'timestamp_backup' rule.

To save space, an incremental backup strategy is used. Inside the `versioned_references/` directory is a symlink called 'latest' that always points to the latest timestamped directory. During each backup, the checksums of the files in 'latest' are compared with those of the newly generated files. If the files are unchanged, hard links to the 'latest' directory are created in the new backup. Otherwise, files will be copied for the new backup. At the end of the 'timestamp_backup' rule, the 'latest' symlink is updated to the newest time-stamped directory.

# The 'species' file

The pipeline is guided by a 'species' file that specifies the URLs for the genome fastas and the gene annotations, and could be augmented to include other things for future rules. The current columns are:

| Column #     | Column name              | Description                                                                                                   |
|----------    |--------------            |-----------------------------------------------------------------------------------------------------------    |
| 1.           | species                  | Common name for the species e.g. mouse, human. This does not have to be unique.                               |
| 2.           | id                       | Unique ID for the species e.g. mm10_gencode or hg38_gencode. This is used for directory and file naming.      |
| 3.           | genome_fasta             | URL for genome fasta.                                                                                         |
| 4.           | gene_gtf                 | URL for GTF.                                                                                                  |
| 5.           | blacklist                | URL for ENCODE blacklist.                                                                                     |
| 6.           | blacklist_id             | ID for renaming blacklist after downloading. Identifies the version of the blacklist e.g. hg38_encode_v3.     |
| 7.           | gatk_resource_bundle     | URL for GATK resource bundle. _Optional_.                                                                     |
| 8.           | gene_basic_gtf           | URL for GTF for 'basic' annotations as defined by GENCODE. Not used for indexing.  _Optional_.                |


The layout of this file and the pipeline allows for stright-forward the specification of different variants of references for the same species. For example, there could be a row for 'human hg38_ensembl' and another row for 'human hg19_ensembl' in the species file. This would result in the creation of two sub-directories: 'hg38_ensembl' and 'hg19_ensembl', with their own genome fastas, annotation files and index files.

One may wish to make a nonstandard variation of the hg38 reference using a customized genome fasta. In this case, one could supply the genome fasta file manually in the appropriately named subdirectory and the Snakemake pipeline would skip the fasta downloading rule. Similarly, one could use a custom GTF and move it manually into the appropriate subdirectory and the pipeline will skip downloading the GTF. One could also place a symlink in place of the fasta or GTF files.

# The workflow

![Workflow](./logs/rulegraph.png)

# Hybrid genomes and spike-in references

Species specified in `species.tsv` can be `cat` together for hybrid reference genomes. Similarly, `bin/spikeins/` contains sequences and annotations for spikein sequences that can be combined with any species in the `species.tsv` file. This feature depends on the `bin/hybrid_genomes.tsv` file. The columns for `hybrid_genomes.tsv`are as follows:

| Column #     | Column name              | Description                                                                                                   |
|----------    |--------------            |-----------------------------------------------------------------------------------------------------------    |
| 1.           | id                       | A unique ID for the hybrid genome. Use the format 'species_plus_spikin'. e.g. mm10_gencode_plus_ERCC92        |
| 2.           | species_ids              | Comma-separated list of species ids to include. Must match 'id' column in 'species.tsv'.                      |
| 3.           | species_prefs            | Comma-separated list of prefices to prepend to species chromosomes. Length must match the # of species in 'species_id' but actually the original chromosome names are retained for the first species.                                                           |
| 4.           | spikein_ids              | Comma-separated list of spike-in IDs. Must match 'spikein_id' column in 'spikeins.tsv'.                       |

Because the chromosome names are not changed for the first species in a hybrid reference, it should retain compatibility with blacklists and other resources for the original reference.
