# prep_bbc_shared

This is a Snakemake pipeline for preparing index files and other shared resources for use in the BBC at VAI. The goal is to automate and standardize the creation of index files for a number of different tools and across different species.

## How to run
`qsub -q bbc src/run_snakemake.sh`

## Directory structure
The file structure is:

* Date of creation
    * hg38_gencode
        * sequence
            * sequence.fa
            * sequence.fa.fai
            * sequence.dict
        * annotation
            * genes.gtf
        * indexes
            * star/
            * bwa/
            * bowtie2/
    * mm10_gencode
        * sequence
            * sequence.fa
            * sequence.fa.fai
            * sequence.dict
        * annotation
            * genes.gtf
        * indexes
            * star/
            * bwa/
            * bowtie2/

## The 'species' file

The pipeline is guided by a 'species' file that specifies the URLs for the genome fastas and the gene annotations, and could be augmented to include other things for future rules. The current columns are:

| Column #     | Column name      | Description                                                                                                   |
|----------    |--------------    |-----------------------------------------------------------------------------------------------------------    |
| 1.           | species          | Common name for the species e.g. mouse, human. This does not have to be unique.                               |
| 2.           | id               | Unique ID for the species e.g. mm10_gencode or hg38_gencode. This is used for directory and file naming.      |
| 3.           | genome_fasta     | URL for genome fasta.                                                                                         |
| 4.           | gene_gtf         | URL for GTF.                                                                                                  |
| 5.           | blacklist        | URL for ENCODE blacklist.                                                                                     |
| 6.           | blacklist_id     | ID for renaming blacklist after downloading. Identifies the version of the blacklist e.g. hg38_encode_v3.     |

The layout of this file and the pipeline allows for stright-forward the specification of different variants of references for the same species. For example, there could be a row for 'human hg38_ensembl' and another row for 'human hg19_ensembl' in the species file. This would result in the creation of two sub-directories: 'hg38_ensembl' and 'hg19_ensembl', with their own genome fastas, annotation files and index files.

One may wish to make a nonstandard variation of the hg38 reference using a customized genome fasta. In this case, one could supply the genome fasta file manually in the appropriately named subdirectory and the Snakemake pipeline would skip the fasta downloading rule. Similarly, one could use a custom GTF and move it manually into the appropriate subdirectory and the pipeline will skip downloading the GTF. One could also place a symlink in place of the fasta or GTF files.

## The workflow

![Workflow](./logs/rulegraph.png)

