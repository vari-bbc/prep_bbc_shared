# prep_bbc_shared

This is a Snakemake pipeline for preparing index files and other shared resources for use in the BBC at VAI. The goal is to automate and standardize the creation of index files for a number of different tools and across different species.

The file structure is:

* Date of creation
 * Human
  * sequence
  * annotation
  * indexes
   * star
   * bwa
 * Mouse
  * sequence
  * annotation
  * indexes
   * star
   * bwa
