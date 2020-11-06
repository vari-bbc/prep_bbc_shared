import pandas as pd
import numpy as np
import os
import re
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.11.0")


##### load config and sample sheets #####

#configfile: "bin/config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

species = pd.read_table("bin/species.tsv", dtype=str).set_index(["id"], drop=False)
validate(species, "schemas/species.schema.yaml")

rule all:
    input:
        #expand("data/{species.id}/sequence/{species.id}.fa", species=species.itertuples()),
        expand("data/{species_id}/annotation/{species_id}.basic.gtf", species_id=species[species.gene_basic_gtf.notnull()]['id'].values),
        expand("data/{species.id}/annotation/{species.id}.transcripts.fasta", species=species.itertuples()),
        expand("data/{species.id}/sequence/{species.id}.fa.fai", species=species.itertuples()),
        expand("data/{species.id}/sequence/{species.id}.dict", species=species.itertuples()),
        expand("data/{species.id}/indexes/star/SA", species=species.itertuples()),
        expand("data/{species.id}/indexes/bwa/{species.id}.bwt", species=species.itertuples()),
        expand("data/{species.id}/indexes/bowtie2/{species.id}.1.bt2", species=species.itertuples()),
        expand("data/{species.id}/indexes/kb_lamanno/{species.id}.idx", species=species[-species["species"].str.contains("coli")].itertuples()),
        expand("data/{species_id}/gatk_resource_bundle/done.txt", species_id=species[species.gatk_resource_bundle.notnull()]['id'].values),
        expand("data/{species.id}/blacklist/{species.blacklist_id}.bed", species=species[(species.replace(np.nan, '', regex=True)["blacklist"] != "") & (species.replace(np.nan, '', regex=True)["blacklist_id"] != "")].itertuples()),

rule download_genome_fasta:
    input: 

    output:
        "data/{species_id}/sequence/{species_id}.fa"
        
    log:
        stdout="logs/download_genome_fasta/{species_id}.o",
        stderr="logs/download_genome_fasta/{species_id}.e",

    benchmark:
        "benchmarks/download_genome_fasta/{species_id}.txt"
    params:
        url=lambda wildcards: species.loc[species.id == wildcards.species_id,'genome_fasta'].values[0]
    threads:1
    resources:
        mem_gb=16
    envmodules:
    shell:
        """
        # download the file
        wget {params.url} -O {output} 
            
        # if gzipped, decompress it. Need to give it a .gz suffix or gnzip will fail 
        if (file {output} | grep -q 'gzip compressed' ) ; then
            mv {output} {output}.gz
            gunzip {output}.gz
        fi
 
        """

rule fai_and_dict:
    input: 
        genome_fa="data/{species_id}/sequence/{species_id}.fa", 
    output:
        fai="data/{species_id}/sequence/{species_id}.fa.fai",
        dict="data/{species_id}/sequence/{species_id}.dict"
        
    log:
        stdout="logs/fai_and_dict/{species_id}.o",
        stderr="logs/fai_and_dict/{species_id}.e",

    benchmark:
        "benchmarks/fai_and_dict/{species_id}.txt"
    params:
        temp="temp/"
    threads:1
    resources:
        mem_gb=30
    envmodules:
        "bbc/samtools/samtools-1.9",
        "bbc/picard/picard-2.21.4-SNAPSHOT"
    shell:
        """
        samtools faidx {input.genome_fa} 
            
        java \
        -Xms8g \
        -Xmx{resources.mem_gb}g \
        -Djava.io.tmpdir={params.temp} \
        -jar $PICARD \
        CreateSequenceDictionary \
        REFERENCE={input.genome_fa} \
        OUTPUT={output.dict} 
 
        """

rule download_genes_gtf:
    input: 

    output:
        "data/{species_id}/annotation/{species_id}.gtf"
        
    log:
        stdout="logs/download_genes_gtf/{species_id}.o",
        stderr="logs/download_genes_gtf/{species_id}.e",
 
    benchmark:
        "benchmarks/download_genes_gtf/{species_id}.txt"
    params:
        url=lambda wildcards: species.loc[species.id == wildcards.species_id,'gene_gtf'].values[0]
    threads:1
    resources:
        mem_gb=16
    envmodules:
    shell:
        """
        # download the file
        wget {params.url} -O {output} 
            
        # if gzipped, decompress it. Need to give it a .gz suffix or gnzip will fail 
        if (file {output} | grep -q 'gzip compressed' ) ; then
            mv {output} {output}.gz
            gunzip {output}.gz
        fi
           
        """

rule download_basic_genes_gtf:
    """
    Download gtf with only the representative transcript for each gene and focuses on protein coding genes.
    """
    input: 

    output:
        "data/{species_id}/annotation/{species_id}.basic.gtf"
        
    log:
        stdout="logs/download_basic_genes_gtf/{species_id}.o",
        stderr="logs/download_basic_genes_gtf/{species_id}.e",
 
    benchmark:
        "benchmarks/download_basic_genes_gtf/{species_id}.txt"
    params:
        url=lambda wildcards: species.loc[species.id == wildcards.species_id,'gene_basic_gtf'].values[0]
    threads:1
    resources:
        mem_gb=16
    envmodules:
    shell:
        """
        # download the file
        wget {params.url} -O {output} 
            
        # if gzipped, decompress it. Need to give it a .gz suffix or gnzip will fail 
        if (file {output} | grep -q 'gzip compressed' ) ; then
            mv {output} {output}.gz
            gunzip {output}.gz
        fi
           
        """

rule download_tx_fasta:
    input: 

    output:
        "data/{species_id}/annotation/{species_id}.transcripts.fasta"
        
    log:
        stdout="logs/download_tx_fasta/{species_id}.o",
        stderr="logs/download_tx_fasta/{species_id}.e",

    benchmark:
        "benchmarks/download_tx_fasta/{species_id}.txt"
    params:
        url=lambda wildcards: species.loc[species.id == wildcards.species_id,'tx_fasta'].values[0]
    threads:1
    resources:
        mem_gb=16
    envmodules:
    shell:
        """
        # download the file
        wget {params.url} -O {output} 
            
        # if gzipped, decompress it. Need to give it a .gz suffix or gnzip will fail 
        if (file {output} | grep -q 'gzip compressed' ) ; then
            mv {output} {output}.gz
            gunzip {output}.gz 
        fi
           
        """


rule star_idx:
    input: 
        genome_fa="data/{species_id}/sequence/{species_id}.fa",
        genes_gtf="data/{species_id}/annotation/{species_id}.gtf"
    output:
        "data/{species_id}/indexes/star/Log.out",
        "data/{species_id}/indexes/star/SA",
        "data/{species_id}/indexes/star/SAindex"
        
    log:
        stdout="logs/star_idx/{species_id}.o",
        stderr="logs/star_idx/{species_id}.e",

    benchmark:
        "benchmarks/star_idx/{species_id}.txt"
    params:
        sjdb_overhang=100,
        outpref="data/{species_id}/indexes/star/"
    threads:16
    resources:
        mem_gb=100
    envmodules:
        "bbc/STAR/STAR-2.7.3a"
    shell:
        """
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --outFileNamePrefix {params.outpref} \
            --genomeDir {params.outpref} \
            --genomeFastaFiles {input.genome_fa} \
            --sjdbGTFfile {input.genes_gtf} \
            --sjdbOverhang {params.sjdb_overhang} 
            
        """

rule bwa_idx:
    input: 
        genome_fa="data/{species_id}/sequence/{species_id}.fa",
    output:
        "data/{species_id}/indexes/bwa/{species_id}.bwt",
        "data/{species_id}/indexes/bwa/{species_id}.sa" 
    log:
        stdout="logs/bwa_idx/{species_id}.o",
        stderr="logs/bwa_idx/{species_id}.e",

    benchmark:
        "benchmarks/bwa_idx/{species_id}.txt"
    params:
        outpref="data/{species_id}/indexes/bwa/{species_id}"
    threads:4
    resources:
        mem_gb=100
    envmodules:
        "bbc/bwa/bwa-0.7.17"
    shell:
        """
        #ln -s "$(pwd)/{input.genome_fa}" {params.outpref}
        
        bwa index \
        -p {params.outpref} \
        {input.genome_fa} 
            
        """

rule bowtie2_idx:
    input: 
        genome_fa="data/{species_id}/sequence/{species_id}.fa",
    output:
        "data/{species_id}/indexes/bowtie2/{species_id}.1.bt2",
        "data/{species_id}/indexes/bowtie2/{species_id}.2.bt2" 
    log:
        stdout="logs/bowtie2_idx/{species_id}.o",
        stderr="logs/bowtie2_idx/{species_id}.e",
    benchmark:
        "benchmarks/bowtie2_idx/{species_id}.txt"
    params:
        outpref="data/{species_id}/indexes/bowtie2/{species_id}"
    threads:8
    resources:
        mem_gb=100
    envmodules:
        "bbc/bowtie2/bowtie2-2.4.1",
        "bbc/python3/python-3.8.1"
    shell:
        """
        bowtie2-build --threads {threads} {input.genome_fa} {params.outpref} 
            
        """

# We need the genome fai as input because we check that the chromosome names are compatible between the blacklist and the genome
rule download_blacklist:
    input: 
        genome_fai="data/{species_id}/sequence/{species_id}.fa.fai"
    output:
        "data/{species_id}/blacklist/{blacklist_id}.bed"
    log:
        stdout="logs/blacklist/{species_id}.{blacklist_id}.o",
        stderr="logs/blacklist/{species_id}.{blacklist_id}.e",

    benchmark:
        "benchmarks/blacklist/{species_id}.{blacklist_id}.txt"
    params:
        url=lambda wildcards: species.loc[species.id == wildcards.species_id, 'blacklist'].values[0],
        temp_chroms="data/{species_id}/blacklist/{blacklist_id}_chroms.temp"
    threads:1
    resources:
        mem_gb=16
    envmodules:
    shell:
        """
        # download the file
        wget {params.url} -O {output} 
            
        # if gzipped, decompress it. Need to give it a .gz suffix or gunzip will fail 
        if (file {output} | grep -q 'gzip compressed' ) ; then
            mv {output} {output}.gz
            gunzip {output}.gz 
        fi
           
        # make sure chromosome names are compatible between genome fasta and blacklist
        paste <(cut -f1 {input.genome_fai} | grep -Pi "^(chr)?[0-9XY]{{1,2}}" | sort) <(cut -f1 {output} | grep -Pi "^(chr)?[0-9XY]{{1,2}}" | sort | uniq) > {params.temp_chroms} 
        perl -lane 'die "Chromosomes in genome fasta and blacklist do not match. See chroms.temp in blacklist directory." unless $F[0] eq $F[1]' {params.temp_chroms} 
        rm {params.temp_chroms}
        """

rule kb_lamanno:
    input:
        genome_fa="data/{species_id}/sequence/{species_id}.fa",
        genes_gtf="data/{species_id}/annotation/{species_id}.gtf"
    output:
        idx="data/{species_id}/indexes/kb_lamanno/{species_id}.idx",
        t2g="data/{species_id}/indexes/kb_lamanno/{species_id}.transcripts_to_genes.txt",
        cdna="data/{species_id}/indexes/kb_lamanno/{species_id}.cdna.fa",
        introns="data/{species_id}/indexes/kb_lamanno/{species_id}.intron.fa",
        cdna_tr2cap="data/{species_id}/indexes/kb_lamanno/{species_id}.cdna_tr2cap.txt",
        introns_tr2cap="data/{species_id}/indexes/kb_lamanno/{species_id}.intron_tr2cap.txt"
    log:
        stdout="logs/kb_lamanno/{species_id}.o",
        stderr="logs/kb_lamanno/{species_id}.e",
    benchmark:
        "benchmarks/kb_lamanno/{species_id}.txt"
    envmodules:
        "bbc/kb-python/kb-python-0.24.4"
    resources:
        mem_gb=160
    shadow: "shallow"
    threads:4
    params:
    shell:
        """
        # because kb ref 0.24.4 makes and uses a tmp/ in the working directory, there will be collisions if you run kb ref for multiple species at the same time. The files in tmp/ are generically named, sorted.fa and sorted.gtf. 
        #The devel version of kb allows user to specify dustom tmp directory location but not in the latest release. 
        # as a workaround we use the shadow rule feature of snakemake. "Shadow rules result in each execution of the rule to be run in isolated temporary directories."
        kb ref \
        -i {output.idx} \
        -g {output.t2g} \
        -f1 {output.cdna} \
        -f2 {output.introns} \
        -c1 {output.cdna_tr2cap} \
        -c2 {output.introns_tr2cap} \
        --lamanno \
        {input.genome_fa} \
        {input.genes_gtf} 
        """

rule download_gatk_resource_bundle:
    input: 
    output:
        touch=touch("data/{species_id}/gatk_resource_bundle/done.txt"),
    log:
        stdout="logs/gatk_resource_bundle/{species_id}.o",
        stderr="logs/gatk_resource_bundle/{species_id}.e",

    benchmark:
        "benchmarks/gatk_resource_bundle/{species_id}.txt"
    params:
        url=lambda wildcards: species[species.id == wildcards.species_id].gatk_resource_bundle.values[0],
        outdir="data/{species_id}/gatk_resource_bundle/"
    threads:4
    resources:
        mem_gb=64
    envmodules:
        "bbc/parallel/parallel-20191122",
        "bbc/gsutil/gsutil-4.52"
    shell:
        """
        # download all files from the directory. Ignore sub-directories.
        gsutil ls $(echo {params.url} | perl -npe 's/https:\/\/storage.googleapis.com\//gs:\/\//') | grep -Pv '\/$' | perl -npe 's/gs:\/\//https:\/\/storage.googleapis.com\//' | parallel -k --will-cite --jobs {threads} 'wget -P {params.outdir} {{}}' 
            
        """

