import pandas as pd
import numpy as np
import os
import re
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.11.0")


##### load config and sample sheets #####

#configfile: "src/config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

species = pd.read_table("src/species.tsv", dtype=str).set_index(["id"], drop=False)
validate(species, "schemas/species.schema.yaml")

rule all:
    input:
        #expand("data/{species.id}/sequence/{species.id}.fa", species=species.itertuples()),
        #expand("data/{species.id}/annotation/{species.id}.gtf", species=species.itertuples()),
        expand("data/{species.id}/sequence/{species.id}.fa.fai", species=species.itertuples()),
        expand("data/{species.id}/sequence/{species.id}.fa.dict", species=species.itertuples()),
        expand("data/{species.id}/indexes/star/SA", species=species.itertuples()),
        expand("data/{species.id}/indexes/bwa/{species.id}.bwt", species=species.itertuples()),
        expand("data/{species.id}/indexes/bowtie2/{species.id}.1.bt2", species=species.itertuples()),

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
        mem_mb=16000
    envmodules:
    shell:
        """
        # download the file
        wget {params.url} -O {output} \
        2>{log.stderr} 1>{log.stdout} 
            
        # if gzipped, decompress it. Need to give it a .gz suffix or gnzip will fail 
        if (file {output} | grep -q 'gzip compressed' ) ; then
            mv {output} {output}.gz
            gunzip {output}.gz 2>>{log.stderr} 1>>{log.stdout}
        fi
 
        """

rule fai_and_dict:
    input: 
        genome_fa="data/{species_id}/sequence/{species_id}.fa", 
    output:
        fai="data/{species_id}/sequence/{species_id}.fa.fai",
        dict="data/{species_id}/sequence/{species_id}.fa.dict"
        
    log:
        stdout="logs/fai_and_dict/{species_id}.o",
        stderr="logs/fai_and_dict/{species_id}.e",

    benchmark:
        "benchmarks/fai_and_dict/{species_id}.txt"
    params:
        temp="temp/"
    threads:1
    resources:
        mem_mb=30000
    envmodules:
        "bbc/samtools/samtools-1.9",
        "bbc/picard/picard-2.21.4-SNAPSHOT"
    shell:
        """
        samtools faidx {input.genome_fa} \
        2>{log.stderr} 1>{log.stdout} 
            
        java \
        -Xms8g \
        -Xmx{resources.mem_mb}m \
        -Djava.io.tmpdir={params.temp} \
        -jar $PICARD \
        CreateSequenceDictionary \
        REFERENCE={input.genome_fa} \
        OUTPUT={output.dict} \
        2>>{log.stderr} 1>>{log.stdout}
 
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
        mem_mb=16000
    envmodules:
    shell:
        """
        # download the file
        wget {params.url} -O {output} \
        2>{log.stderr} 1>{log.stdout} 
            
        # if gzipped, decompress it. Need to give it a .gz suffix or gnzip will fail 
        if (file {output} | grep -q 'gzip compressed' ) ; then
            mv {output} {output}.gz
            gunzip {output}.gz 2>>{log.stderr} 1>>{log.stdout}
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
        mem_mb=100000
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
            --sjdbOverhang {params.sjdb_overhang} \
            2>{log.stderr} 1>{log.stdout} 
            
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
        mem_mb=100000
    envmodules:
        "bbc/bwa/bwa-0.7.17"
    shell:
        """
        #ln -s "$(pwd)/{input.genome_fa}" {params.outpref}
        
        bwa index \
        -p {params.outpref} \
        {input.genome_fa} \
        2>{log.stderr} 1>{log.stdout} 
            
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
        mem_mb=100000
    envmodules:
        "bbc/bowtie2/bowtie2-2.4.1",
        "bbc/python3/python-3.8.1"
    shell:
        """
        bowtie2-build --threads {threads} {input.genome_fa} {params.outpref} \
        2>{log.stderr} 1>{log.stdout}
            
        """


