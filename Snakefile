import pandas as pd
import numpy as np
import os
import re
import time
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("6.1.0")


##### load config and sample sheets #####

configfile: "bin/config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

species = pd.read_table("bin/species.tsv", dtype=str).set_index(["id"], drop=False)
validate(species, "schemas/species.schema.yaml")

# _plus_ is reserved to indicate hybrid genomes or genomes plus spikein
assert (not any(species.id.str.contains('_plus_'))), "_plus_ is reserved to indicate hybrid genomes or genomes plus spikein."

spikeins = pd.read_table("bin/spikeins.tsv", dtype=str)
validate(spikeins, "schemas/spikeins.schema.yaml")

hybrid_genomes = pd.read_table("bin/hybrid_genomes.tsv", dtype=str)
validate(hybrid_genomes, "schemas/hybrid_genomes.schema.yaml")

timestr = time.strftime("%Y-%m-%d_%H.%M.%S")
timestamp_dir = config["timestamp_dir"]
latest_symlink = timestamp_dir + 'latest'

# If there is a latest run (previous runs already present in the target directory) then add one to version number.
ref_version = int(re.findall(r'(?<=v)\d+$', os.readlink(latest_symlink))[0]) + 1 if os.path.exists(latest_symlink) else '1'
rule all:
    input:
       "{timestamp_dir}{timestr}_v{version}/rsync.done".format(timestamp_dir=timestamp_dir, timestr=timestr, version = ref_version)

rule timestamp_backup:
    input:
        #expand("data/{species.id}/sequence/{species.id}.fa", species=species.itertuples()),
        expand("data/{species_id}/annotation/{species_id}.basic.gtf", species_id=species[species.gene_basic_gtf.notnull()]['id'].values),
        expand("data/{species.id}/annotation/{species.id}.transcripts.fasta", species=species.itertuples()),
        expand("data/{species.id}/sequence/{species.id}.fa.fai", species=species.itertuples()),
        expand("data/{species.id}/sequence/{species.id}.dict", species=species.itertuples()),
        expand("data/{species.id}/indexes/star/SA", species=species.itertuples()),
        expand("data/{species.id}/indexes/bwa/{species.id}.bwt", species=species.itertuples()),
        expand("data/{species.id}/indexes/bowtie2/{species.id}.1.bt2", species=species.itertuples()),
        #expand("data/{species.id}/indexes/kb_lamanno/{species.id}.idx", species=species[-species["species"].str.contains("coli")].itertuples()),
        expand("data/{species.id}/indexes/bismark/{species.id}.fa", species=species.itertuples()),
        expand("data/{species.id}/indexes/kallisto/{species.id}", species=species.itertuples()),
        expand("data/{species.id}/indexes/salmon/{species.id}", species=species.itertuples()),
        #expand("data/{species.id}/annotation/{species.id}.gentrome.fa", species=species.itertuples()),
        expand("data/{species_id}/gatk_resource_bundle/done.txt", species_id=species[species.gatk_resource_bundle.notnull()]['id'].values),
        expand("data/{species.id}/blacklist/{species.blacklist_id}.bed", species=species[(species.replace(np.nan, '', regex=True)["blacklist"] != "") & (species.replace(np.nan, '', regex=True)["blacklist_id"] != "")].itertuples()),
        expand("data/{hybrid_genome_id}/indexes/bwa/{hybrid_genome_id}.bwt", hybrid_genome_id=hybrid_genomes.id),
        expand("data/{hybrid_genome_id}/sequence/{hybrid_genome_id}.{ext}", hybrid_genome_id=hybrid_genomes.id, ext=["fa.fai","dict"]),
        expand("data/{hybrid_genome_id}/indexes/star/SA", hybrid_genome_id=hybrid_genomes.id),
        expand("data/{hybrid_genome_id}/indexes/bowtie2/{hybrid_genome_id}.1.bt2", hybrid_genome_id=hybrid_genomes.id),
    output:
        flag=touch("{timestamp_dir}{{timestr}}/rsync.done".format(timestamp_dir=timestamp_dir)),
        outdir=directory("{timestamp_dir}{{timestr}}".format(timestamp_dir=timestamp_dir))
    log:
        stdout="logs/timestamp_backup/{timestr}.o",
        stderr="logs/timestamp_backup/{timestr}.e",

    benchmark:
        "benchmarks/timestamp_backup/{timestr}.txt"
    params:
        latest_link=lambda wildcards: "{timestamp_dir}/latest".format(timestamp_dir=timestamp_dir),
        sourceDir=config["sourceDir"],
        #backupPath=lambda wildcards: "{timestamp_dir}{timestr}".format(timestamp_dir=timestamp_dir, timestr=wildcards.timestr)
    threads:1
    resources:
        mem_gb=64
    envmodules:
        config["python3"]
    shell:
        """
        mkdir -p "{output.outdir}"

        python3 ./bin/python_scripts/main.py  -b "{params.latest_link}"  -s "{params.sourceDir}"  -d "{output.outdir}" 2>logs/python_script.log

#        rsync -rlDv \
#          -H \
#          --checksum \
#          --link-dest "{params.latest_link}" \
#          "{params.sourceDir}/" \
#          "{output.outdir}"
        
        rm -f "{params.latest_link}"
        ln -s "{output.outdir}" "{params.latest_link}" 
        """


# need this because both rules produce the same output (in Snakemake terminology, ambiguous rules)
ruleorder: cat_hybrid_seq > download_genome_fasta
ruleorder: cat_hybrid_gtf > download_genes_gtf

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
        dict="data/{species_id}/sequence/{species_id}.dict",
        temp=temp(directory("data/{species_id}/temp/")),

    log:
        stdout="logs/fai_and_dict/{species_id}.o",
        stderr="logs/fai_and_dict/{species_id}.e",

    benchmark:
        "benchmarks/fai_and_dict/{species_id}.txt"
    params:
    threads:1
    resources:
        mem_gb=30
    envmodules:
        config["samtools"],
        config["picard"]
    shell:
        """
        samtools faidx {input.genome_fa} 
            
        java \
        -Xms8g \
        -Xmx{resources.mem_gb}g \
        -Djava.io.tmpdir={output.temp} \
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
        config["STAR"]
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
        #genome_fa=lambda wildcards: "data/{species_id}/sequence/{species_id}.fa" if ({wildcards.species_id} in species['id'].values) else "data/withSpikein_{species_id}/sequence/{species_id}.fa",
        genome_fa="data/{species_id}/sequence/{species_id}.fa"
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
        config["bwa"]
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
        config["bowtie2"],
        config["python3"]
    shell:
        """
        bowtie2-build --threads {threads} {input.genome_fa} {params.outpref} 
            
        """


bismark_threads=16 # must be even number
rule bismark_idx:
    input:
        genome_fa="data/{species_id}/sequence/{species_id}.fa",
    output:
        "data/{species_id}/indexes/bismark/{species_id}.fa",
        directory("data/{species_id}/indexes/bismark/Bisulfite_Genome")
    log:
        stdout="logs/bismark_idx/{species_id}.o",
        stderr="logs/bismark_idx/{species_id}.e",
    benchmark:
        "benchmarks/bismark_idx/{species_id}.txt"
    params:
        outdir="data/{species_id}/indexes/bismark/",
        parallel=int(bismark_threads/2)
    threads:bismark_threads
    resources:
        mem_gb=100
    envmodules:
        config["bismark"]
    shell:
        """
        ln -rs {input.genome_fa} {output[0]}

        bismark_genome_preparation --bowtie2 --parallel {params.parallel} {params.outdir}
            
        """

rule kallisto_idx:
    input:
        tx_fa="data/{species_id}/annotation/{species_id}.transcripts.fasta"
    output:
        "data/{species_id}/indexes/kallisto/{species_id}"
    log:
        stdout="logs/kallisto_idx/{species_id}.o",
        stderr="logs/kallisto_idx/{species_id}.e",
    benchmark:
        "benchmarks/kallisto_idx/{species_id}.txt",
    params:
        out_idx="data/{species_id}/indexes/kallisto/{species_id}",
    threads: 4,
    resources:
        mem_gb=50,
    envmodules:
        config["kallisto"],
    shell:
        """
        kallisto index -i {params.out_idx} {input.tx_fa}
        """

rule build_salmon_gentrome_fa:
    input:
        genome_fa = "data/{species_id}/sequence/{species_id}.fa",
        tx_fa="data/{species_id}/annotation/{species_id}.transcripts.fasta",
    output:
        gentrome_fa ="data/{species_id}/annotation/{species_id}.gentrome.fa",
    log:
        stdout="logs/build_salmon_tx_fa/{species_id}.o",
        stderr="logs/build_salmon_tx_fa/{species_id}.e",
    benchmark:
        "benchmarks/build_salmon_tx_fa/{species_id}.txt",
    params:
        decoys_txt = "data/{species_id}/annotation/{species_id}.decoys.txt",
    threads: 1,
    resources:
        mem_gb=10,
    shell:
        """
        grep "^>" <{input.genome_fa} | cut -d " " -f 1 > {params.decoys_txt}
        sed -i.bak -e "s/>//g" {params.decoys_txt}
        cat {input.tx_fa} {input.genome_fa} > {output.gentrome_fa}
        """

# ruleorder: salmon_idx_from_gencode > salmon_idx_from_not_gencode

## if a transcript.fa from gencode references
# def get_species_from_gencode(wildcards):
#     input_list=list()
#     species_ids_list = species[species.tx_fasta.notnull()]['id'].values
#
#     for i in species_ids_list:
#         if "gencode" in i:
#             result = "data/{species_id}/annotation/{species_id}.gentrome.fa".format(species_id=i)
#             input_list.append(result)
#     return input_list

# rule salmon_idx_from_gencode:
rule salmon_idx:
    input:
        gentrome_fa="data/{species_id}/annotation/{species_id}.gentrome.fa"
        # get_species_from_gencode,
    output:
        directory("data/{species_id}/indexes/salmon/{species_id}")
    log:
        stdout="logs/salmon_idx/{species_id}.o",
        stderr="logs/salmon_idx/{species_id}.e",
    benchmark:
        "benchmarks/salmon_idx/{species_id}.txt",
    params:
        out_idx="data/{species_id}/indexes/salmon/{species_id}",
        decoys_txt= "data/{species_id}/annotation/{species_id}.decoys.txt",
    threads: 30,
    resources:
        mem_gb=50,
    envmodules:
        config["salmon"],
    shell:
        """        
        if [[ "{input.gentrome_fa}" == *"gencode"* ]]
        then 
            salmon index -t {input} -i {params.out_idx} --gencode -d {params.decoys_txt} -p {threads}
        else 
            salmon index -t {input} -i {params.out_idx} -d {params.decoys_txt} -p {threads} 
        fi
        """

# def get_species_from_not_gencode(wildcards):
#     input_list=list()
#     species_ids_list = species[species.tx_fasta.notnull()]['id'].values
#
#     for i in species_ids_list:
#         if "gencode" not in i:
#             result = "data/{species_id}/annotation/{species_id}.gentrome.fa".format(species_id=i)
#             input_list.append(result)
#     return input_list

# rule salmon_idx_from_not_gencode:
#     input:
#         get_species_from_not_gencode,
#     output:
#         directory("data/{species_id}/indexes/salmon/{species_id}")
#     log:
#         stdout="logs/salmon_idx/{species_id}.o",
#         stderr="logs/salmon_idx/{species_id}.e",
#     benchmark:
#         "benchmarks/salmon_idx/{species_id}.txt",
#     params:
#         out_idx="data/{species_id}/indexes/salmon/{species_id}",
#         decoys_txt= "data/{species_id}/annotation/{species_id}.decoys.txt",
#     threads: 30,
#     resources:
#         mem_gb=50,
#     envmodules:
#         config["salmon"],
#     shell:
#         """
#         salmon index -t {input} -i {params.out_idx} -d {params.decoys_txt} -p {threads}
#         """

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

        # Store the chromosome names in the fai file in a variable
        export fai_chroms=$(cut -f1 {input.genome_fai} | perl -lne 'push(@out_line, $_); END{{print join("\t", @out_line); print "\n"}}')
       
        # Go through the blacklist line by line. If a blacklist chrom name is not among the chrom names from the fai file, either add 'chr' prefix or take it away.
        perl -i -F"\\t" -lanpe  'BEGIN{{$fai_chroms=$ENV{{"fai_chroms"}}}} if($fai_chroms !~ /\\b$F[0]\\b/){{ if (/^chr/) {{s/^chr//;}} else {{$_="chr".$_; }}}}' {output}
        
        # Check that every chrom name in the blacklist is in the fai file
        perl -F"\\t" -ane  'BEGIN{{$fai_chroms=$ENV{{"fai_chroms"}}}} die "$F[0] not in .fai file" unless ($fai_chroms =~ /\\b$F[0]\\b/)' {output}
           
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
        config["kb-python"]
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
        config["parallel"],
        config["gsutil"]
    shell:
        """
        # download all files from the directory. Ignore sub-directories.
        gsutil ls $(echo {params.url} | perl -npe 's/https:\/\/storage.googleapis.com\//gs:\/\//') | grep -Pv '\/$' | perl -npe 's/gs:\/\//https:\/\/storage.googleapis.com\//' | parallel -k --will-cite --jobs {threads} 'wget -P {params.outdir} {{}}' 
            
        """


rule add_species_prefs_for_hybrid_fa:
    """
    Add prefix to species fasta to allow concatenation without duplicate sequence names.
    """
    input:
        "data/{species_id}/sequence/{species_id}.fa",
    output:
        "data/{species_id}/add_species_prefs_for_hybrid/{species_id}.{species_pref}_prefixed.fa",
    log:
        stdout="logs/add_species_prefs_for_hybrid_fa/{species_id}.{species_pref}.o",
        stderr="logs/add_species_prefs_for_hybrid_fa/{species_id}.{species_pref}.e"
    benchmark:
        "benchmarks/add_species_prefs_for_hybrid_fa/{species_id}.{species_pref}.txt"
    params:
        pref="{species_pref}"
    threads:1
    resources:
        mem_gb=1
    envmodules:
    shell:
        """
        # add species prefix to fasta
        perl -npe 's/^>/>{params.pref}_/' {input} > {output}
        """

rule add_species_prefs_for_hybrid_gtf:
    """
    Add prefix to species gtf to allow concatenation without duplicate sequence names.
    """
    input:
        "data/{species_id}/annotation/{species_id}.gtf"
    output:
        "data/{species_id}/add_species_prefs_for_hybrid/{species_id}.{species_pref}_prefixed.gtf"
    log:
        stdout="logs/add_species_prefs_for_hybrid_gtf/{species_id}.{species_pref}.o",
        stderr="logs/add_species_prefs_for_hybrid_gtf/{species_id}.{species_pref}.e"
    benchmark:
        "benchmarks/add_species_prefs_for_hybrid_gtf/{species_id}.{species_pref}.txt"
    params:
        pref="{species_pref}"
    threads:1
    resources:
        mem_gb=1
    envmodules:
    shell:
        """
        # add species prefix to gtf
        perl -npe 'next if /^#/; $_ = "{params.pref}_".$_' {input} > {output}
        """

def get_species_and_spikein_seqs(wildcards):
    # get the comma-separated species ids and prefixes
    species_ids = hybrid_genomes[hybrid_genomes.id == wildcards.hybrid_id]['species_ids'].values[0]
    species_prefs = hybrid_genomes[hybrid_genomes.id == wildcards.hybrid_id]['species_prefs'].values[0]

    # split into a list
    species_ids_list = species_ids.split(",")

    if len(species_ids_list) > 1:
        species_prefs_list = species_prefs.split(",")

        # get the original files for the first species
        first_species_paths = ["data/{id}/sequence/{id}.fa".format(id=species_ids_list[0])]

        # get the species-prefixed files for the 2nd and up to the last species
        other_species_paths = expand("data/{id}/add_species_prefs_for_hybrid/{id}.{pref}_prefixed.fa", zip, id=species_ids_list[1:], pref=species_prefs_list[1:])

        species_paths = first_species_paths + other_species_paths
    else:
        # else get the original species files
        species_paths = ["data/{id}/sequence/{id}.fa".format(id=species_ids)]

    # get spikein ids
    spikein_ids = hybrid_genomes[hybrid_genomes.id == wildcards.hybrid_id]['spikein_ids'].values[0]

    if not pd.isnull(spikein_ids):
        # split into a list
        spikein_ids_list = spikein_ids.split(",")
        # get the corresponding filenames
        spikein_filenames = [spikeins[spikeins.spikein_id == spikein_id]['spikein_fasta'].values[0] for spikein_id in spikein_ids_list]
        # setup paths
        spikein_paths = expand("bin/spikeins/sequence/{spikein_file}", spikein_file = spikein_filenames)

        input_paths = species_paths + spikein_paths
    else:
        input_paths = species_paths

    return input_paths

rule cat_hybrid_seq:
    input:
        get_species_and_spikein_seqs
    output:
        "data/{hybrid_id, .+_plus_.+}/sequence/{hybrid_id}.fa",
    log:
        stdout="logs/cat_hybrid_seq/{hybrid_id}.o",
        stderr="logs/cat_hybrid_seq/{hybrid_id}.e",

    benchmark:
        "benchmarks/cat_hybrid_seq/{hybrid_id}.txt"
    params:
    threads:1
    resources:
        mem_gb=1
    envmodules:
        config["seqtk"]
    shell:
        """
        cat {input} | seqtk seq -l 60 > {output}
        """


def get_species_and_spikein_gtfs(wildcards):
    # get the comma-separated species ids and prefixes
    species_ids = hybrid_genomes[hybrid_genomes.id == wildcards.hybrid_id]['species_ids'].values[0]
    species_prefs = hybrid_genomes[hybrid_genomes.id == wildcards.hybrid_id]['species_prefs'].values[0]

    # split into a list
    species_ids_list = species_ids.split(",")

    if len(species_ids_list) > 1:
        species_prefs_list = species_prefs.split(",")

        # get the original files for the first species
        first_species_paths = ["data/{id}/annotation/{id}.gtf".format(id=species_ids_list[0])]

        # get the species-prefixed files for the 2nd and up to the last species
        other_species_paths = expand("data/{id}/add_species_prefs_for_hybrid/{id}.{pref}_prefixed.gtf", zip, id=species_ids_list[1:], pref=species_prefs_list[1:])

        species_paths = first_species_paths + other_species_paths

    else:
        # else get the original species files
        species_paths = ["data/{id}/annotation/{id}.gtf".format(id=species_ids)]

    # get spikein ids
    spikein_ids = hybrid_genomes[hybrid_genomes.id == wildcards.hybrid_id]['spikein_ids'].values[0]

    if not pd.isnull(spikein_ids):
        # split into a list
        spikein_ids_list = spikein_ids.split(",")
        # get the corresponding filenames
        spikein_filenames = [spikeins[spikeins.spikein_id == spikein_id]['spikein_gtf'].values[0] for spikein_id in spikein_ids_list]
        # setup paths
        spikein_paths = ["bin/spikeins/annotation/{filename}".format(filename=x) for x in spikein_filenames if str(x) != 'nan']

        input_paths = species_paths + spikein_paths
    else:
        input_paths = species_paths

    return input_paths

# need to keep cat_hybrid_seq and cat_hybrid_gtf separate because some spikeins don't have a GTF file.
rule cat_hybrid_gtf:
    input:
        get_species_and_spikein_gtfs
    output:
        "data/{hybrid_id, .+_plus_.+}/annotation/{hybrid_id}.gtf",
    log:
        stdout="logs/cat_hybrid_gtf/{hybrid_id}.o",
        stderr="logs/cat_hybrid_gtf/{hybrid_id}.e",

    benchmark:
        "benchmarks/cat_hybrid_gtf/{hybrid_id}.txt"
    params:
    threads:1
    resources:
        mem_gb=1
    envmodules:
    shell:
        """
        cat {input} > {output}
        """


