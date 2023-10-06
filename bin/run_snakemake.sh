#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J prep_bbc_shared_workflow
#SBATCH -o logs/prep_bbc_shared_workflow.o
#SBATCH -e logs/prep_bbc_shared_workflow.e
#SBATCH --ntasks 1
#SBATCH --time 100:00:00
#SBATCH --mem=8G
#SBATCH --partition=long

cd $SLURM_SUBMIT_DIR

snakemake_module="bbc2/snakemake/snakemake-7.25.0"

module load $snakemake_module

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
if [ ! -d "./logs/runs" ]
then
    mkdir -p ./logs/runs
fi
snakemake --use-conda -n > logs/runs/prep_bbc_shared_workflow_${TIME}.txt
snakemake --dag | dot -Tpng > logs/dag.png
snakemake --filegraph | dot -Tpng > logs/filegraph.png
snakemake --rulegraph | dot -Tpng > logs/rulegraph.png


echo "Start snakemake workflow." >&1                   
echo "Start snakemake workflow." >&2     

snakemake \
-p \
--latency-wait 20 \
--snakefile 'Snakefile' \
--use-envmodules \
--jobs 50 \
--cluster "mkdir -p logs/{rule}; sbatch \
-p ${SLURM_JOB_PARTITION} \
--export=ALL \
--nodes 1 \
--ntasks-per-node {threads} \
--mem={resources.mem_gb}G \
-t 100:00:00 \
-o logs/{rule}/{resources.log_prefix}.o \
-e logs/{rule}/{resources.log_prefix}.e" # SLURM hangs if output dir does not exist, so we create it before running sbatch on the snakemake jobs.
#--slurm \
#--default-resources slurm_account=${SLURM_JOB_USER} slurm_partition=${SLURM_JOB_PARTITION}

echo "snakemake workflow done." >&1                   
echo "snakemake workflow done." >&2                

