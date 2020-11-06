#PBS -l walltime=100:00:00
#PBS -l mem=8gb
#PBS -N prep_bbc_shared_workflow
#PBS -o logs/prep_bbc_shared_workflow.o
#PBS -e logs/prep_bbc_shared_workflow.e
#PBS -W umask=0022

cd ${PBS_O_WORKDIR}

module load bbc/snakemake/snakemake-5.17.0

# make temp directory for tools that need it
if [ ! -d "./temp/" ]
then
    mkdir ./temp/
fi


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

# run snakemake
snakemake \
--printshellcmds \
--use-envmodules \
--jobs 100 \
--cluster "qsub \
-V \
-q bbc \
-l nodes=1:ppn={threads} \
-l mem={resources.mem_gb}gb \
-l walltime=100:00:00 \
-o {log.stdout} \
-e {log.stderr}"

