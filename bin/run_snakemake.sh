#PBS -l walltime=100:00:00
#PBS -l mem=8gb
#PBS -N prep_bbc_shared_workflow
#PBS -o logs/prep_bbc_shared_workflow.o
#PBS -e logs/prep_bbc_shared_workflow.e
#PBS -W umask=002

cd ${PBS_O_WORKDIR}

snakemake_module="bbc/snakemake/snakemake-6.1.0"

module load $snakemake_module

# make temp directory for tools that need it
#if [ ! -d "./temp/" ]
#then
#    mkdir ./temp/
#fi


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
-p \
--latency-wait 20 \
--use-envmodules \
--jobs 50 \
--cluster "ssh ${PBS_O_LOGNAME}@submit 'module load ${snakemake_module}; cd ${PBS_O_WORKDIR}; qsub \
-V \
-q bbc \
-l nodes=1:ppn={threads} \
-l mem={resources.mem_gb}gb \
-l walltime=100:00:00 \
-W umask=002 \
-o {log.stdout} \
-e {log.stderr}'" \
--conda-frontend conda
