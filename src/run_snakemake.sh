#PBS -l walltime=100:00:00
#PBS -l mem=8gb
#PBS -N prep_bbc_shared_workflow
#PBS -o logs/prep_bbc_shared_workflow.o
#PBS -e logs/prep_bbc_shared_workflow.e

cd ${PBS_O_WORKDIR}

# load snakemake module which is actually a python virtualenv that also loads pandas (and its own python3)
#module load bbc/snakemake/snakemake-5.10.0
#module load bbc/snakemake/snakemake-5.8.2
. /secondary/projects/bbc/tools/kin_miniconda2/miniconda2/etc/profile.d/conda.sh
#condame
conda activate snakemake

# make temp directory for tools that need it
if [ ! -d "./temp/" ]
then
    mkdir ./temp/
fi

# irecotry for storing the job o and e files
if [ ! -d "./error_files/" ]
then
    mkdir ./error_files/
fi

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
if [ ! -d "./logs/" ]
then
    mkdir -p ./logs/
fi
snakemake --use-conda -n > logs/prep_bbc_shared_workflow_${TIME}.txt
snakemake --dag | dot -Tpng > logs/dag.png
snakemake --filegraph | dot -Tpng > logs/filegraph.png
snakemake --rulegraph | dot -Tpng > logs/rulegraph.png

# run snakemake
snakemake \
--use-envmodules \
--jobs 100 \
--cluster "qsub \
-q bbc \
-l nodes=1:ppn={threads} \
-l mem={resources.mem_mb}mb \
-l walltime=100:00:00 \
-o error_files/ \
-e error_files/"

