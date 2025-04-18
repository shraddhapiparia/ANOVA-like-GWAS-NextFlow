#!/bin/bash -l
#SBATCH --job-name=gwas_nextflow
#SBATCH --output=path/gwas_output_%j.txt
#SBATCH --error=path/gwas_error_%j.txt
#SBATCH --mem=x
#SBATCH --cpus-per-task=x
#SBATCH --mail-user=emailid
#SBATCH --mail-type=ALL
#SBATCH --partition=partition_name
#SBATCH --nodelist=node_name

# Set temporary directory (optional)
export TMPDIR="/d/tmp/username"

# Load modules
source /etc/profile.d/modules.sh
module purge
module load java/15
module load R/4.3.2
module load nextflow 
# Navigate to the folder with your main.nf and config
cd /path/Nextflow

# Run Nextflow
nextflow run new_pipeline.nf -params-file params-file.yaml -resume
