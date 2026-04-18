#!/bin/bash -l
# TEMPLATE: Copy this file and fill in all <PLACEHOLDER> values before submitting.
# Do not submit this file directly.

#SBATCH --job-name=gwas_nextflow
#SBATCH --output=<PATH_TO_LOG_DIR>/gwas_output_%j.txt
#SBATCH --error=<PATH_TO_LOG_DIR>/gwas_error_%j.txt
#SBATCH --mem=<MEMORY_GB>G
#SBATCH --cpus-per-task=<NUM_CPUS>
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL
#SBATCH --partition=<PARTITION_NAME>
#SBATCH --nodelist=<NODE_NAME>   # Remove this line if you don't need a specific node

# Set temporary directory (required on most HPC systems)
export TMPDIR="<PATH_TO_TMP_DIR>"

# Load modules (adjust module names/versions for your cluster)
source /etc/profile.d/modules.sh
module purge
module load java/17
module load R/4.3.2
module load nextflow

# Navigate to the workflow directory
cd <PATH_TO_REPO>/workflow

# Run Nextflow
nextflow run new_pipeline.nf -params-file params-file.yaml -resume
