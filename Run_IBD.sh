#!/bin/bash
#SBATCH --job-name=IBD_segment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --out=Run_IBD_%j.out

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"

module load nextflow

export NXF_OFFLINE=true

CONFIG_FILE="configs/YOUR_CONFIG_FILE"
WORK_DIRECTORY="work/WORKING_DIRECTORY"

echo "Configuration: ${CONFIG_FILE}"
echo "Work directory: ${WORK_DIRECTORY}"

nextflow run IBD_Pipeline.nf \
    -c "${CONFIG_FILE}" \
    -work-dir "${WORK_DIRECTORY}" \
    -resume
