#!/bin/bash
#SBATCH --job-name=IBD_segment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00

module load nextflow
nextflow run /path/to/IBD_pipeline.nf -w /path/to/working/directory
