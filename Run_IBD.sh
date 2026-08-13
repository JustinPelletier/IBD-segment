#!/bin/bash
#SBATCH --job-name=IBD_segment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --out=Run_IBD_%j%.out

module load nextflow
nextflow run IBD_Pipeline.nf
#nextflow run IBD_Pipeline.nf -resume

