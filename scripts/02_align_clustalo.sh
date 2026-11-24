#!/bin/bash
#SBATCH --job-name=meat_alignment
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:0:00
#SBATCH --mem=1000M
#SBATCH --account=bisc040764

# =============================================================================
# 02_align_clustalo.sh
# Align all 100 labelled sequences with Clustal Omega on Blue Pebble.
# Submit on Blue Pebble with:  sbatch ~/job_data/02_align_clustalo.sh
# Input : ~/job_data/all_sequences.fasta   (made by 01_prepare_data.sh)
# Output: ~/job_data/all_sequences_aligned.fasta  (aligned FASTA)
# =============================================================================

module load clustal-omega/1.2.4-skoy

echo "Alignment job started"
clustalo -i ~/job_data/all_sequences.fasta -o ~/job_data/all_sequences_aligned.fasta --outfmt=fasta --force

# record some potentially useful details about the job:
echo "Clustal Omega version: $(clustalo --version)"
echo "Running on host $(hostname)"
echo "Directory is $(pwd)"
echo "Slurm job ID is ${SLURM_JOBID}"
echo "This job ran on the following machines:"
echo "${SLURM_JOB_NODELIST}"
echo "Ended on: $(date)"
