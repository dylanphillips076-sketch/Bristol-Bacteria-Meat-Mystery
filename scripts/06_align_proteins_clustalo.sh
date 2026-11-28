#!/bin/bash
#SBATCH --job-name=meat_protein_alignment
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:0:00
#SBATCH --mem=1000M
#SBATCH --account=bisc040764

# =============================================================================
# 06_align_proteins_clustalo.sh
# Step 6 (part 2): align the representative proteins, references and outgroup
# with Clustal Omega on Blue Pebble.
# Submit on Blue Pebble with:  sbatch ~/job_data/06_align_proteins_clustalo.sh
# Input : ~/job_data/tree_input_proteins.fasta     (from 06_prepare_tree_input.py)
# Output: ~/job_data/tree_proteins_aligned.fasta   (aligned FASTA)
# =============================================================================

module load clustal-omega/1.2.4-skoy

echo "Protein alignment job started"
clustalo -i ~/job_data/tree_input_proteins.fasta \
         -o ~/job_data/tree_proteins_aligned.fasta \
         --seqtype=Protein --outfmt=fasta --force

# record some potentially useful details about the job:
echo "Clustal Omega version: $(clustalo --version)"
echo "Running on host $(hostname)"
echo "Slurm job ID is ${SLURM_JOBID}"
echo "This job ran on the following machines:"
echo "${SLURM_JOB_NODELIST}"
echo "Ended on: $(date)"