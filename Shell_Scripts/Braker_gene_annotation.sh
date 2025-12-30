#!/bin/bash
#SBATCH --job-name=ET_Braker.%j.job
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../Braker.%j.out
#SBATCH --error=../Braker.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($GENOME, fastqpath)

source config.txt

#load modules
module load BRAKER/3.0.8-foss-2022a
module load GeneMark-ET/4.72-GCCcore-12.3.0
module load AUGUSTUS
module load SAMtools

# #extract unannotated gene regions from genome
# samtools faidx "$GENOME" CM002241.1:1931560-1965456 > subregion_genome.fasta
#
# #index bam
# samtools index /scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles/SRR7970598/SRR7970598_Aligned.sortedByCoord.out.bam

#set AUGUSTUS directory w/ write permissions
export AUGUSTUS_CONFIG_PATH=/scratch/evt82290/BRAKER/augustus_config

#Run BRAKER
braker.pl --genome /home/evt82290/Research/GCA_000182925.2_NC12_clean.fna --bam subregion.bam --softmasking --gff3 --species Neurospora_crassa_cac3_new --workingdir /scratch/evt82290/BRAKER
