#!/bin/bash
#SBATCH --job-name=ET_callpeaks_macs2
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=24:00:00
#SBATCH --output=../callpeaks_macs2.%j.out
#SBATCH --error=../callpeaks_macs2.%j.err

#load modules
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4 BWA/0.7.17-GCC-8.3.0 SAMtools/1.16.1-GCC-11.3.0 MACS2/2.2.7.1-foss-2019b-Python-3.7.4

#Make directory where you can access your bam files AND make a subfolder for macs2 output
OUTDIR="/scratch/evt82290/Run133"
OUTDIR2="scratch/evt82290/Run131"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

##Run callpeak on all your samples
macs2 callpeak -t $OUTDIR/SortedBamFiles/133-43_ChIP_WT_H3K27me3_24hr_Rep1_S40_L001_R1_001_val_1.fq.gz.bam -c $OUTDIR2/SortedBamFiles/129-43_ChIP_WT_input_S42_L001_R1_001_val_1.fq.gz.bam -f BAM -n WT_H3K27me3 --outdir $OUTDIR/macs2
macs2 callpeak -t $OUTDIR/SortedBamFiles/133-44_ChIP_cac-1_H3K27me3_24hr_Rep1_S41_L001_R1_001_val_1.fq.gz.bam -c $OUTDIR2/SortedBamFiles/129-44_ChIP_cac-1_input_S43_L001_R1_001_val_1.fq.gz.bam -f BAM -n cac-1_H3K27me3 --outdir $OUTDIR/macs2
macs2 callpeak -t $OUTDIR/SortedBamFiles/133-45_ChIP_cac-2_H3K27me3_24hr_Rep1_S42_L001_R1_001_val_1.fq.gz.bam -c $OUTDIR2/SortedBamFiles/129-45_ChIP_cac-2_input_S44_L001_R1_001_val_1.fq.gz.bam -f BAM -n cac-2_H3K27me3 --outdir $OUTDIR/macs2
macs2 callpeak -t $OUTDIR/SortedBamFiles/133-46_ChIP_set-7_H3K27me3_24hr_Rep1_S43_L001_R1_001_val_1.fq.gz.bam -c $OUTDIR2/SortedBamFiles/129-46_ChIP_cac-3_input_S45_L001_R1_001_val_1.fq.gz.bam -f BAM -n cac-3_H3K27me3 --outdir $OUTDIR/macs2



# 133-43_ChIP_WT_H3K27me3_24hr_Rep1_S40_L001_R1_001_val_1.fq.gz.bam

# 133-44_ChIP_cac-1_H3K27me3_24hr_Rep1_S41_L001_R1_001_val_1.fq.gz.bam

# 133-45_ChIP_cac-2_H3K27me3_24hr_Rep1_S42_L001_R1_001_val_1.fq.gz.bam

# 133-46_ChIP_set-7_H3K27me3_24hr_Rep1_S43_L001_R1_001_val_1.fq.gz.bam

#
#

# 129-33_ChIP_WT_K27me3_AM_Rep_2_S32_L001_R1_001_val_1.fq.gz.bam
# 129-33_ChIP_WT_K27me3_AM_Rep_2_S32_L001_R1_001_val_1.fq.gz.bam.bai
# 129-34_ChIP_cac-1_K27me3_AM_Rep_2_S33_L001_R1_001_val_1.fq.gz.bam
# 129-34_ChIP_cac-1_K27me3_AM_Rep_2_S33_L001_R1_001_val_1.fq.gz.bam.bai
# 129-35_ChIP_cac-2_K27me3_AM_Rep_2_S34_L001_R1_001_val_1.fq.gz.bam
# 129-35_ChIP_cac-2_K27me3_AM_Rep_2_S34_L001_R1_001_val_1.fq.gz.bam.bai
# 129-36_ChIP_cac-3_K27me3_AM_Rep_2_S35_L001_R1_001_val_1.fq.gz.bam
# 129-36_ChIP_cac-3_K27me3_AM_Rep_2_S35_L001_R1_001_val_1.fq.gz.bam.bai
# 129-37_ChIP_set-7_K27me3_AM_Rep_2_S36_L001_R1_001_val_1.fq.gz.bam
# 129-37_ChIP_set-7_K27me3_AM_Rep_2_S36_L001_R1_001_val_1.fq.gz.bam.bai
# 129-38_ChIP_WT_K27me3_AbC_Rep_1_S37_L001_R1_001_val_1.fq.gz.bam
# 129-38_ChIP_WT_K27me3_AbC_Rep_1_S37_L001_R1_001_val_1.fq.gz.bam.bai
# 129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38_L001_R1_001_val_1.fq.gz.bam
# 129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38_L001_R1_001_val_1.fq.gz.bam.bai
# 129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39_L001_R1_001_val_1.fq.gz.bam
# 129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39_L001_R1_001_val_1.fq.gz.bam.bai
# 129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40_L001_R1_001_val_1.fq.gz.bam
# 129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40_L001_R1_001_val_1.fq.gz.bam.bai
# 129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41_L001_R1_001_val_1.fq.gz.bam
# 129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41_L001_R1_001_val_1.fq.gz.bam.bai
# 129-43_ChIP_WT_input_S42_L001_R1_001_val_1.fq.gz.bam
# 129-43_ChIP_WT_input_S42_L001_R1_001_val_1.fq.gz.bam.bai
# 129-44_ChIP_cac-1_input_S43_L001_R1_001_val_1.fq.gz.bam
# 129-44_ChIP_cac-1_input_S43_L001_R1_001_val_1.fq.gz.bam.bai
# 129-45_ChIP_cac-2_input_S44_L001_R1_001_val_1.fq.gz.bam
# 129-45_ChIP_cac-2_input_S44_L001_R1_001_val_1.fq.gz.bam.bai
# 129-46_ChIP_cac-3_input_S45_L001_R1_001_val_1.fq.gz.bam
# 129-46_ChIP_cac-3_input_S45_L001_R1_001_val_1.fq.gz.bam.bai
# 129-47_ChIP_set-7_input_S46_L001_R1_001_val_1.fq.gz.bam
# 129-47_ChIP_set-7_input_S46_L001_R1_001_val_1.fq.gz.bam.bai
# 129-48_ChIP-seq_Dim-5_input_S47_L001_R1_001_val_1.fq.gz.bam
# 129-48_ChIP-seq_Dim-5_input_S47_L001_R1_001_val_1.fq.gz.bam.bai
# 129-90_ChIP_WT_H3K36me3_Rep1_S71_L001_R1_001_val_1.fq.gz.bam
# 129-90_ChIP_WT_H3K36me3_Rep1_S71_L001_R1_001_val_1.fq.gz.bam.bai
# 129-91_ChIP_cac-1_H3K36me3_Rep1_S72_L001_R1_001_val_1.fq.gz.bam
# 129-91_ChIP_cac-1_H3K36me3_Rep1_S72_L001_R1_001_val_1.fq.gz.bam.bai
# 129-92_ChIP_cac-2_H3K36me3_Rep1_S73_L001_R1_001_val_1.fq.gz.bam
# 129-92_ChIP_cac-2_H3K36me3_Rep1_S73_L001_R1_001_val_1.fq.gz.bam.bai
# 129-93_ChIP_cac-3_H3K36me3_Rep1_S74_L001_R1_001_val_1.fq.gz.bam
# 129-93_ChIP_cac-3_H3K36me3_Rep1_S74_L001_R1_001_val_1.fq.gz.bam.bai
# 129-94_ChIP_set-7_H3K36me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bam
# 129-94_ChIP_set-7_H3K36me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bam.bai
