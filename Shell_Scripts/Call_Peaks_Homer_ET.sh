#!/bin/bash
#SBATCH --job-name=Call_Peaks_Homer_ET
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../Homer/logs/CallPeak.%j.out
#SBATCH --error=../Homer/logs/CallPeak.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

#HOMER

#source config.txt

#Make  Directories
OUTDIR="/scratch/evt82290/Peaks/qa-suz12"
P126DIR="/scratch/evt82290/MappingOutputs/Run126/bamFiles"
P129DIR="/scratch/evt82290/MappingOutputs/Run129/bamFiles"
P131DIR="/scratch/evt82290/MappingOutputs/Run131/bamFiles"
P133DIR="/scratch/evt82290/MappingOutputs/Run133/bamFiles"
P136DIR="/scratch/evt82290/MappingOutputs/Run136/bamFiles"
P137DIR="/scratch/evt82290/MappingOutputs/Run137/bamFiles"
P138DIR="/scratch/evt82290/MappingOutputs/Run138/bamFiles"
P139DIR="/scratch/evt82290/MappingOutputs/Run139/bamFiles"
P141DIR="/scratch/evt82290/MappingOutputs/Run141/bamFiles"
P144DIR="/scratch/evt82290/MappingOutputs/Run144/bamFiles"
P145DIR="/scratch/evt82290/MappingOutputs/Run145/bamFiles"
P146DIR="/scratch/evt82290/MappingOutputs/Run146/bamFiles"
P148DIR="/scratch/evt82290/MappingOutputs/Run148/bamFiles"
TAGDIR="/scratch/evt82290/Peaks/qa-suz12/homertags"
#PEAKDIR="/scratch/evt82290/Run137/Peaks"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
###

ml SAMtools
ml BWA
# #
# bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIRtempReps $OUTDIR "$bam" -
# samtools index "$bam"
#
# #samtools view -b -q 30 $bam > "$QualityBam"
# #samtools index "$QualityBam"
#
# ############################
# # # #deeptools

ml deepTools

#HOMER
module load Homer/4.11-foss-2022a

#Make Tag Directories

###QA-SUZ12###
#H3K27me3 24 hr
#Rep 1 qa-suz12
# makeTagDirectory $TAGDIR/WT_0hr_rep1 $P137DIR/137-66_ChIP_WT_H3K27me3_Rep1_S63.bam
# makeTagDirectory $TAGDIR/qasuz12_0hr_rep1 $P137DIR/137-67_ChIP_qa-suz12_H3K27me3_Rep1_S64.bam
# makeTagDirectory $TAGDIR/qasuz12_4hr_rep1 $P137DIR/137-69_ChIP_qa-suz12_H3K27me3_Rep1_S66.bam
# makeTagDirectory $TAGDIR/qasuz12_8hr_rep1 $P137DIR/137-71_ChIP_qa-suz12_H3K27me3_Rep1_S68.bam
# makeTagDirectory $TAGDIR/qasuz12_12hr_rep1 $P137DIR/137-73_ChIP_qa-suz12_H3K27me3_Rep1_S70.bam
# makeTagDirectory $TAGDIR/qasuz12_24hr_rep1 $P137DIR/137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72.bam
# makeTagDirectory $TAGDIR/WT_24hr_rep1 $P137DIR/137-74_ChIP_WT_H3K27me3_Rep1_S71.bam
# makeTagDirectory $TAGDIR/WT_0hr_input_rep1 $P139DIR/139-29_ChIP_WT_Input__S29.bam
# makeTagDirectory $TAGDIR/qasuz12_0hr_input_rep1 $P139DIR/139-30_ChIP_qa-suz12_Input__S30.bam
# makeTagDirectory $TAGDIR/qasuz12_4hr_input_rep1 $P139DIR/139-32_ChIP_qa-suz12_Input__S32.bam
# makeTagDirectory $TAGDIR/qasuz12_8hr_input_rep1 $P139DIR/139-34_ChIP_qa-suz12_Input__S34.bam
# makeTagDirectory $TAGDIR/qasuz12_12hr_input_rep1 $P139DIR/139-36_ChIP_qa-suz12_Input__S36.bam
# makeTagDirectory $TAGDIR/qasuz12_24hr_input_rep1 $P139DIR/139-38_ChIP_qa-suz12_Input__S38.bam
# makeTagDirectory $TAGDIR/WT_24hr_input_rep1 $P139DIR/139-37_ChIP_WT_Input__S37.bam

# #Rep 2 qa-suz12
# makeTagDirectory $TAGDIR/WT_0hr_rep2 $P144DIR/144-44_ChIP_WT_0hr_H3K27me3_Rep4_S44.bam
# makeTagDirectory $TAGDIR/qasuz12_0hr_rep2 $P144DIR/144-48_ChIP_qa-suz12_0hr_H3K27me3_Rep4_S48.bam
# makeTagDirectory $TAGDIR/qasuz12_4hr_rep2 $P144DIR/144-49_ChIP_qa-suz12_4hr_H3K27me3_Rep3_S49.bam
# makeTagDirectory $TAGDIR/qasuz12_8hr_rep2 $P144DIR/144-50_ChIP_qa-suz12_8hr_H3K27me3_Rep3_S50.bam
# makeTagDirectory $TAGDIR/qasuz12_12hr_rep2 $P144DIR/144-51_ChIP_qa-suz12_12hr_H3K27me3_Rep3_S51.bam
# makeTagDirectory $TAGDIR/qasuz12_24hr_rep2 $P144DIR/144-52_ChIP_qa-suz12_24hr_H3K27me3_Rep4_S52.bam
# makeTagDirectory $TAGDIR/WT_24hr_rep2 $P144DIR/144-45_ChIP_WT_24hr_H3K27me3_Rep3_S45.bam
# makeTagDirectory $TAGDIR/WT_0hr_input_rep2 $P144DIR/144-62_ChIP_WT_Input__S62.bam
# makeTagDirectory $TAGDIR/qasuz12_0hr_input_rep2 $P144DIR/144-71_ChIP_qa-suz12_Input__S71.bam
# makeTagDirectory $TAGDIR/qasuz12_4hr_input_rep2 $P144DIR/144-73_ChIP_qa-suz12_Input__S73.bam
# makeTagDirectory $TAGDIR/qasuz12_8hr_input_rep2 $P144DIR/144-79_ChIP_qa-suz12_Input__S79.bam
# #makeTagDirectory $TAGDIR/qasuz12_12hr_input_rep2 $P144DIR/
# makeTagDirectory $TAGDIR/qasuz12_24hr_input_rep2 $P144DIR/144-132_ChIP_qa-suz12_Input__S132.bam
# makeTagDirectory $TAGDIR/WT_24hr_input_rep2 $INPDIR/

# #Rep 3 qa-suz12
# makeTagDirectory $TAGDIR/WT_0hr_rep3 $P146DIR/
# makeTagDirectory $TAGDIR/qasuz12_0hr_rep3 $P146DIR/146-125_ChIP_qa-suz12_0hr_H3K27me3_Rep1_S144.bam
# makeTagDirectory $TAGDIR/qasuz12_4hr_rep3 $P146DIR/146-126_ChIP_qa-suz12_4hr_H3K27me3_Rep1_S145.bam
# makeTagDirectory $TAGDIR/qasuz12_8hr_rep3 $P146DIR/146-127_ChIP_qa-suz12_8hr_H3K27me3_Rep1_S146.bam
# makeTagDirectory $TAGDIR/qasuz12_12hr_rep3 $P146DIR/146-128_ChIP_qa-suz12_12hr_H3K27me3_Rep1_S147.bam
# makeTagDirectory $TAGDIR/qasuz12_24hr_rep3 $P146DIR/146-129_ChIP_qa-suz12_24hr_H3K27me3_Rep1_S148.bam
# makeTagDirectory $TAGDIR/WT_24hr_rep3 $P146DIR/146-6_ChIP_WT_H3K27me3_Rep1_S6.bam
# makeTagDirectory $TAGDIR/WT_0hr_input_rep3 $P139DIR/
# makeTagDirectory $TAGDIR/qasuz12_0hr_input_rep3 $P139DIR/
# makeTagDirectory $TAGDIR/qasuz12_4hr_input_rep3 $P139DIR/
# makeTagDirectory $TAGDIR/qasuz12_8hr_input_rep3 $P139DIR/
# makeTagDirectory $TAGDIR/qasuz12_12hr_input_rep3 $P139DIR/
# makeTagDirectory $TAGDIR/qasuz12_24hr_input_rep3 $P139DIR/
# makeTagDirectory $TAGDIR/WT_24hr_input_rep3 $P146DIR/146-18_ChIP_WT_input__S18.bam

#H3K27me3 96 hr
#Rep 1
makeTagDirectory $TAGDIR/qasuz12_0hr_rep4 $P146DIR/146-39_ChIP_qa-suz12_0hr_H3K27me3_Rep1_S39.bam
makeTagDirectory $TAGDIR/qasuz12_24hr_rep4 $P146DIR/146-40_ChIP_qa-suz12_24hr_H3K27me3_Rep1_S40.bam
makeTagDirectory $TAGDIR/qasuz12_48hr_rep1 $P146DIR/146-41_ChIP_qa-suz12_48hr_H3K27me3_Rep1_S41.bam
makeTagDirectory $TAGDIR/qasuz12_72hr_rep1 $P146DIR/146-42_ChIP_qa-suz12_72hr_H3K27me3_Rep1_S42.bam
makeTagDirectory $TAGDIR/qasuz12_96hr_rep1 $P146DIR/146-43_ChIP_qa-suz12_96hr_H3K27me3_Rep1_S43.bam
makeTagDirectory $TAGDIR/WT_96hr_rep1 $P146DIR/146-44_ChIP_WT_H3K27me3_Rep1_S44.bam

# 146-39_ChIP_qa-suz12_0hr_H3K27me3_Rep1_S39.bam
# 146-40_ChIP_qa-suz12_24hr_H3K27me3_Rep1_S40.bam
# 146-41_ChIP_qa-suz12_48hr_H3K27me3_Rep1_S41.bam
# 146-42_ChIP_qa-suz12_72hr_H3K27me3_Rep1_S42.bam
# 146-43_ChIP_qa-suz12_96hr_H3K27me3_Rep1_S43.bam
# 146-44_ChIP_WT_H3K27me3_Rep1_S44.bam

###CAF-1###

#H3K27me3 abcam
# makeTagDirectory $TAGDIR/WT_H3K27me3_abc_rep1 $P129DIR/129-38_ChIP_WT_K27me3_AbC_Rep_1_S37.bam
# makeTagDirectory $TAGDIR/WT_H3K27me3_abc_rep2 $P136DIR/6147_136-1_ChIP_WT_H3K27me3_abcam_Rep2_S1.bam
# # makeTagDirectory $TAGDIR/WT_H3K27me3_abc_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/cac-1_H3K27me3_abc_rep1 $P129DIR/129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38.bam
# makeTagDirectory $TAGDIR/cac-1_H3K27me3_abc_rep2 $P136DIR/6147_136-2_ChIP_cac-1_H3K27me3_abcam_Rep2_S2.bam
# # makeTagDirectory $TAGDIR/cac-1_H3K27me3_abc_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/cac-2_H3K27me3_abc_rep1 $P129DIR/129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39.bam
# makeTagDirectory $TAGDIR/cac-2_H3K27me3_abc_rep2 $P136DIR/6147_136-3_ChIP_cac-2_H3K27me3_abcam_Rep2_S3.bam
# # makeTagDirectory $TAGDIR/cac-2_H3K27me3_abc_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/cac-3_H3K27me3_abc_rep1 $P129DIR/129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40.bam
# makeTagDirectory $TAGDIR/cac-3_H3K27me3_abc_rep2 $P136DIR/6147_136-4_ChIP_cac-3_H3K27me3_abcam_Rep2_S4.bam
# # makeTagDirectory $TAGDIR/cac-3_H3K27me3_abc_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/set-7_H3K27me3_abc_rep1 $P129DIR/129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41.bam
# makeTagDirectory $TAGDIR/set-7_H3K27me3_abc_rep2 $P136DIR/6147_136-77_ChIP_set-7_H3K27me3_abcam_Rep3_S76.bam
# # makeTagDirectory $TAGDIR/set-7_H3K27me3_abc_rep3 $138DIR/
#
# #H3K27me3 Cell Signaling
# makeTagDirectory $TAGDIR/WT_H3K27me3_CS_rep1 $P136DIR/6147_136-6_ChIP_WT_H3K27me3_CS_Rep1_S6.bam
# makeTagDirectory $TAGDIR/WT_H3K27me3_CS_rep2 $P136DIR/6147_136-78_ChIP_WT_H3K27me3_CS_Rep2_S77.bam
# makeTagDirectory $TAGDIR/WT_H3K27me3_CS_rep3 $P138DIR/138-57_ChIP_WT_H3K27me3_Rep3_6252_S56.bam
#
# makeTagDirectory $TAGDIR/cac-1_H3K27me3_CS_rep1 $P136DIR/6147_136-7_ChIP_cac-1_H3K27me3_CS_Rep1_S7.bam
# makeTagDirectory $TAGDIR/cac-1_H3K27me3_CS_rep2 $P136DIR/6147_136-79_ChIP_cac-1_H3K27me3_CS_Rep2_S78.bam
# makeTagDirectory $TAGDIR/cac-1_H3K27me3_CS_rep3 $P138DIR/138-58_ChIP_cac-1_H3K27me3_Rep3_6252_S57.bam
#
# makeTagDirectory $TAGDIR/cac-2_H3K27me3_CS_rep1 $P136DIR/6147_136-8_ChIP_cac-2_H3K27me3_CS_Rep1_S8.bam
# makeTagDirectory $TAGDIR/cac-2_H3K27me3_CS_rep2 $P136DIR/6147_136-80_ChIP_cac-2_H3K27me3_CS_Rep2_S79.bam
# makeTagDirectory $TAGDIR/cac-2_H3K27me3_CS_rep3 $P138DIR/138-59_ChIP_cac-2_H3K27me3_Rep3_6252_S58.bam
#
# makeTagDirectory $TAGDIR/cac-3_H3K27me3_CS_rep1 $P136DIR/6147_136-9_ChIP_cac-3_H3K27me3_CS_Rep1_S9.bam
# makeTagDirectory $TAGDIR/cac-3_H3K27me3_CS_rep2 $P136DIR/6147_136-81_ChIP_cac-3_H3K27me3_CS_Rep2_S80.bam
# makeTagDirectory $TAGDIR/cac-3_H3K27me3_CS_rep3 $P138DIR/138-60_ChIP_cac-3_H3K27me3_Rep3_6252_S59.bam

# makeTagDirectory $TAGDIR/set-7_H3K27me3_CS_rep1 $P126DIR/126-58_ChIP_set-7_H3K27me3_CS_Rep3_S52.bam
# makeTagDirectory $TAGDIR/set-7_H3K27me3_CS_rep2 $P136DIR/6147_136-83_ChIP_set-7_H3K27me3_CS_Rep2_S82.bam
# makeTagDirectory $TAGDIR/set-7_H3K27me3_CS_rep3 $P138DIR/138-61_ChIP_set-7_H3K27me3_Rep3_6252_S60.bam

#H3K27me2
# makeTagDirectory $TAGDIR/WT_H3K27me2_rep1 $P141DIR/141-88_ChIP_WT_H3K27me2_Rep1_S82.bam
# # makeTagDirectory $TAGDIR/WT_H3K27me2_rep2 $141DIR/
# # makeTagDirectory $TAGDIR/WT_H3K27me2_rep3 $141DIR/
#
# makeTagDirectory $TAGDIR/cac-1_H3K27me2_rep1 $P141DIR/141-89_ChIP_cac-1_H3K27me2_Rep1_S83.bam
# # makeTagDirectory $TAGDIR/cac-1_H3K27me2_rep2 $141DIR/
# # makeTagDirectory $TAGDIR/cac-1_H3K27me2_rep3 $141DIR/
#
# makeTagDirectory $TAGDIR/cac-2_H3K27me2_rep1 $P141DIR/141-90_ChIP_cac-2_H3K27me2_Rep1_S84.bam
# # makeTagDirectory $TAGDIR/cac-2_H3K27me2_rep2 $141DIR/
# # makeTagDirectory $TAGDIR/cac-2_H3K27me2_rep3 $141DIR/
#
# makeTagDirectory $TAGDIR/cac-3_H3K27me2_rep1 $P141DIR/141-91_ChIP_cac-3_H3K27me2_Rep1_S85.bam
# # makeTagDirectory $TAGDIR/cac-3_H3K27me2_rep2 $141DIR/
# # makeTagDirectory $TAGDIR/cac-3_H3K27me2_rep3 $141DIR/
#
# makeTagDirectory $TAGDIR/set-7_H3K27me2_rep1 $P141DIR/141-92_ChIP_set-7_H3K27me2_Rep1_S86.bam
# # makeTagDirectory $TAGDIR/set-7_H3K27me2_rep2 $138DIR/
# # makeTagDirectory $TAGDIR/set-7_H3K27me2_rep3 $138DIR/
#
# #H3K36me3
# makeTagDirectory $TAGDIR/WT_H3K36me3_rep1 $P129DIR/129-90_ChIP_WT_H3K36me3_Rep1_S71.bam
# makeTagDirectory $TAGDIR/WT_H3K36me3_rep2 $P138DIR/138-62_ChIP_WT_H3K36me3_Rep2_6252_S61.bam
# # makeTagDirectory $TAGDIR/WT_H3K36me3_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/cac-1_H3K36me3_rep1 $P129DIR/129-91_ChIP_cac-1_H3K36me3_Rep1_S72.bam
# makeTagDirectory $TAGDIR/cac-1_H3K36me3_rep2 $P138DIR/138-63_ChIP_cac-1_H3K36me3_Rep2_6252_S62.bam
# # makeTagDirectory $TAGDIR/cac-1_H3K36me3_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/cac-2_H3K36me3_rep1 $P129DIR/129-92_ChIP_cac-2_H3K36me3_Rep1_S73.bam
# makeTagDirectory $TAGDIR/cac-2_H3K36me3_rep2 $P138DIR/138-64_ChIP_cac-2_H3K36me3_Rep2_6252_S63.bam
# # makeTagDirectory $TAGDIR/cac-2_H3K36me3_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/cac-3_H3K36me3_rep1 $P129DIR/129-93_ChIP_cac-3_H3K36me3_Rep1_S74.bam
# makeTagDirectory $TAGDIR/cac-3_H3K36me3_rep2 $P138DIR/138-65_ChIP_cac-3_H3K36me3_Rep2_6252_S64.bam
# # makeTagDirectory $TAGDIR/cac-3_H3K36me3_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/set-7_H3K36me3_rep1 $P129DIR/129-94_ChIP_set-7_H3K36me3_Rep1_S75.bam
# makeTagDirectory $TAGDIR/set-7_H3K36me3_rep2 $P138DIR/138-66_ChIP_set-7_H3K36me3_Rep2_6252_S65.bam
# # makeTagDirectory $TAGDIR/set-7_H3K36me3_rep3 $138DIR/
#
# #H3K9me3
# makeTagDirectory $TAGDIR/WT_H3K9me3_rep1 $P133DIR/133-23_ChIP_WT_H3K9me3_Rep2_S21.bam
# # makeTagDirectory $TAGDIR/WT_H3K9me3_rep2 $133DIR/
# # makeTagDirectory $TAGDIR/WT_H3K9me3_rep3 $133DIR/
#
# makeTagDirectory $TAGDIR/cac-1_H3K9me3_rep1 $P133DIR/133-24_ChIP_cac-1_H3K9me3_Rep2_S22.bam
# # makeTagDirectory $TAGDIR/cac-1_H3K9me3_rep2 $133DIR/
# # makeTagDirectory $TAGDIR/cac-1_H3K9me3_rep3 $133DIR/
#
# makeTagDirectory $TAGDIR/cac-2_H3K9me3_rep1 $P133DIR/133-25_ChIP_cac-2_H3K9me3_Rep2_S23.bam
# # makeTagDirectory $TAGDIR/cac-2_H3K9me3_rep2 $133DIR/
# # makeTagDirectory $TAGDIR/cac-2_H3K9me3_rep3 $133DIR/
#
# # makeTagDirectory $TAGDIR/cac-3_H3K9me3_rep1 $133DIR/
# # makeTagDirectory $TAGDIR/cac-3_H3K9me3_rep2 $133DIR/
# # makeTagDirectory $TAGDIR/cac-3_H3K9me3_rep3 $133DIR/
#
# makeTagDirectory $TAGDIR/set-7_H3K9me3_rep1 $P133DIR/133-26_ChIP_set-7_H3K9me3_Rep2_S24.bam
# # makeTagDirectory $TAGDIR/set-7_H3K9me3_rep2 $138DIR/
# # makeTagDirectory $TAGDIR/set-7_H3K9me3_rep3 $138DIR/

# #H3K4me2
# makeTagDirectory $TAGDIR/WT_H3K4me2_rep1 $P146DIR/146-19_ChIP_WT_H3K4me2_Rep1_S19.bam
# makeTagDirectory $TAGDIR/WT_H3K4me2_rep2 $P146DIR/146-29_ChIP_WT_H3K4me2_Rep2_S29.bam
# makeTagDirectory $TAGDIR/WT_H3K4me2_rep3 $146DIR/
#
# makeTagDirectory $TAGDIR/cac-1_H3K4me2_rep1 $P146DIR/146-20_ChIP_cac-1_H3K4me2_Rep1_S20.bam
# makeTagDirectory $TAGDIR/cac-1_H3K4me2_rep2 $P146DIR/146-30_ChIP_cac-1_H3K4me2_Rep2_S30.bam
# makeTagDirectory $TAGDIR/cac-1_H3K4me2_rep3 $146DIR/
#
# makeTagDirectory $TAGDIR/cac-2_H3K4me2_rep1 $P146DIR/146-21_ChIP_cac-2_H3K4me2_Rep1_S21.bam
# makeTagDirectory $TAGDIR/cac-2_H3K4me2_rep2 $P146DIR/146-31_ChIP_cac-2_H3K4me2_Rep2_S31.bam
# makeTagDirectory $TAGDIR/cac-2_H3K4me2_rep3 $146DIR/
#
# makeTagDirectory $TAGDIR/cac-3_H3K4me2_rep1 $P146DIR/146-22_ChIP_cac-3_H3K4me2_Rep1_S22.bam
# makeTagDirectory $TAGDIR/cac-3_H3K4me2_rep2 $P146DIR/146-32_ChIP_cac-3_H3K4me2_Rep2_S32.bam
# makeTagDirectory $TAGDIR/cac-3_H3K4me2_rep3 $146DIR/
#
# makeTagDirectory $TAGDIR/set-7_H3K4me2_rep1 $P146DIR/146-23_ChIP_set-7_H3K4me2_Rep1_S23.bam
# makeTagDirectory $TAGDIR/set-7_H3K4me2_rep2 $P146DIR/146-33_ChIP_set-7_H3K4me2_Rep2_S33.bam
# makeTagDirectory $TAGDIR/set-7_H3K4me2_rep3 $146DIR/
#
# #H4K20me3
# makeTagDirectory $TAGDIR/WT_H4K20me3_rep1 $P133DIR/133-19_ChIP_WT_H4K20me3_Rep1_S17.bam
# makeTagDirectory $TAGDIR/WT_H4K20me3_rep2 $P138DIR/138-67_ChIP_WT_H4K20me3_Rep2_6252_S66.bam
# # makeTagDirectory $TAGDIR/WT_H4K20me3_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/cac-1_H4K20me3_rep1 $P133DIR/133-20_ChIP_cac-1_H4K20me3_Rep1_S18.bam
# makeTagDirectory $TAGDIR/cac-1_H4K20me3_rep2 $P138DIR/138-68_ChIP_cac-1_H4K20me3_Rep2_6252_S67.bam
# # makeTagDirectory $TAGDIR/cac-1_H4K20me3_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/cac-2_H4K20me3_rep1 $P133DIR/133-21_ChIP_cac-2_H4K20me3_Rep1_S19.bam
# makeTagDirectory $TAGDIR/cac-2_H4K20me3_rep2 $P138DIR/138-69_ChIP_cac-2_H4K20me3_Rep2_6252_S68.bam
# # makeTagDirectory $TAGDIR/cac-2_H4K20me3_rep3 $138DIR/
#
# # makeTagDirectory $TAGDIR/cac-3_H4K20me3_rep1 $138DIR/
# makeTagDirectory $TAGDIR/cac-3_H4K20me3_rep2 $P138DIR/138-70_ChIP_cac-3_H4K20me3_Rep2_6252_S69.bam
# # makeTagDirectory $TAGDIR/cac-3_H4K20me3_rep3 $138DIR/
#
# makeTagDirectory $TAGDIR/set-7_H4K20me3_rep1 $P133DIR/133-22_ChIP_set-7_H4K20me3_Rep1_S20.bam
# makeTagDirectory $TAGDIR/set-7_H4K20me3_rep2 $P138DIR/138-71_ChIP_set-7_H4K20me3_Rep2_6252_S70.bam
# makeTagDirectory $TAGDIR/set-7_H4K20me3_rep3 $138DIR/

####DO NOT USE RUN131 INPUT
#input
#WT
# makeTagDirectory $TAGDIR/WT_input_1 $P129DIR/129-43_ChIP_WT_input_S42.bam
# makeTagDirectory $TAGDIR/WT_input_2 $P131DIR/131-37_ChIP_WT_input_Rep1_S27.bam
# makeTagDirectory $TAGDIR/WT_input_5 $P138DIR/138-72_ChIP_WT_input__6252_S71.bam
# makeTagDirectory $TAGDIR/WT_input_4 $P136DIR/6147_136-84_ChIP_WT_input_S83.bam
# makeTagDirectory $TAGDIR/WT_input_3 $P136DIR/6147_136-11_ChIP_WT_input_S11.bam
# makeTagDirectory $TAGDIR/WT_input_6 $P146DIR/146-34_ChIP_WT_input__S34.bam
#cac-1
# makeTagDirectory $TAGDIR/cac-1_input_1 $P129DIR/129-44_ChIP_cac-1_input_S43.bam
# makeTagDirectory $TAGDIR/cac-1_input_2 $P131DIR/131-38_ChIP_cac-1_input_Rep1_S28.bam
# makeTagDirectory $TAGDIR/cac-1_input_5 $P138DIR/138-73_ChIP_cac-1_input__6252_S72.bam
# makeTagDirectory $TAGDIR/cac-1_input_4 $P136DIR/6147_136-85_ChIP_cac-1_input_S84.bam
# makeTagDirectory $TAGDIR/cac-1_input_3 $P136DIR/6147_136-12_ChIP_cac-1_input_S12.bam
# makeTagDirectory $TAGDIR/cac-1_input_6 $P146DIR/146-35_ChIP_cac-1_input__S35.bam
#cac-2
# makeTagDirectory $TAGDIR/cac-2_input_1 $P129DIR/129-45_ChIP_cac-2_input_S44.bam
# makeTagDirectory $TAGDIR/cac-2_input_2 $P131DIR/131-39_ChIP_cac-2_input_Rep1_S29.bam
# makeTagDirectory $TAGDIR/cac-2_input_5 $P138DIR/138-74_ChIP_cac-2_input__6252_S73.bam
# makeTagDirectory $TAGDIR/cac-2_input_4 $P136DIR/6147_136-89_ChIP_cac-2_input_S88.bam
# makeTagDirectory $TAGDIR/cac-2_input_3 $P136DIR/6147_136-13_ChIP_cac-2_input_S13.bam
# makeTagDirectory $TAGDIR/cac-2_input_6 $P146DIR/146-36_ChIP_cac-2_input__S36.bam
#cac-3
# makeTagDirectory $TAGDIR/cac-3_input_1 $P129DIR/129-46_ChIP_cac-3_input_S45.bam
# makeTagDirectory $TAGDIR/cac-3_input_2 $P131DIR/131-40_ChIP_cac-3_input_Rep1_S30.bam
# makeTagDirectory $TAGDIR/cac-3_input_4 $P138DIR/138-75_ChIP_cac-3_input__6252_S74.bam
# makeTagDirectory $TAGDIR/cac-3_input_3 $P136DIR/6147_136-14_ChIP_cac-3_input_S14.bam
# makeTagDirectory $TAGDIR/cac-3_input_5 $P146DIR/146-37_ChIP_cac-3_input__S37.bam
#set-7
# makeTagDirectory $TAGDIR/set-7_input_1 $P129DIR/129-47_ChIP_set-7_input_S46.bam
# makeTagDirectory $TAGDIR/set-7_input_2 $P131DIR/131-41_ChIP_set-7_input_Rep1_S31.bam
# makeTagDirectory $TAGDIR/set-7_input_4 $P138DIR/138-76_ChIP_set-7_input__6252_S75.bam
# makeTagDirectory $TAGDIR/set-7_input_3 $P136DIR/6147_136-92_ChIP_set-7_input_S91.bam
# makeTagDirectory $TAGDIR/set-7_input_5 $P146DIR/146-38_ChIP_set-7_input__S38.bam

#call Peaks

#qa-suz12 #24 hr
#rep 1
# findPeaks $TAGDIR/WT_0hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_0hr_input_rep1
# findPeaks $TAGDIR/qasuz12_0hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_0hr_input_rep1
# findPeaks $TAGDIR/qasuz12_4hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_4hr_input_rep1
# findPeaks $TAGDIR/qasuz12_8hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_8hr_input_rep1
# findPeaks $TAGDIR/qasuz12_12hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_12hr_input_rep1
# findPeaks $TAGDIR/qasuz12_24hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_24hr_input_rep1
# findPeaks $TAGDIR/WT_24hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_0hr_input_rep1
# #rep 2
# findPeaks $TAGDIR/WT_0hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_0hr_input_rep2
# findPeaks $TAGDIR/qasuz12_0hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_0hr_input_rep2
# findPeaks $TAGDIR/qasuz12_4hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_4hr_input_rep2
# findPeaks $TAGDIR/qasuz12_8hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_8hr_input_rep2
# findPeaks $TAGDIR/qasuz12_12hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_12hr_input_rep1
# findPeaks $TAGDIR/qasuz12_24hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_24hr_input_rep2
# findPeaks $TAGDIR/WT_24hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_0hr_input_rep1
#rep 3
# findPeaks $TAGDIR/qasuz12_0hr_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_0hr_input_rep1
# findPeaks $TAGDIR/qasuz12_4hr_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_4hr_input_rep1
# findPeaks $TAGDIR/qasuz12_8hr_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_8hr_input_rep1
# findPeaks $TAGDIR/qasuz12_12hr_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_12hr_input_rep1
# findPeaks $TAGDIR/qasuz12_24hr_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_24hr_input_rep1
# findPeaks $TAGDIR/WT_24hr_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_24hr_input_rep3

#qa-suz12 #96 hr

findPeaks $TAGDIR/qasuz12_0hr_rep4 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_0hr_input_rep1
findPeaks $TAGDIR/qasuz12_24hr_rep4 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_4hr_input_rep1
findPeaks $TAGDIR/qasuz12_48hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_8hr_input_rep1
findPeaks $TAGDIR/qasuz12_72hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_12hr_input_rep1
findPeaks $TAGDIR/qasuz12_96hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_24hr_input_rep1
findPeaks $TAGDIR/WT_96hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_24hr_input_rep3

#CAF-1
#H3K27me3

#CS
#rep 1
# findPeaks $TAGDIR/WT_H3K27me3_CS_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_3
# findPeaks $TAGDIR/cac-1_H3K27me3_CS_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_3
# findPeaks $TAGDIR/cac-2_H3K27me3_CS_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_3
# findPeaks $TAGDIR/cac-3_H3K27me3_CS_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_3
# findPeaks $TAGDIR/set-7_H3K27me3_CS_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_3
# #rep 2
# findPeaks $TAGDIR/WT_H3K27me3_CS_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_4
# findPeaks $TAGDIR/cac-1_H3K27me3_CS_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_4
# findPeaks $TAGDIR/cac-2_H3K27me3_CS_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_4
# findPeaks $TAGDIR/cac-3_H3K27me3_CS_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_3
# findPeaks $TAGDIR/set-7_H3K27me3_CS_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_4
# #rep 3
# findPeaks $TAGDIR/WT_H3K27me3_CS_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_5
# findPeaks $TAGDIR/cac-1_H3K27me3_CS_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_5
# findPeaks $TAGDIR/cac-2_H3K27me3_CS_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_5
# findPeaks $TAGDIR/cac-3_H3K27me3_CS_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_4
# findPeaks $TAGDIR/set-7_H3K27me3_CS_rep3 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_4

# #Abcam
# #rep 1
# findPeaks $TAGDIR/WT_H3K27me3_abc_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_1
# findPeaks $TAGDIR/cac-1_H3K27me3_abc_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_1
# findPeaks $TAGDIR/cac-2_H3K27me3_abc_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_1
# findPeaks $TAGDIR/cac-3_H3K27me3_abc_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_1
# findPeaks $TAGDIR/set-7_H3K27me3_abc_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_1
# #rep 2
# findPeaks $TAGDIR/WT_H3K27me3_abc_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_3
# findPeaks $TAGDIR/cac-1_H3K27me3_abc_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_3
# findPeaks $TAGDIR/cac-2_H3K27me3_abc_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_3
# findPeaks $TAGDIR/cac-3_H3K27me3_abc_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_3
# findPeaks $TAGDIR/set-7_H3K27me3_abc_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_3

# #H3K36me3
# #rep 1
# findPeaks $TAGDIR/WT_H3K36me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_2
# findPeaks $TAGDIR/cac-1_H3K36me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_2
# findPeaks $TAGDIR/cac-2_H3K36me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_2
# findPeaks $TAGDIR/cac-3_H3K36me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_2
# findPeaks $TAGDIR/set-7_H3K36me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_2
# #rep 2
# findPeaks $TAGDIR/WT_H3K36me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_5
# findPeaks $TAGDIR/cac-1_H3K36me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_5
# findPeaks $TAGDIR/cac-2_H3K36me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_5
# findPeaks $TAGDIR/cac-3_H3K36me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_4
# findPeaks $TAGDIR/set-7_H3K36me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_4
#
# #H4K20me3
# #rep 1
# findPeaks $TAGDIR/WT_H4K20me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_2
# findPeaks $TAGDIR/cac-1_H4K20me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_2
# findPeaks $TAGDIR/cac-2_H4K20me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_2
# findPeaks $TAGDIR/cac-3_H4K20me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_2
# findPeaks $TAGDIR/set-7_H4K20me3_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_2
# #rep 2
# findPeaks $TAGDIR/WT_H4K20me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_5
# findPeaks $TAGDIR/cac-1_H4K20me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_5
# findPeaks $TAGDIR/cac-2_H4K20me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_5
# findPeaks $TAGDIR/cac-3_H4K20me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_4
# findPeaks $TAGDIR/set-7_H4K20me3_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_4

# #H3K4me2
# #rep 1
# findPeaks $TAGDIR/WT_H3K4me2_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_6
# findPeaks $TAGDIR/cac-1_H3K4me2_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_6
# findPeaks $TAGDIR/cac-2_H3K4me2_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_6
# findPeaks $TAGDIR/cac-3_H3K4me2_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_5
# findPeaks $TAGDIR/set-7_H3K4me2_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_5
# #rep 2
# findPeaks $TAGDIR/WT_H3K4me2_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_input_6
# findPeaks $TAGDIR/cac-1_H3K4me2_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-1_input_6
# findPeaks $TAGDIR/cac-2_H3K4me2_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-2_input_6
# findPeaks $TAGDIR/cac-3_H3K4me2_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/cac-3_input_5
# findPeaks $TAGDIR/set-7_H3K4me2_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/set-7_input_5

#Find Motifs
# findMotifsGenome.pl qasuz12_4hr_peaks.bed /scratch/evt82290/Foxy_Ncrassa_merged.fasta $OUTDIRmotifs/4hr -size given
# findMotifsGenome.pl qasuz12_12hr_peaks.bed /scratch/evt82290/Foxy_Ncrassa_merged.fasta $OUTDIRmotifs/12hr -size given
# findMotifsGenome.pl qasuz12_8hr_no_telo_macs_peaks.bed /home/evt82290/Research/GCA_000182925.2_NC12_genomic_wTetO_at_his3_CLEAN.fasta $OUTDIRmotifs/nucsites -size given -len 8,9,10,11,12,13,14,15

#bedtools

# module load BEDTools
#
# Combining all overlapping peaks
#  bedtools multiinter -header -i ${OUTDIR3}/2024_04_23_WT_peaks.bed \
#                                ${OUTDIR3}/2024_04_23_136_abcam_cac-1_peaks.bed \
#                                ${OUTDIR3}/2024_04_23_136_abcam_cac-2_peaks.bed \
#                                ${OUTDIR3}/2024_04_23_136_abcam_cac-3_peaks.bed \
#                                ${OUTDIR3}/2024_04_23_24hr_peaks.bed > ${OUTDIR3}/merge_peaks.txt
#
#
# bedtools sort -i ${OUTDIR3}/merged_sorted.bed > ${OUTDIR3}/merged_sorted_2.bed
# bedtools merge -i ${OUTDIR3}/merged_sorted_2.bed > ${OUTDIR3}/merged_file.txt
#
# determining which peaks overlap across peak files
# bedtools intersect -wa -wb \
#     -a ${OUTDIR3}/2024_04_23_WT_peaks.bed \
#     -b ${OUTDIR3}/2024_04_23_136_abcam_cac-1_peaks.bed ${OUTDIR3}/2024_04_23_136_abcam_cac-2_peaks.bed ${OUTDIR3}/2024_04_23_136_abcam_cac-3_peaks.bed ${OUTDIR3}/2024_04_23_24hr_peaks.bed \
#     -names cac-1 cac-2 cac-3 24hr \
#     -sorted > ${OUTDIR3}/intersect_peaks.txt

#Run129 bams
# 129-33_ChIP_WT_K27me3_AM_Rep_2_S32.bam
# 129-34_ChIP_cac-1_K27me3_AM_Rep_2_S33.bam
# 129-35_ChIP_cac-2_K27me3_AM_Rep_2_S34.bam
# 129-36_ChIP_cac-3_K27me3_AM_Rep_2_S35.bam
# 129-37_ChIP_set-7_K27me3_AM_Rep_2_S36.bam
# 129-38_ChIP_WT_K27me3_AbC_Rep_1_S37.bam
# 129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38.bam
# 129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39.bam
# 129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40.bam
# 129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41.bam
# 129-43_ChIP_WT_input_S42.bam
# 129-44_ChIP_cac-1_input_S43.bam
# 129-45_ChIP_cac-2_input_S44.bam
# 129-46_ChIP_cac-3_input_S45.bam
# 129-47_ChIP_set-7_input_S46.bam
# 129-48_ChIP-seq_Dim-5_input_S47.bam
# 129-90_ChIP_WT_H3K36me3_Rep1_S71.bam
# 129-91_ChIP_cac-1_H3K36me3_Rep1_S72.bam
# 129-92_ChIP_cac-2_H3K36me3_Rep1_S73.bam
# 129-93_ChIP_cac-3_H3K36me3_Rep1_S74.bam
# 129-94_ChIP_set-7_H3K36me3_Rep1_S75.bam

#Run131 bams
# 131-22_ChIP_WT_H3K27me2me3_Rep1_S20.bam
# 131-23_ChIP_cac-1_H3K27me2me3_Rep1_S21.bam
# 131-24_ChIP_cac-2_H3K27me2me3_Rep1_S22.bam
# 131-25_ChIP_cac-3_H3K27me2me3_Rep1_S23.bam
# 131-26_ChIP_set-7_H3K27me2me3_Rep1_S24.bam
# 131-27_ChIP_WT_H3K36me2_Rep1_S25.bam
# 131-28_ChIP_cac-1_H3K36me2_Rep1_S26.bam
# 131-37_ChIP_WT_input_Rep1_S27.bam
# 131-38_ChIP_cac-1_input_Rep1_S28.bam
# 131-39_ChIP_cac-2_input_Rep1_S29.bam
# 131-40_ChIP_cac-3_input_Rep1_S30.bam
# 131-41_ChIP_set-7_input_Rep1_S31.bam
# 131-82_ChIP_WT_H3K27me3_Rep1_S70.bam
# 131-83_ChIP_cac-1_H3K27me3_Rep1_S71.bam
# 131-84_ChIP_cac-2_H3K27me3_Rep1_S72.bam
# 131-85_ChIP_cac-3_H3K27me3_Rep1_S73.bam
# 131-86_ChIP_set-7_H3K27me3_Rep1_S74.bam
# 131-87_ChIP_WT_H3K27me2me3_Rep3_S75.bam
# 131-88_ChIP_cac-1_H3K27me2me3_Rep3_S76.bam
# 131-89_ChIP_cac-2_H3K27me2me3_Rep3_S77.bam
# 131-90_ChIP_cac-3_H3K27me2me3_Rep3_S78.bam
# 131-91_ChIP_set-7_H3K27me2me3_Rep3_S79.bam

#Run133 bams
# 133-10_ChIP_cac-1_H3K27me2me3_17hr_Rep1_S10.bam
# 133-11_ChIP_cac-2_H3K27me2me3_17hr_Rep1_S11.bam
# 133-14_ChIP_set-7_H3K27me2me3_17hr_Rep1_S12.bam
# 133-15_ChIP_WT_H3K27me2me3_24hr_Rep1_S13.bam
# 133-16_ChIP_cac-1_H3K27me2me3_24hr_Rep1_S14.bam
# 133-17_ChIP_cac-2_H3K27me2me3_24hr_Rep1_S15.bam
# 133-18_ChIP_set-7_H3K27me2me3_24hr_Rep1_S16.bam
# 133-19_ChIP_WT_H4K20me3_Rep1_S17.bam
# 133-1_ChIP_WT_H3K27me2me3_4hr_Rep1_S1.bam
# 133-20_ChIP_cac-1_H4K20me3_Rep1_S18.bam
# 133-21_ChIP_cac-2_H4K20me3_Rep1_S19.bam
# 133-22_ChIP_set-7_H4K20me3_Rep1_S20.bam
# 133-23_ChIP_WT_H3K9me3_Rep2_S21.bam
# 133-24_ChIP_cac-1_H3K9me3_Rep2_S22.bam
# 133-25_ChIP_cac-2_H3K9me3_Rep2_S23.bam
# 133-26_ChIP_set-7_H3K9me3_Rep2_S24.bam
# 133-27_ChIP_WT_H3K4me2_Rep1_S25.bam
# 133-28_ChIP_cac-1_H3K4me2_Rep2_S26.bam
# 133-29_ChIP_cac-2_H3K4me2_Rep1_S27.bam
# 133-2_ChIP_cac-1_H3K27me2me3_4hr_Rep1_S2.bam
# 133-30_ChIP_set-7_H3K4me2_Rep1_S28.bam
# 133-31_ChIP_WT_H3K27me3_4hr_Rep1_S29.bam
# 133-32_ChIP_cac-1_H3K27me3_4hr_Rep1_S30.bam
# 133-33_ChIP_cac-2_H3K27me3_4hr_Rep1_S31.bam
# 133-34_ChIP_set-7_H3K27me3_4hr_Rep1_S32.bam
# 133-35_ChIP_WT_H3K27me3_8hr_Rep1_S33.bam
# 133-36_ChIP_cac-1_H3K27me3_8hr_Rep1_S34.bam
# 133-37_ChIP_cac-2_H3K27me3_8hr_Rep1_S35.bam
# 133-38_ChIP_set-7_H3K27me3_8hr_Rep1_S36.bam
# 133-3_ChIP_cac-2_H3K27me2me3_4hr_Rep1_S3.bam
# 133-40_ChIP_cac-1_H3K27me3_17hr_Rep1_S37.bam
# 133-41_ChIP_cac-2_H3K27me3_17hr_Rep1_S38.bam
# 133-42_ChIP_set-7_H3K27me3_17hr_Rep1_S39.bam
# 133-43_ChIP_WT_H3K27me3_24hr_Rep1_S40.bam
# 133-44_ChIP_cac-1_H3K27me3_24hr_Rep1_S41.bam
# 133-45_ChIP_cac-2_H3K27me3_24hr_Rep1_S42.bam
# 133-46_ChIP_set-7_H3K27me3_24hr_Rep1_S43.bam
# 133-47_ChIP_WT_Acetylated-Lysine_Rep1_S44.bam
# 133-48_ChIP_cac-1_Acetylated-Lysine_Rep1_S45.bam
# 133-49_ChIP_cac-2_Acetylated-Lysine_Rep1_S46.bam
# 133-4_ChIP_set-7_H3K27me2me3_4hr_Rep1_S4.bam
# 133-5_ChIP_WT_H3K27me2me3_8hr_Rep1_S5.bam
# 133-6_ChIP_cac-1_H3K27me2me3_8hr_Rep1_S6.bam
# 133-7_ChIP_cac-2_H3K27me2me3_8hr_Rep1_S7.bam
# 133-8_ChIP_set-7_H3K27me2me3_8hr_Rep1_S8.bam
# 133-95_ChIP_cac-3_Acetylated-Lysine_Rep1_S92.bam
# 133-96_ChIP_set-7_Acetylated-Lysine_Rep1_S93.bam
# 133-9_ChIP_WT_H3K27me2me3_17hr_Rep1_S9.bam

#Run136 bams
# 6147_136-10_ChIP_cac-1-2_H3K27me3_CS_Rep1_S10.bam
# 6147_136-11_ChIP_WT_input_S11.bam
# 6147_136-12_ChIP_cac-1_input_S12.bam
# 6147_136-13_ChIP_cac-2_input_S13.bam
# 6147_136-14_ChIP_cac-3_input_S14.bam
# 6147_136-1_ChIP_WT_H3K27me3_abcam_Rep2_S1.bam
# 6147_136-21_ChIP_cac-1-2_input_S21.bam
# 6147_136-22_ChIP_WT_H3K27me3_CS_0hr_Rep1_S22.bam
# 6147_136-23_ChIP_qa-suz12_H3K27me3_CS_0hr_Rep1_S23.bam
# 6147_136-24_ChIP_qa-suz12_cac-1_H3K27me3_CS_0hr_Rep1_S24.bam
# 6147_136-25_ChIP_WT_H3K27me3_CS_8hr_Rep1_S25.bam
# 6147_136-26_ChIP_qa-suz12_H3K27me3_CS_8hr_Rep1_S26.bam
# 6147_136-27_ChIP_qa-suz12_cac-1_H3K27me3_CS_8hr_Rep1_S27.bam
# 6147_136-28_ChIP_WT_H3K27me3_CS_24hr_Rep1_S28.bam
# 6147_136-29_ChIP_qa-suz12_H3K27me3_CS_24hr_Rep1_S29.bam
# 6147_136-2_ChIP_cac-1_H3K27me3_abcam_Rep2_S2.bam
# 6147_136-39_ChIP_qa-suz12_cac-1_H3K27me3_CS_24hr_Rep1_S39.bam
# 6147_136-3_ChIP_cac-2_H3K27me3_abcam_Rep2_S3.bam
# 6147_136-40_ChIP_WT_H3K27me3_abcam_6hr_Rep1_S40.bam
# 6147_136-41_ChIP_qa-suz12_H3K27me3_abcam_6hr_Rep1_S41.bam
# 6147_136-42_ChIP_qa-suz12_H3K27me3_abcam_6hr_Rep1_S42.bam
# 6147_136-43_ChIP_qa-suz12_cac-1_H3K27me3_abcam_6hr_Rep1_S43.bam
# 6147_136-44_ChIP_qa-suz12_cac-1_H3K27me3_abcam_6hr_Rep1_S44.bam
# 6147_136-45_ChIP_set-7_H3K27me3_abcam_6hr_Rep1_S45.bam
# 6147_136-4_ChIP_cac-3_H3K27me3_abcam_Rep2_S4.bam
# 6147_136-58_ChIP_WT_H3K27me3_abcam_Rep3_S58.bam
# 6147_136-5_ChIP_cac-1-2_H3K27me3_abcam_Rep2_S5.bam
# 6147_136-6_ChIP_WT_H3K27me3_CS_Rep1_S6.bam
# 6147_136-71_ChIP_cac-1_H3K27me3_abcam_Rep3_S70.bam
# 6147_136-72_ChIP_cac-2_H3K27me3_abcam_Rep3_S71.bam
# 6147_136-75_ChIP_cac-3_H3K27me3_abcam_Rep3_S74.bam
# 6147_136-76_ChIP_cac-1-2_H3K27me3_abcam_Rep3_S75.bam
# 6147_136-77_ChIP_set-7_H3K27me3_abcam_Rep3_S76.bam
# 6147_136-78_ChIP_WT_H3K27me3_CS_Rep2_S77.bam
# 6147_136-79_ChIP_cac-1_H3K27me3_CS_Rep2_S78.bam
# 6147_136-7_ChIP_cac-1_H3K27me3_CS_Rep1_S7.bam
# 6147_136-80_ChIP_cac-2_H3K27me3_CS_Rep2_S79.bam
# 6147_136-81_ChIP_cac-3_H3K27me3_CS_Rep2_S80.bam
# 6147_136-82_ChIP_cac-1-2_H3K27me3_CS_Rep2_S81.bam
# 6147_136-83_ChIP_set-7_H3K27me3_CS_Rep2_S82.bam
# 6147_136-84_ChIP_WT_input_S83.bam
# 6147_136-85_ChIP_cac-1_input_S84.bam
# 6147_136-89_ChIP_cac-2_input_S88.bam
# 6147_136-8_ChIP_cac-2_H3K27me3_CS_Rep1_S8.bam
# 6147_136-91_ChIP_cac-1-2_input_S90.bam
# 6147_136-92_ChIP_set-7_input_S91.bam
# 6147_136-9_ChIP_cac-3_H3K27me3_CS_Rep1_S9.bam

#Run137 bams
# 137-60_ChIP_WT_HA_Rep1_S57.bam
# 137-61_ChIP_ETX51-cac-1-HA_HA_Rep1_S58.bam
# 137-62_ChIP_ETX52-cac-1-HA_HA_Rep1_S59.bam
# 137-63_ChIP_WT_input__S60.bam
# 137-64_ChIP_ETX51-cac-1-HA_input__S61.bam
# 137-65_ChIP_ETX52-cac-1-HA_input__S62.bam
# 137-66_ChIP_WT_H3K27me3_Rep1_S63.bam
# 137-67_ChIP_qa-suz12_H3K27me3_Rep1_S64.bam
# 137-68_ChIP_WT_H3K27me3_Rep1_S65.bam
# 137-69_ChIP_qa-suz12_H3K27me3_Rep1_S66.bam
# 137-70_ChIP_WT_H3K27me3_Rep1_S67.bam
# 137-71_ChIP_qa-suz12_H3K27me3_Rep1_S68.bam
# 137-72_ChIP_WT_H3K27me3_Rep1_S69.bam
# 137-73_ChIP_qa-suz12_H3K27me3_Rep1_S70.bam
# 137-74_ChIP_WT_H3K27me3_Rep1_S71.bam
# 137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72.bam

#Run138 bams
# 138-57_ChIP_WT_H3K27me3_Rep3_6252_S56.bam
# 138-58_ChIP_cac-1_H3K27me3_Rep3_6252_S57.bam
# 138-59_ChIP_cac-2_H3K27me3_Rep3_6252_S58.bam
# 138-60_ChIP_cac-3_H3K27me3_Rep3_6252_S59.bam
# 138-61_ChIP_set-7_H3K27me3_Rep3_6252_S60.bam
# 138-62_ChIP_WT_H3K36me3_Rep2_6252_S61.bam
# 138-63_ChIP_cac-1_H3K36me3_Rep2_6252_S62.bam
# 138-64_ChIP_cac-2_H3K36me3_Rep2_6252_S63.bam
# 138-65_ChIP_cac-3_H3K36me3_Rep2_6252_S64.bam
# 138-66_ChIP_set-7_H3K36me3_Rep2_6252_S65.bam
# 138-67_ChIP_WT_H4K20me3_Rep2_6252_S66.bam
# 138-68_ChIP_cac-1_H4K20me3_Rep2_6252_S67.bam
# 138-69_ChIP_cac-2_H4K20me3_Rep2_6252_S68.bam
# 138-70_ChIP_cac-3_H4K20me3_Rep2_6252_S69.bam
# 138-71_ChIP_set-7_H4K20me3_Rep2_6252_S70.bam
# 138-72_ChIP_WT_input__6252_S71.bam
# 138-73_ChIP_cac-1_input__6252_S72.bam
# 138-74_ChIP_cac-2_input__6252_S73.bam
# 138-75_ChIP_cac-3_input__6252_S74.bam
# 138-76_ChIP_set-7_input__6252_S75.bam
# 138-77_ChIP_WT_0hr_H3K27me3_Rep1_6252_S76.bam
# 138-78_ChIP_tetO-cac-1_0hr_H3K27me3_Rep1_6252_S77.bam
# 138-79_ChIP_tetO-cac-2_0hr_H3K27me3_Rep1_6252_S78.bam
# 138-80_ChIP_tetO-cac-3_0hr_H3K27me3_Rep1_6252_S79.bam
# 138-81_ChIP_tetO_0hr_H3K27me3_Rep1_6252_S80.bam
# 138-82_ChIP_WT_6hr_H3K27me3_Rep1_6252_S81.bam
# 138-83_ChIP_tetO-cac-1_6hr_H3K27me3_Rep1_6252_S82.bam
# 138-84_ChIP_tetO-cac-2_6hr_H3K27me3_Rep1_6252_S83.bam
# 138-85_ChIP_tetO-cac-3_6hr_H3K27me3_Rep1_6252_S84.bam
# 138-86_ChIP_tetO_6hr_H3K27me3_Rep1_6252_S85.bam
# 138-87_ChIP_WT_12hr_H3K27me3_Rep1_6252_S86.bam
# 138-88_ChIP_tetO-cac-1_12hr_H3K27me3_Rep1_6252_S87.bam
# 138-89_ChIP_tetO-cac-2_12hr_H3K27me3_Rep1_6252_S88.bam
# 138-90_ChIP_tetO-cac-3_12hr_H3K27me3_Rep1_6252_S89.bam
# 138-91_ChIP_tetO_12hr_H3K27me3_Rep1_6252_S90.bam
# 138-92_ChIP_WT_24hr_H3K27me3_Rep1_6252_S91.bam
# 138-93_ChIP_tetO-cac-1_24hr_H3K27me3_Rep1_6252_S92.bam
# 138-94_ChIP_tetO-cac-2_24hr_H3K27me3_Rep1_6252_S93.bam
# 138-95_ChIP_tetO-cac-3_24hr_H3K27me3_Rep1_6252_S94.bam
# 138-96_ChIP_tetO_24hr_H3K27me3_Rep1_6252_S95.bam

#Run139 bam
# 139-29_ChIP_WT_Input__S29.bam
# 139-30_ChIP_qa-suz12_Input__S30.bam
# 139-31_ChIP_WT_Input__S31.bam
# 139-32_ChIP_qa-suz12_Input__S32.bam
# 139-33_ChIP_WT_Input__S33.bam
# 139-34_ChIP_qa-suz12_Input__S34.bam
# 139-35_ChIP_WT_Input__S35.bam
# 139-36_ChIP_qa-suz12_Input__S36.bam
# 139-37_ChIP_WT_Input__S37.bam
# 139-38_ChIP_qa-suz12_Input__S38.bam
# 139-49_ChIP_qa-suz12_H3K27me3_Rep3_S40.bam
# 139-50_ChIP_qa-suz12-cac-1_H3K27me3_Rep2_S41.bam
# 139-51_ChIP_qa-suz12_H3K27me3_Rep1_S42.bam
# 139-52_ChIP_qa-suz12-cac-1_H3K27me3_Rep1_S43.bam
# 139-53_ChIP_qa-suz12_H3K27me3_Rep1_S44.bam
# 139-54_ChIP_qa-suz12-cac-1_H3K27me3_Rep1_S45.bam
# 139-55_ChIP_qa-suz12_H3K27me3_Rep3_S46.bam
# 139-56_ChIP_qa-suz12-cac-1_H3K27me3_Rep2_S47.bam
# 139-57_ChIP_WT_H3K27me3_Rep1_S48.bam
# 139-58_ChIP_qa-suz12_H3K27me3_Rep1_S49.bam
# 139-59_ChIP_qa-suz12-cac-1_H3K27me3_Rep1_S50.bam
# 139-60_ChIP_WT_H3K27me2me3_Rep1_S51.bam
# 139-61_ChIP_qa-suz12_H3K27me2me3_Rep1_S52.bam
# 139-62_ChIP_qa-suz12-cac-1_H3K27me2me3_Rep1_S53.bam
# 139-63_ChIP_qa-suz12_H3K27me2me3_Rep1_S54.bam
# 139-64_ChIP_qa-suz12-cac-1_H3K27me2me3_Rep1_S55.bam
# 139-65_ChIP_qa-suz12_H3K27me2me3_Rep1_S56.bam
# 139-66_ChIP_qa-suz12-cac-1_H3K27me2me3_Rep1_S57.bam
# 139-67_ChIP_qa-suz12_H3K27me2me3_Rep1_S58.bam
# 139-68_ChIP_qa-suz12-cac-1_H3K27me2me3_Rep1_S59.bam
# 139-69_ChIP_WT_H3K27me2me3_Rep1_S60.bam
# 139-70_ChIP_qa-suz12_H3K27me2me3_Rep1_S61.bam

#Run141 bams
# 141-58_ChIP_WT_H3K23me1_Rep1_S52.bam
# 141-59_ChIP_cac-1_H3K23me1_Rep1_S53.bam
# 141-60_ChIP_cac-2_H3K23me1_Rep1_S54.bam
# 141-61_ChIP_cac-3_H3K23me1_Rep1_S55.bam
# 141-62_ChIP_set-7_H3K23me1_Rep1_S56.bam
# 141-63_ChIP_WT_teto_P0_H3K27me3_Rep1_S57.bam
# 141-64_ChIP_WT_teto_P1_H3K27me3_Rep1_S58.bam
# 141-65_ChIP_WT_teto_P2_H3K27me3_Rep1_S59.bam
# 141-66_ChIP_WT_teto_P3_H3K27me3_Rep1_S60.bam
# 141-67_ChIP_cac-1_teto_P0_H3K27me3_Rep1_S61.bam
# 141-68_ChIP_cac-1_teto_P1_H3K27me3_Rep1_S62.bam
# 141-69_ChIP_cac-1_teto_P2_H3K27me3_Rep1_S63.bam
# 141-70_ChIP_cac-1_teto_P3_H3K27me3_Rep1_S64.bam
# 141-88_ChIP_WT_H3K27me2_Rep1_S82.bam
# 141-89_ChIP_cac-1_H3K27me2_Rep1_S83.bam
# 141-90_ChIP_cac-2_H3K27me2_Rep1_S84.bam
# 141-91_ChIP_cac-3_H3K27me2_Rep1_S85.bam
# 141-92_ChIP_set-7_H3K27me2_Rep1_S86.bam
# 141-93_ChIP_WT_H4K12ac_Rep1_S87.bam
# 141-94_ChIP_cac-1_H4K12ac_Rep1_S88.bam
# 141-95_ChIP_cac-2_H4K12ac_Rep1_S89.bam
# 141-96_ChIP_set-7_H4K12ac_Rep1_S90.bam
# 141-N11_ATAC_WT__Rep1_S100.bam
# 141-N12_ATAC_cac-1__Rep1_S101.bam
# 141-N13_ATAC_cac-2__Rep1_S102.bam
# 141-N14_ATAC_cac-3__Rep1_S103.bam
# 141-N15_ATAC_set-7__Rep1_S104.bam


# Run144 bams
# 144-10_ChIP_tetO_p0_H3K27me3_Rep3_S10.bam
# 144-11_ChIP_tetO_p1_H3K27me3_Rep3_S11.bam
# 144-122_ChIP_tetO_p3_input__S122.bam
# 144-124_ChIP_tetO_p1_input__S124.bam
# 144-125_ChIP_tetO_p2_input__S125.bam
# 144-126_ChIP_tetO_p3_input__S126.bam
# 144-127_ChIP_tetO_p0_tet_input__S127.bam
# 144-12_ChIP_tetO_p2_H3K27me3_Rep3_S12.bam
# 144-132_ChIP_qa-suz12_Input__S132.bam
# 144-13_ChIP_tetO_p3_H3K27me3_Rep3_S13.bam
# 144-141_ChIP_qa-suz12_24hr_H3K27me2_Rep1_S141.bam
# 144-14_ChIP_tetO_p0_H3K27me3_Rep4_S14.bam
# 144-15_ChIP_tetO_p1_H3K27me3_Rep4_S15.bam
# 144-16_ChIP_tetO_p2_H3K27me3_Rep4_S16.bam
# 144-17_ChIP_tetO_p3_H3K27me3_Rep4_S17.bam
# 144-18_ChIP_tetO_p0_H3K36me3_Rep1_S18.bam
# 144-19_ChIP_tetO_p1_H3K36me3_Rep1_S19.bam
# 144-1_ChIP_tetO_p0_H3K27me3_Rep2_S1.bam
# 144-20_ChIP_tetO_p2_H3K36me3_Rep1_S20.bam
# 144-21_ChIP_tetO_p3_H3K36me3_Rep1_S21.bam
# 144-22_ChIP_tetO_p0_H3K36me3_Rep2_S22.bam
# 144-23_ChIP_tetO_p1_H3K36me3_Rep2_S23.bam
# 144-24_ChIP_tetO_p2_H3K36me3_Rep2_S24.bam
# 144-25_ChIP_tetO_p3_H3K36me3_Rep2_S25.bam
# 144-26_ChIP_tetO_p0_tet_H3K27me3_Rep2_S26.bam
# 144-27_ChIP_tetO_p0_tet_H3K36me3_Rep2_S27.bam
# 144-28_ChIP_WT_0hr_H3K27me3_Rep3_S28.bam
# 144-29_ChIP_WT_24hr_H3K27me3_Rep2_S29.bam
# 144-2_ChIP_tetO_p1_H3K27me3_Rep2_S2.bam
# 144-30_ChIP_qa-suz12_0hr_H3K27me3_Rep3_S30.bam
# 144-31_ChIP_qa-suz12_4hr_H3K27me3_Rep2_S31.bam
# 144-32_ChIP_qa-suz12_8hr_H3K27me3_Rep2_S32.bam
# 144-33_ChIP_qa-suz12_12hr_H3K27me3_Rep2_S33.bam
# 144-34_ChIP_qa-suz12_24hr_H3K27me3_Rep3_S34.bam
# 144-35_ChIP_qa-suz12_48hr_H3K27me3_Rep2_S35.bam
# 144-37_ChIP_suz12_24hr_H3K27me3_Rep1_S37.bam
# 144-38_ChIP_qa-suz12_cac-1_0hr_H3K27me3_Rep2_S38.bam
# 144-39_ChIP_qa-suz12_cac-1_4hr_H3K27me3_Rep1_S39.bam
# 144-3_ChIP_tetO_p2_H3K27me3_Rep2_S3.bam
# 144-40_ChIP_qa-suz12_cac-1_8hr_H3K27me3_Rep1_S40.bam
# 144-41_ChIP_qa-suz12_cac-1_12hr_H3K27me3_Rep2_S41.bam
# 144-42_ChIP_qa-suz12_cac-1_24hr_H3K27me3_Rep2_S42.bam
# 144-43_ChIP_qa-suz12_cac-1_48hr_H3K27me3_Rep2_S43.bam
# 144-44_ChIP_WT_0hr_H3K27me3_Rep4_S44.bam
# 144-45_ChIP_WT_24hr_H3K27me3_Rep3_S45.bam
# 144-46_ChIP_suz12_0hr_H3K27me3_Rep2_S46.bam
# 144-47_ChIP_suz12_24hr_H3K27me3_Rep2_S47.bam
# 144-48_ChIP_qa-suz12_0hr_H3K27me3_Rep4_S48.bam
# 144-49_ChIP_qa-suz12_4hr_H3K27me3_Rep3_S49.bam
# 144-4_ChIP_tetO_p3_H3K27me3_Rep2_S4.bam
# 144-50_ChIP_qa-suz12_8hr_H3K27me3_Rep3_S50.bam
# 144-51_ChIP_qa-suz12_12hr_H3K27me3_Rep3_S51.bam
# 144-52_ChIP_qa-suz12_24hr_H3K27me3_Rep4_S52.bam
# 144-53_ChIP_WT_0hr_H3K27me2_Rep1_S53.bam
# 144-54_ChIP_WT_24hr_H3K27me2_Rep1_S54.bam
# 144-55_ChIP_suz12_0hr_H3K27me2_Rep1_S55.bam
# 144-56_ChIP_suz12_24hr_H3K27me2_Rep1_S56.bam
# 144-57_ChIP_qa-suz12_0hr_H3K27me2_Rep1_S57.bam
# 144-58_ChIP_qa-suz12_4hr_H3K27me2_Rep1_S58.bam
# 144-59_ChIP_qa-suz12_8hr_H3K27me2_Rep1_S59.bam
# 144-5_ChIP_tetO_p0_tet_H3K27me3_Rep2_S5.bam
# 144-60_ChIP_qa-suz12_12hr_H3K27me2_Rep1_S60.bam
# 144-62_ChIP_WT_Input__S62.bam
# 144-69_ChIP_suz12_Input__S69.bam
# 144-6_ChIP_tetO_cac-1_p0_H3K27me3_Rep2_S6.bam
# 144-71_ChIP_qa-suz12_Input__S71.bam
# 144-73_ChIP_qa-suz12_Input__S73.bam
# 144-79_ChIP_qa-suz12_Input__S79.bam
# 144-7_ChIP_tetO_cac-1_p1_H3K27me3_Rep2_S7.bam
# 144-8_ChIP_tetO_cac-1_p2_H3K27me3_Rep2_S8.bam
# 144-9_ChIP_tetO_cac-1_p3_H3K27me3_Rep2_S9.bam

#Run145 bams
# 145-10_ChIP_csr1_LGVI_2_nuc_H3K27me3_Rep1_S10.bam
# 145-11_ChIP_WT_H3K27me3_nuc_S11.bam
# 145-15_ChIP_tetO_0hr_H3K27me3_Rep1_S12.bam
# 145-16_ChIP_tetO_6hr_H3K27me3_Rep1_S13.bam
# 145-17_ChIP_tetO_12hr_H3K27me3_Rep1_S14.bam
# 145-18_ChIP_tetO_18hr_H3K27me3_Rep1_S15.bam
# 145-19_ChIP_tetO_24hr_H3K27me3_Rep1_S16.bam
# 145-1_ChIP_csr1_LGIII_2_nuc_H3K27me3_Rep1_S1.bam
# 145-20_ChIP_tetO_30hr_H3K27me3_Rep1_S17.bam
# 145-21_ChIP_tetO_36hr_H3K27me3_Rep1_S18.bam
# 145-22_ChIP_tetO_cac-1_0hr_H3K27me3_Rep1_S19.bam
# 145-23_ChIP_tetO_cac-1_6hr_H3K27me3_Rep1_S20.bam
# 145-24_ChIP_tetO_cac-1_12hr_H3K27me3_Rep1_S21.bam
# 145-25_ChIP_tetO_cac-1_18hr_H3K27me3_Rep1_S22.bam
# 145-26_ChIP_tetO_cac-1_24hr_H3K27me3_Rep1_S23.bam
# 145-27_ChIP_tetO_cac-1_30hr_H3K27me3_Rep1_S24.bam
# 145-28_ChIP_tetO_cac-1_36hr_H3K27me3_Rep1_S25.bam
# 145-2_ChIP_csr1_LGIII_2_nuc_H3K27me3_Rep1_S2.bam
# 145-3_ChIP_csr1_LGV_1_nuc_H3K27me3_Rep1_S3.bam
# 145-4_ChIP_csr1_LGVI_1_nuc_H3K27me3_Rep1_S4.bam
# 145-5_ChIP_csr1_LGIII_1_nuc_H3K27me3_Rep1_S5.bam
# 145-60_ChIP_WT_0hr_input__S57.bam
# 145-6_ChIP_csr1_LGIII_2_nuc_H3K27me3_Rep1_S6.bam
# 145-7_ChIP_csr1_LGV_1_nuc_H3K27me3_Rep1_S7.bam
# 145-8_ChIP_csr1_LGVI_1_nuc_H3K27me3_Rep1_S8.bam
# 145-9_ChIP_csr1_LGVI_1_nuc_H3K27me3_Rep1_S9.bam

#Run146 bams
# 146-106_ChIP_epr-1GFP_GFP-trap_Rep1_S125.bam
# 146-10_ChIP_qa-suz12_72hr_H3K27me2_Rep1_S10.bam
# 146-117_ChIP_qa-suz12_96hr_H3K27me2_Rep1_S136.bam
# 146-118_ChIP_epr-1_GFP_H3K27me3_Rep1_S137.bam
# 146-11_ChIP_qa-suz12_96hr_H3K27me2_Rep1_S11.bam
# 146-123_ChIP_WT_H3K27me2_Rep1_S142.bam
# 146-124_ChIP_suz12_0hr_H3K27me3_Rep1_S143.bam
# 146-125_ChIP_qa-suz12_0hr_H3K27me3_Rep1_S144.bam
# 146-126_ChIP_qa-suz12_4hr_H3K27me3_Rep1_S145.bam
# 146-127_ChIP_qa-suz12_8hr_H3K27me3_Rep1_S146.bam
# 146-128_ChIP_qa-suz12_12hr_H3K27me3_Rep1_S147.bam
# 146-129_ChIP_qa-suz12_24hr_H3K27me3_Rep1_S148.bam
# 146-12_ChIP_WT_H3K27me2_Rep1_S12.bam
# 146-130_ChIP_suz12_0hr_H3K36me3_Rep1_S149.bam
# 146-131_ChIP_qa-suz12_0hr_H3K36me3_Rep1_S150.bam
# 146-132_ChIP_qa-suz12_4hr_H3K36me3_Rep1_S151.bam
# 146-133_ChIP_qa-suz12_8hr_H3K36me3_Rep1_S152.bam
# 146-134_ChIP_qa-suz12_12hr_H3K36me3_Rep1_S153.bam
# 146-135_ChIP_qa-suz12_24hr_H3K36me3_Rep1_S154.bam
# 146-136_ChIP_cac-2-epr-1-GFP_GFP-trap_Rep1_S155.bam
# 146-137_ChIP_cac-2-epr-1-GFP_GFP-trap_Rep1_S156.bam
# 146-13_ChIP_qa-suz12_0hr_input__S13.bam
# 146-14_ChIP_qa-suz12_24hr_input__S14.bam
# 146-15_ChIP_qa-suz12_48hr_input__S15.bam
# 146-16_ChIP_qa-suz12_72hr_input__S16.bam
# 146-17_ChIP_qa-suz12_96hr_input__S17.bam
# 146-18_ChIP_WT_input__S18.bam
# 146-19_ChIP_WT_H3K4me2_Rep1_S19.bam
# 146-1_ChIP_qa-suz12_0hr_H3K27me3_Rep4_S1.bam
# 146-20_ChIP_cac-1_H3K4me2_Rep1_S20.bam
# 146-21_ChIP_cac-2_H3K4me2_Rep1_S21.bam
# 146-22_ChIP_cac-3_H3K4me2_Rep1_S22.bam
# 146-23_ChIP_set-7_H3K4me2_Rep1_S23.bam
# 146-24_ChIP_WT_H3K27me2_Rep2_S24.bam
# 146-25_ChIP_cac-1_H3K27me2_Rep2_S25.bam
# 146-26_ChIP_cac-2_H3K27me2_Rep2_S26.bam
# 146-27_ChIP_cac-3_H3K27me2_Rep2_S27.bam
# 146-28_ChIP_set-7_H3K27me2_Rep2_S28.bam
# 146-29_ChIP_WT_H3K4me2_Rep2_S29.bam
# 146-2_ChIP_qa-suz12_24hr_H3K27me3_Rep4_S2.bam
# 146-30_ChIP_cac-1_H3K4me2_Rep2_S30.bam
# 146-31_ChIP_cac-2_H3K4me2_Rep2_S31.bam
# 146-32_ChIP_cac-3_H3K4me2_Rep2_S32.bam
# 146-33_ChIP_set-7_H3K4me2_Rep2_S33.bam
# 146-34_ChIP_WT_input__S34.bam
# 146-35_ChIP_cac-1_input__S35.bam
# 146-36_ChIP_cac-2_input__S36.bam
# 146-37_ChIP_cac-3_input__S37.bam
# 146-38_ChIP_set-7_input__S38.bam
# 146-39_ChIP_qa-suz12_0hr_H3K27me3_Rep1_S39.bam
# 146-3_ChIP_qa-suz12_48hr_H3K27me3_Rep1_S3.bam
# 146-40_ChIP_qa-suz12_24hr_H3K27me3_Rep1_S40.bam
# 146-41_ChIP_qa-suz12_48hr_H3K27me3_Rep1_S41.bam
# 146-42_ChIP_qa-suz12_72hr_H3K27me3_Rep1_S42.bam
# 146-43_ChIP_qa-suz12_96hr_H3K27me3_Rep1_S43.bam
# 146-44_ChIP_WT_H3K27me3_Rep1_S44.bam
# 146-45_ChIP_qa-suz12_0hr_H3K27me2_Rep1_S45.bam
# 146-46_ChIP_qa-suz12_24hr_H3K27me2_Rep1_S46.bam
# 146-47_ChIP_qa-suz12_48hr_H3K27me2_Rep1_S47.bam
# 146-48_ChIP_qa-suz12_72hr_H3K27me2_Rep1_S48.bam
# 146-4_ChIP_qa-suz12_72hr_H3K27me3_Rep1_S4.bam
# 146-5_ChIP_qa-suz12_96hr_H3K27me3_Rep1_S5.bam
# 146-6_ChIP_WT_H3K27me3_Rep1_S6.bam
# 146-7_ChIP_qa-suz12_0hr_H3K27me2_Rep1_S7.bam
# 146-8_ChIP_qa-suz12_24hr_H3K27me2_Rep1_S8.bam
# 146-9_ChIP_qa-suz12_48hr_H3K27me2_Rep1_S9.bam
# 146-N1_ATAC_WT__Rep2_S97.bam
# 146-N3_ATAC_cac-1__Rep2_S98.bam
# 146-N4_ATAC_cac-2__Rep2_S99.bam
# 146-N5_ATAC_set-7__Rep2_S100.bam
# 146-N6_ATAC_cac-3__Rep2_S101.bam
