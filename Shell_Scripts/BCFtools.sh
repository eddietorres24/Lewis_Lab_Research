#!/bin/bash
#SBATCH --job-name=BCFtools
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=08:00:00
#SBATCH --output=../BCFtools/logs/%x.out
#SBATCH --error=../BCFtools/logs/%x.err

#change directory & load BCFtools
cd $SLURM_SUBMIT_DIR
ml BCFtools

#set working & output directories
OUTDIR="/scratch/evt82290/BCFtools/cac_analysis"
MAPDIR="/scratch/evt82290/MappingOutputs"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#Run bcftools
bcftools view --call -s -O z -o WT_v_cac

/scratch/evt82290/MappingOutputs/Run138/bamFiles

#
bcftools mpileup -O u -f /home/evt82290/Research/GCA_000182925.2_NC12_genomic_wTetO_at_his3_CLEAN.fasta $MAPDIR/Run129/bamFiles/129-43_ChIP_WT_input_S42.bam $MAPDIR/Run138/bamFiles/138-72_ChIP_WT_input__6252_S71.bam  $MAPDIR/Run136/bamFiles/6147_136-11_ChIP_WT_input_S11.bam $MAPDIR/Run136/bamFiles/6147_136-84_ChIP_WT_input_S83.bam $MAPDIR/Run131/bamFiles/131-37_ChIP_WT_input_Rep1_S27.bam $MAPDIR/Run138/bamFiles/138-73_ChIP_cac-1_input__6252_S72.bam $MAPDIR/Run131/bamFiles/131-38_ChIP_cac-1_input_Rep1_S28.bam $MAPDIR/Run129/bamFiles/129-44_ChIP_cac-1_input_S43.bam $MAPDIR/Run136/bamFiles/6147_136-12_ChIP_cac-1_input_S12.bam $MAPDIR/Run136/bamFiles/6147_136-85_ChIP_cac-1_input_S84.bam | bcftools call -m -v -o $OUTDIR/raw_variants.vcf
bcftools filter -i 'GT[0-2]="RR" && GT[3-5]="AA"' $OUTDIR/raw_variants.vcf > $OUTDIR/filtered_variants.vcf
