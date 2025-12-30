#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --time=8:00:00
#SBATCH --output=../RNAseqMap/logs/featureCounts.%j.out
#SBATCH --error=../RNAseqMap/logs/featureCounts.%j.err

#This code will only be used for building a count file across multiple Samples/Strains

cd $SLURM_SUBMIT_DIR

THREADS=2

#Make & Assign Directories
BAMDIR="/scratch/evt82290/RNAseq/qa-suz12/bamFiles"
BAMDIR2="/scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles"
B149DIR="/scratch/evt82290/MappingOutputs/Run149/RNA/bamFiles"
B150DIR="/scratch/evt82290/MappingOutputs/Run150/bamFiles"
OUTDIR="/scratch/evt82290/RNAseq/qa-suz12/counts"

# #if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
###
# /scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles/SRR7970630/SRR7970630_Aligned.sortedByCoord.out.bam
# /scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles/SRR7970630/SRR7970630_Aligned.sortedByCoord.out.bam
#Run featureCounts
module load Subread/

featureCounts -T $THREADS \
  -p \
  -t exon \
  -g gene_name \
  -s 2 \
  --primary \
  -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  -o $OUTDIR/readcounts_qa_paper_new_TEST.txt \
$BAMDIR/SRR9027658/SRR9027658_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR9027689/SRR9027689_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR9027759/SRR9027759_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-85_RNA_WT__Rep2_S85/149-85_RNA_WT__Rep2_S85_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-86_RNA_WT__Rep3_S86/149-86_RNA_WT__Rep3_S86_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-87_RNA_WT_24hrqa__Rep1_S87/149-87_RNA_WT_24hrqa__Rep1_S87_Aligned.sortedByCoord.out.bam \
    $B149DIR/149-89_RNA_WT_24hrqa__Rep3_S89/149-89_RNA_WT_24hrqa__Rep3_S89_Aligned.sortedByCoord.out.bam \
    $B150DIR/150-108_RNA_WT_24hr__Rep4_S123/150-108_RNA_WT_24hr__Rep4_S123_Aligned.sortedByCoord.out.bam \
    $B150DIR/150-109_RNA_WT_24hr__Rep5_S124/150-109_RNA_WT_24hr__Rep5_S124_Aligned.sortedByCoord.out.bam \
      $B149DIR/149-91_RNA_qa-suz12__Rep2_S91/149-91_RNA_qa-suz12__Rep2_S91_Aligned.sortedByCoord.out.bam \
      $B150DIR/150-112_RNA_qa_0hr__Rep5_S127/150-112_RNA_qa_0hr__Rep5_S127_Aligned.sortedByCoord.out.bam \
      $B149DIR/149-90_RNA_qa-suz12__Rep1_S90/149-90_RNA_qa-suz12__Rep1_S90_Aligned.sortedByCoord.out.bam \
          $B149DIR/149-95_RNA_qa-suz12_24hrqa__Rep3_S95/149-95_RNA_qa-suz12_24hrqa__Rep3_S95_Aligned.sortedByCoord.out.bam \
          $B150DIR/150-114_RNA_qa_24hr__Rep5_S129/150-114_RNA_qa_24hr__Rep5_S129_Aligned.sortedByCoord.out.bam \
          $B150DIR/150-115_RNA_qa_24hr__Rep6_S130/150-115_RNA_qa_24hr__Rep6_S130_Aligned.sortedByCoord.out.bam

###NOTE: 149-87 is actually 0 hr 149-84 is actually 24 hr

#24hr qa
# $B149DIR/149-93_RNA_qa-suz12_24hrqa__Rep1_S93/149-93_RNA_qa-suz12_24hrqa__Rep1_S93_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-94_RNA_qa-suz12_24hrqa__Rep2_S94/149-94_RNA_qa-suz12_24hrqa__Rep2_S94_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-95_RNA_qa-suz12_24hrqa__Rep3_S95/149-95_RNA_qa-suz12_24hrqa__Rep3_S95_Aligned.sortedByCoord.out.bam \

#0hr qa
# $B149DIR/149-90_RNA_qa-suz12__Rep1_S90/149-90_RNA_qa-suz12__Rep1_S90_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-92_RNA_qa-suz12__Rep3_S92/149-92_RNA_qa-suz12__Rep3_S92_Aligned.sortedByCoord.out.bam \

# $B150DIR/150-111_RNA_qa_0hr__Rep4_S126/150-111_RNA_qa_0hr__Rep4_S126_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-92_RNA_qa-suz12__Rep3_S92/149-92_RNA_qa-suz12__Rep3_S92_Aligned.sortedByCoord.out.bam \

#24hr WT
# $B149DIR/149-84_RNA_WT__Rep1_S84/149-84_RNA_WT__Rep1_S84_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-88_RNA_WT_24hrqa__Rep2_S88/149-88_RNA_WT_24hrqa__Rep2_S88_Aligned.sortedByCoord.out.bam \
# $B150DIR/150-110_RNA_WT_24hr__Rep6_S125/150-110_RNA_WT_24hr__Rep6_S125_Aligned.sortedByCoord.out.bam \

#0hr WT
# $B149DIR/149-86_RNA_WT__Rep3_S86/149-86_RNA_WT__Rep3_S86_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-87_RNA_WT_24hrqa__Rep1_S87/149-87_RNA_WT_24hrqa__Rep1_S87_Aligned.sortedByCoord.out.bam \

# $B150DIR/150-106_RNA_WT_0hr__Rep4_S121/150-106_RNA_WT_0hr__Rep4_S121_Aligned.sortedByCoord.out.bam \
# $B150DIR/150-107_RNA_WT_0hr__Rep5_S122/150-107_RNA_WT_0hr__Rep5_S122_Aligned.sortedByCoord.out.bam \

###WT FROM SCREEN
# $BAMDIR2/SRR8269825/SRR8269825_Aligned.sortedByCoord.out.bam \
# $BAMDIR2/SRR8269775/SRR8269775_Aligned.sortedByCoord.out.bam \
# $BAMDIR2/SRR8269782/SRR8269782_Aligned.sortedByCoord.out.bam \
# $BAMDIR2/SRR8269810/SRR8269810_Aligned.sortedByCoord.out.bam \
