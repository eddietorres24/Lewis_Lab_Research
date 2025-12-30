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

# Build a count file across multiple Samples/Strains (PE first, then SE, then merge)

cd "$SLURM_SUBMIT_DIR"

THREADS=2

# Make & Assign Directories
BAMDIR="/scratch/evt82290/MappingOutputs/RNAseq/bamFiles"
B149DIR="/scratch/evt82290/MappingOutputs/Run149/RNA/bamFiles"
OUTDIR="/scratch/evt82290/MappingOutputs/RNAseq/counts"

# create output dir if missing
[ -d "$OUTDIR" ] || mkdir -p "$OUTDIR"

# Modules
module load Subread
module load SAMtools
module load R

############################################
# Paired-end (PE) featureCounts  (-p)
############################################
featureCounts -T "$THREADS" \
  -p \
  -t CDS \
  -g gene_name \
  -s 0 --primary \
  -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  -o "$OUTDIR/readcounts_All_CAF1paper.txt" \
  $BAMDIR/SRR8444037/SRR8444037_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR8444038/SRR8444038_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR8444043/SRR8444043_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7970629/SRR7970629_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7970630/SRR7970630_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7970631/SRR7970631_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7970598/SRR7970598_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7970599/SRR7970599_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7970600/SRR7970600_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR8269830/SRR8269830_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR8269647/SRR8269647_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR8269650/SRR8269650_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR10916318/SRR10916318_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR10916319/SRR10916319_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR10916320/SRR10916320_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7970603/SRR7970603_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7970606/SRR7970606_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7970610/SRR7970610_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR9027727/SRR9027727_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR9027728/SRR9027728_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR9027730/SRR9027730_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR8269825/SRR8269825_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR8269775/SRR8269775_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR8269782/SRR8269782_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR8269810/SRR8269810_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR10916163/SRR10916163_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR10916164/SRR10916164_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR10916165/SRR10916165_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-115_RNA_cac-1__Rep1_S139/149-115_RNA_cac-1__Rep1_S139_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-116_RNA_cac-1__Rep2_S140/149-116_RNA_cac-1__Rep2_S140_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-117_RNA_cac-1__Rep3_S141/149-117_RNA_cac-1__Rep3_S141_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-65_RNA_cac-1_cac-2___S65/149-65_RNA_cac-1_cac-2___S65_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-66_RNA_cac-1_cac-2___S66/149-66_RNA_cac-1_cac-2___S66_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-67_RNA_cac-1_cac-2___S67/149-67_RNA_cac-1_cac-2___S67_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-96_RNA_cac-1-suz-12__Rep1_S96/149-96_RNA_cac-1-suz-12__Rep1_S96_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-97_RNA_cac-1-suz-12__Rep2_S121/149-97_RNA_cac-1-suz-12__Rep2_S121_Aligned.sortedByCoord.out.bam \
  $B149DIR/149-98_RNA_cac-1-suz-12__Rep3_S122/149-98_RNA_cac-1-suz-12__Rep3_S122_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR9027658/SRR9027658_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR9027759/SRR9027759_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR9027689/SRR9027689_Aligned.sortedByCoord.out.bam

############################################
# Single-end (SE) featureCounts  (no -p)
############################################
featureCounts -T "$THREADS" \
  -t CDS \
  -g gene_name \
  -s 0 --primary \
  -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  -o "$OUTDIR/readcounts_All_CAF1paper.SE.txt" \
  $BAMDIR/SRR7690267/SRR7690267_Aligned.sortedByCoord.out.bam \
  $BAMDIR/SRR7690268/SRR7690268_Aligned.sortedByCoord.out.bam

############################################
# Merge PE + SE (SE appended as last 2 cols)
############################################
PE="$OUTDIR/readcounts_All_CAF1paper.txt"
SE="$OUTDIR/readcounts_All_CAF1paper.SE.txt"
MERGED="$OUTDIR/readcounts_All_CAF1paper.merged.txt"

Rscript - "$PE" "$SE" "$MERGED" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
pe <- args[1]; se <- args[2]; out <- args[3]

read_fc <- function(f) {
  x <- read.table(f, header=TRUE, sep="\t", quote="",
                  comment.char="#", check.names=FALSE, stringsAsFactors=FALSE)
  colnames(x)[1] <- sub("^\ufeff", "", colnames(x)[1])  # strip BOM if present
  x
}

pe_df <- read_fc(pe)
se_df <- read_fc(se)

# Align SE rows to PE order, then append SE sample columns at the end
se_df <- se_df[match(pe_df$Geneid, se_df$Geneid), ]
if (!all(pe_df$Geneid == se_df$Geneid)) stop("Gene sets do not match between PE and SE.")

merged <- cbind(pe_df[, 1:6, drop=FALSE],
                pe_df[, 7:ncol(pe_df), drop=FALSE],
                se_df[, 7:ncol(se_df), drop=FALSE])

write.table(merged, out, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
RS

# quick header check
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) printf "%3d\t%s\n", i,$i; exit}' "$MERGED"
