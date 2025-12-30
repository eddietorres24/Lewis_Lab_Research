#!/bin/bash
#SBATCH --job-name=zl_mapChIPseq
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=72:00:00
#SBATCH --output=../MapCutAndRun.%j.out
#SBATCH --error=../MapCutAndRun.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

#Make Output Directory
OUTDIR="/scratch/evt82290/Run144"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
###

# #process reads using trimGalore
#
ml Trim_Galore SAMtools BWA

bwa index $GENOME
# trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
# #
FILES="${FASTQ}/*fastq\.gz" #Don't forget the *
# #

#/137-66_ChIP_WT_H3K27me3_Rep1_S63_R1_001_val_1.fq.gz.bam

 mkdir "${OUTDIR}/SortedBamFiles"
 mkdir "${OUTDIR}/BigWigs"
 mkdir "${OUTDIR}/Peaks"
#mkdir "$OUTDIR/HomerTagDirectories"
#mkdir "$OUTDIR/TdfFiles"
#
#Iterate over the files
for f in $FILES
do
#
# 	#Examples to Get Different parts of the file name
# 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#${string//substring/replacement}
		# dir=${f%/*}

	file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
	name=${file/%_S[1-12]*_L002_R1_001_val_1.fq.gz/}


# 137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72_R1_001_val_1.fq.gz
# 137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72_R2_001_val_2.fq.gz
# /scratch/evt82290/Run137/SortedBamFiles/137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72_R1_001_val_1.fq.gz.bam

#
# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
	#variable for naming bam file
 	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	#variable name for bigwig output
	bigwig="${OUTDIR}/BigWigs/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"

ml SAMtools
ml BWA
#
bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"

#samtools view -b -q 30 $bam > "$QualityBam"
#samtools index "$QualityBam"

############################
# # #deeptools

ml deepTools
#Plot all reads
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes (don't need to do for ChIP)
#bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

#call Peaks
module load MACS3

#using --nolambda paramenter to call peaks without control
# macs3 callpeak -t "${bam}" -f BAMPE -n "${name}" --broad -g 41037538 --broad-cutoff 0.1 --outdir "${OUTDIR}/Peaks" --min-length 800 --max-gap 500 --nolambda
 done
