#!/bin/bash
#SBATCH --job-name=ET_mapChIPseq
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=24:00:00
#SBATCH --output=../MapCutAndRun.%j.out
#SBATCH --error=../MapCutAndRun.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

#Make output directory
OUTDIR="/scratch/evt82290/Run137"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
###


# #process reads using trimGalore
#
 ml Trim_Galore BWA SAMtools deepTools
 # trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
#
FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
#
 mkdir "${OUTDIR}/SortedBamFiles"
 mkdir "${OUTDIR}/BigWigs"
 mkdir "${OUTDIR}/Peaks"
#mkdir "$OUTDIR/HomerTagDirectories"
#mkdir "$OUTDIR/TdfFiles"


#Iterate over the files
for f in $FILES
do
#
# 	#Examples to Get Different parts of the file name
# 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#${string//substring/replacement}
# 		#dir=${f%/*}

	file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
	name=${file/%_S[1-12]*_L001_R1_001_val_1.fq.gz/}

#
# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
	#variable for naming bam file
 	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	#variable name for bigwig output
	bigwig="${OUTDIR}/BigWigs/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#

ml SAMtools
#
bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"

#samtools view -b -q 30 $bam > "$QualityBam"
#samtools index "$QualityBam"

############################
# # #deeptools

ml deepTools
# #use these parameters for ChIP data
bamCoverage -p $THREADS $MNase -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}${MN}.bw"

done
