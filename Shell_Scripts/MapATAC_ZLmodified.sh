#!/bin/bash
#SBATCH --job-name=ET_ATAC.%j.job
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapATAC.%j.out
#SBATCH --error=../ATAC.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config_bwt.txt

OUTDIR=${OutputFolderName}
# mkdir ${OUTDIR}


if [ ! -e "${BWA_GENOME}.bwt" ]; then
  echo "ERROR: BWA index not found for ${BWA_GENOME}"
  exit 1
fi

# #process reads using trimGalore
#
ml Trim_Galore/0.6.7-GCCcore-11.2.0
trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
# #
FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
#
 # mkdir "${OUTDIR}/SortedBamFiles"
 # mkdir "${OUTDIR}/BigWigs"
 # mkdir "${OUTDIR}/hmmr_Peaks"
 # mkdir "${OUTDIR}/reg_Peaks"
 # mkdir "${OUTDIR}/Beds"
 # mkdir "${OUTDIR}/Histograms"
 # mkdir "${OUTDIR}/Metaplots"
mkdir "${OUTDIR}/Stats"

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
	name=${file/%_S[1-12]*_L003_R1_001_val_1.fq.gz/}

# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
	#variable for naming bam file
 	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	shifted="${OUTDIR}/SortedBamFiles/${name}.shifted.bam"
  	deduped="${OUTDIR}/SortedBamFiles/${name}_deduped.bam"
  	bed="${OUTDIR}/Beds/${name}.bed"
	#variable name for bigwig output
	bigwig="${OUTDIR}/BigWigs/${name}"
  meta="${OUTDIR}/Metaplots/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#
ml SAMtools
ml BWA

bwa mem -M -v 3 -a -t $THREADS $BWA_GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"
#

#### Collect insert sizes using picardtools and plot histogram
 module load picard

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
     I=$bam \
     O=${OUTDIR}/Histograms/${name}.fraglen.stats \
     H=${OUTDIR}/Histograms/${name}.fraglen.pdf \
     M=0.5 \
    VALIDATION_STRINGENCY=LENIENT QUIET=TRUE VERBOSITY=ERROR

#### Mark duplicate sequences to remove from analysis and plot insert sizes before and after marking duplicates
### Install Picard in home directory (replace "/home/ad45368/apps/picard/3.4.0/picard.jar" with path to your Picard installation)
### For installation directions, see here: https://wiki.gacrc.uga.edu/wiki/Installing_Applications_Sapelo2

module load Java
module load R
ml SAMtools

## Reformat bam file to be compatible with Picard - this next step is optional; include if Picard throws out an error

#samtools addreplacerg -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o ${OUTDIR}/SortedBamFiles/${name}.tmp.bam $bam

### Plot size histogram BEFORE removing MarkDuplicates
java -Xmx20g -classpath "/home/evt82290/Research/picard/3.4.0" -jar /home/ad45368/apps/picard/3.4.0/picard.jar CollectInsertSizeMetrics \
    --Histogram_FILE ${OUTDIR}/Histograms/${name}.fraglen.pdf \
     --INPUT $bam \
     --OUTPUT ${OUTDIR}/Stats/${name}.fraglen.stats \
     -M 0.5 \
    --VALIDATION_STRINGENCY SILENT --QUIET TRUE


### Remove duplicates
java -jar -Xmx20g -classpath "/home/evt82290/Research/picard/3.4.0" -jar /home/ad45368/apps/picard/3.4.0/picard.jar MarkDuplicates \
     --INPUT ${OUTDIR}/SortedBamFiles/${name}.tmp.bam \
     --OUTPUT ${deduped} \
     -M ${OUTDIR}/Stats/${name}.marked_dup_metrics.txt \
     --REMOVE_DUPLICATES TRUE --QUIET TRUE

#uncomment if you had to run samtools addreplacerg above to delete intermediate bam file
#rm ${OUTDIR}/SortedBamFiles/${name}.tmp.bam

### Removes 50% of reads --- not using deduplicated files.
### Plot sizes again with duplicates removed
java -Xmx20g -classpath "/home/evt82290/Research/picard/3.4.0" -jar /home/ad45368/apps/picard/3.4.0/picard.jar CollectInsertSizeMetrics \
         --Histogram_FILE ${OUTDIR}/Histograms/${name}.deduped.fraglen.pdf \
          --INPUT $deduped \
          --OUTPUT ${OUTDIR}/Stats/${name}.deduped.fraglen.stats \
          -M 0.5 \
         --VALIDATION_STRINGENCY SILENT --QUIET TRUE

#index again
samtools index -@ $THREADS ${deduped}


perl ./shiftTn5_BAM_2_BED.pl "${bam}" > "${name}.bed"
#
# ############################
# #deeptools
module load deepTools/3.5.2-foss-2022a
alignmentSieve -p $THREADS --ATACshift --bam ${bam} -o ${name}.tmp.bam

# the bam file needs to be sorted again
samtools sort -@ $THREADS -O bam -o ${shifted} ${name}.tmp.bam
samtools index -@ $THREADS ${shifted}
rm ${name}.tmp.bam
#
# #Plot all reads
bamCoverage -p $THREADS --Offset 1 3 -bs 3 --smoothLength 6 --normalizeUsing BPM  -of bigwig -b ${shifted} -o "${bigwig}.ATAC_bin_3.smooth_6_Bulk.bw"
bamCoverage -p $THREADS --Offset 1 3 -bs 1 --smoothLength 6 --normalizeUsing BPM  -of bigwig -b ${shifted} -o "${bigwig}.ATAC_bin_1.smooth_6_Bulk.bw"

#plot mononucleosomes
#bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4

macs3 hmmratac -b $shifted --outdir ${OUTDIR}/hmmr_Peaks -n $name
macs3 callpeak -t $shifted --format BAMPE --outdir ${OUTDIR}/reg_Peaks -n $name


done
