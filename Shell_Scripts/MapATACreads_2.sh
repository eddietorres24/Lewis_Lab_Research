#!/bin/bash
#SBATCH --job-name=ATACMap.%j.job
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=./logs/bwtMapATAC.%j.out
#SBATCH --error=./logs/ATAC.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config_bwt.txt

OUTDIR=${OutputFolderName}
# mkdir ${OUTDIR}


# #process reads using trimGalore
 # ml Trim_Galore/0.6.7-GCCcore-11.2.0
 # trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz

FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
#

 # mkdir "${OUTDIR}/TrimmedReads"
 # mkdir "${OUTDIR}/SortedBamFiles"
 # mkdir "${OUTDIR}/ShiftedBamFiles"
 # mkdir "${OUTDIR}/FilteredBamFiles"
 #
 # mkdir "${OUTDIR}/BigWigs"
 # mkdir "${OUTDIR}/Peaks"
 # mkdir "${OUTDIR}/logs"


mkdir "$OUTDIR/HomerTagDirectories"
mkdir "$OUTDIR/TdfFiles"
#
#Iterate over the files
for f in $FILES
do
#
# 	#Examples to Get Different parts of the file name
		# See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		# ${string//substring/replacement}
		# dir=${f%/*}

	file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
	name=${file/%_R1_001_val_1.fq.gz/}

# 	# File Vars

###***Differences b/w the older version of the script here. changed the folders where the variables assigned to***

# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$file" | sed 's/_R1_001_val_1\.fq\.gz/_R2_001_val_2\.fq\.gz/g')
  	#variable for naming bam file
  bwt_bam="${OUTDIR}/SortedBamFiles/${name}.bam"

  deduped1="${OUTDIR}/SortedBamFiles/${name}_deduped.bam"
  markeddupes="${OUTDIR}/SortedBamFiles/${name}_marked_dup_metrics.txt"

	shifted="${OUTDIR}/ShiftedBamFiles/${name}_shifted.bam"
  deduped2="${OUTDIR}/FilteredBamFiles/${name}_shifted_deduped.bam"

    nfr="${OUTDIR}/SortedFilteredBamFiles/${name}_nfr.bam"
    mono="${OUTDIR}/SortedFilteredBamFiles/${name}_mono.bam"
    di="${OUTDIR}/SortedFilteredBamFiles/${name}_di.bam"
    tri="${OUTDIR}/SortedFilteredBamFiles/${name}_tri.bam"
	#variable name for bigwig output
	bwdir="${OUTDIR}/BigWigs"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"

############ BWA Mapping #####################

ml SAMtools/1.16.1-GCC-11.3.0
ml BWA/0.7.17-GCCcore-11.3.0


### bowtie2 alignment -- works ###
module load Bowtie2/2.5.2-GCC-11.3.0
#
bowtie2 -p $THREADS -q --local --very-sensitive -x $BWT_GENOME -1 ${OUTDIR}/TrimmedReads/$file -2 ${OUTDIR}/TrimmedReads/$read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/Bowtie2/SortedBamFiles/tempReps -o "$bwt_bam"
samtools index "$bwt_bam"

##removing duplicates ##
module load picard/2.27.5-Java-15
# #
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
     -I ${bwt_bam} \
     -O ${deduped1} \
     -M ${markeddupes} \
     --REMOVE_DUPLICATES
samtools index "$deduped1"

###bwt_deduped -> deduped1

#deeptools
module load deepTools/3.5.2-foss-2022a
alignmentSieve -p $THREADS --ATACshift --bam $deduped1 -o ${OUTDIR}/ShiftedBamFiles/${name}_tmp.bam

#the bam file needs to be sorted again
samtools sort -@ $THREADS -O bam -o ${shifted} ${OUTDIR}/ShiftedBamFiles/${name}_tmp.bam
samtools index -@ $THREADS ${shifted}

#rm ${name}.tmp.bam


## filter bam to keep only unique mapping reads ###
module load Sambamba/0.8.2-GCC-11.3.0

sambamba view -h -t $THREADS -f bam \
-F "[XS] == null and not unmapped and not duplicate" \
${shifted} > ${deduped2}
samtools index -@ $THREADS ${deduped2}

#second filtering step to sort by fragment sizes#
sambamba view --format \
  bam --nthreads $THREADS \
  -F "((template_length > 0 and template_length < 100) or (template_length < 0 and template_length > -100))" $deduped2 | samtools view -b > $nfr
  samtools index -@ $THREADS ${nfr}


  sambamba view --format \
    bam --nthreads $THREADS \
    -F "((template_length > 100 and template_length < 200) or (template_length < -100 and template_length > -200))" $deduped2 | samtools view -b > $mono
    samtools index -@ $THREADS ${mono}

    sambamba view --format \
      bam --nthreads $THREADS \
      -F "((template_length > 0 and template_length < 400) or (template_length < -200 and template_length > -400))" $deduped2 | samtools view -b > $di
      samtools index -@ $THREADS ${di}

      sambamba view --format \
        bam --nthreads $THREADS \
        -F "((template_length > 0 and template_length < 600) or (template_length < -400 and template_length > -600))" $deduped2 | samtools view -b > $tri
        samtools index -@ $THREADS ${tri}


#Plot all reads
bamCoverage -p $THREADS --Offset 1 3 -bs 3 --smoothLength 6 --minMappingQuality 20 --normalizeUsing BPM  -of bigwig -b ${nfr} -o "${bwdir}/${name}_nfr_ATAC_bin_3_smooth_6_Bulk.bw"
bamCoverage -p $THREADS --Offset 1 3 -bs 3 --smoothLength 6 --minMappingQuality 20 --normalizeUsing BPM  -of bigwig -b ${mono} -o "${bwdir}/${name}_mono_ATAC_bin_3_smooth_6_Bulk.bw"
bamCoverage -p $THREADS --Offset 1 3 -bs 3 --smoothLength 6 --minMappingQuality 20 --normalizeUsing BPM  -of bigwig -b ${di} -o "${bwdir}/${name}_di_ATAC_bin_3_smooth_6_Bulk.bw"
bamCoverage -p $THREADS --Offset 1 3 -bs 3 --smoothLength 6 --minMappingQuality 20 --normalizeUsing BPM  -of bigwig -b ${tri} -o "${bwdir}/${name}_tri_ATAC_bin_3_smooth_6_Bulk.bw"
bamCoverage -p $THREADS --Offset 1 3 -bs 3 --smoothLength 6 --minMappingQuality 20 --normalizeUsing BPM  -of bigwig -b  ${deduped2} -o "${bwdir}/${name}_all_ATAC_bin_3_smooth_6_Bulk.bw"
#

module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4
#
#
macs3 hmmratac -b $deduped2 --outdir ${OUTDIR}/Peaks -n $name

#plot mononucleosomes
# ?\bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"
### call atac peaks ###
done
### merging replicates ###
