#!/bin/bash
#SBATCH --job-name=Lewislab_ATAC.%j.job
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
ml GCC
ml GCCcore

source config_Nc.txt

OUTDIR=../${OutputFolderName}
mkdir ${OUTDIR}
mkdir ${OUTDIR}/TrimmedReads

# #process reads using trimGalore
#
ml Trim_Galore
trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
# # #
FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
#
 #mkdir "${OUTDIR}/SortedBamFiles"
 mkdir "${OUTDIR}/BigWigs"
 mkdir "${OUTDIR}/Peaks/hmmratac"
 mkdir "${OUTDIR}/Peaks/BroadPeaks"
# mkdir "${OUTDIR}/Histograms"
# mkdir "${OUTDIR}/Metaplots"
 mkdir "${OUTDIR}/Stats"
# mkdir "${OUTDIR}/DedupedBamFiles"
bamdir="${OUTDIR}/BamFiles"
mkdir ${bamdir}
mkdir "${bamdir}/RawBams"
mkdir "${bamdir}/DownsampledBams"
mkdir "${bamdir}/MarkedDupes"
mkdir "${bamdir}/Deduped"
mkdir "${bamdir}/Shifted"

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
# 		#dir=${f%/*}

	file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
  name=$(echo "$file" | sed -E 's/(Rep[0-9]+).*/\1/')

#
# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
	#variable for naming bam FILES

 	bam="${bamdir}/RawBams/${name}.bam"

  downsampled_bam="${bamdir}/DownsampledBams/${name}_12M.bam"

  marked_dupes="${bamdir}/MarkedDupes/${name}_marked.dupe.bam"
  downsize_marked_dupes="${bamdir}/MarkedDupes/${name}_12M_marked.dupe.bam"

  deduped="${bamdir}/Deduped/${name}_dedupe.bam"
  downsize_deduped="${bamdir}/Deduped/${name}_12M_dedupe.bam"

  shifted_marked_dupes="${bamdir}/Shifted/${name}_shifted.marked.dupe.bam"
  shifted_downsized_marked_dupes="${bamdir}/Shifted/${name}_shifted.marked.dupe._12M.bam"
  shifted_dedupes="${bamdir}/Shifted/${name}_shifted.dedupe.bam"
  shifted_downsized_dedupes="${bamdir}/Shifted/${name}_shifted.marked.dupe.12M.bam"
  shifted_raw="${bamdir}/Shifted/${name}_raw.shifted.bam"
  shifted_downsize="${bamdir}/Shifted/${name}_shifted.12M.bam"

  	#variable name for bigwig files. will make names here since all will be mapped with mnase + bulk in different folders
  bwdir="${OUTDIR}/BigWigs"
  mkdir ${bwdir}
  mkdir ${bwdir}/MNase
  mkdir ${bwdir}/Bin15_Smooth45

  bw_mnase="${bwdir}/MNase/${name}"
  bw_bulk="${bwdir}/Bin15_Smooth45/${name}"

	bigwig="${OUTDIR}/BigWigs/${name}"
  meta="${OUTDIR}/Metaplots/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#
mkdir ${OUTDIR}/tmp
tmp="$OUTDIR/tmp/${name}"

 ml SAMtools
 ml BWA
 ml picard
 ml R
# # #
#
bwa mem -M -v 3 -a -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/BamFiles/tempReps -o "$bam" -
samtools index "$bam"

###### Processing RAW bam file ############
# ##mark dupes
 samtools addreplacerg -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o ${tmp}_raw.bam $bam
# # #
# # #
 java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      --INPUT ${tmp}_raw.bam \
      --OUTPUT ${marked_dupes} \
      -M ${OUTDIR}/Stats/${name}.marked_dup_metrics.txt \
      --QUIET TRUE

rm ${tmp}_raw.bam
 #index again
 samtools index -@ $THREADS ${marked_dupes}
#  ## collect stats
  java -jar $EBROOTPICARD/picard.jar CollectMultipleMetrics \
       --INPUT $marked_dupes \
       --OUTPUT ${OUTDIR}/Stats/${name}.multi \
       -R /home/ad45368/Af_Genomes/Afum_A1163_edit6.fasta \
     --VALIDATION_STRINGENCY SILENT --QUIET TRUE
#
#shift file
  module load deepTools
  alignmentSieve -p $THREADS --ATACshift --bam ${marked_dupes} -o ${tmp}_mark.tmp.bam
  samtools sort -@ $THREADS -O bam -o ${shifted_marked_dupes} ${tmp}_mark.tmp.bam
  samtools index -@ $THREADS ${shifted_marked_dupes}
  rm ${tmp}_mark.tmp.bam

#make BigWigs

 # bamCoverage -p $THREADS --Offset 1 3 -bs 15 --smoothLength 45 --minMappingQuality 25 --effectiveGenomeSize 317634770 --normalizeUsing BPM  -of bigwig -b ${shifted_downsize} -o "${bw_bulk}_marked.12M.bin15_smooth45_Bulk.bw" --blackListFileName /lustre2/scratch/ad45368/Af_ATAC/Run152/MapATAC_Output/BroadPeaks/152_N6_Genomic_CEA10_WT__Rep_2_S102_L006_R1_001_val_1.fq.gz_peaks.broadPeak
 # bamCoverage -p $THREADS --Offset 1 3 -bs 1 --smoothLength 6 --MNase --minMappingQuality 25 --effectiveGenomeSize 317634770 --normalizeUsing BPM  -of bigwig -b ${shifted_downsize} -o "${bw_mnase}_marked.12M.bin1_smooth6_MNase.bw" --blackListFileName /lustre2/scratch/ad45368/Af_ATAC/Run152/MapATAC_Output/BroadPeaks/152_N6_Genomic_CEA10_WT__Rep_2_S102_L006_R1_001_val_1.fq.gz_peaks.broadPeak
# ##setting fragsize
 bamCoverage -p $THREADS --Offset 1 3 -bs 15 --smoothLength 45 --minMappingQuality 25 --effectiveGenomeSize 41037538 --minFragmentLength 30 --maxFragmentLength 200 --normalizeUsing BPM  -of bigwig -b ${shifted_marked_dupes} -o "${bw_bulk}_marked.min30.max200.bin15_smooth45_Bulk.bw"
 bamCoverage -p $THREADS --Offset 1 3 -bs 1 --smoothLength 6 --MNase --minMappingQuality 25 --effectiveGenomeSize 41037538 --normalizeUsing BPM --minFragmentLength 30 --maxFragmentLength 200  -of bigwig -b ${shifted_marked_dupes} -o "${bw_mnase}_marked.min30.max200.bin1_smooth6_MNase.bw"

#control="/lustre2/scratch/ad45368/Af_ATAC/Run152/MapATAC_Output/SortedBamFiles"
module load MACS3
macs3 callpeak -t $shifted_marked_dupes --nolambda --format BAMPE --outdir ${OUTDIR}/Peaks/BroadPeaks -n ${name}_bampe --broad-cutoff 0.05 --broad --keep-dup all -q 0.01 -g 41037538
# macs3 hmmratac -i ${shifted_downsized_marked_dupes} -f BAMPE -n ${name}_hmmr --outdir Peaks/hmmratac -e /lustre2/scratch/ad45368/Af_ATAC/Run152/MapATAC_Output/BroadPeaks/152_N6_Genomic_CEA10_WT__Rep_2_S102_L006_R1_001_val_1.fq.gz_peaks.broadPeak



done

module load MultiQC

multiqc --no-ai ${OUTDIR} -i "Nc_ATAC_qc"
