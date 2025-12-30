#!/bin/bash
#SBATCH --job-name=RNAseq_Map
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --time=8:00:00
#SBATCH --output=../RNAseqMap/logs/MapRNAseq.%j.out
#SBATCH --error=../RNAseqMap/logs/MapRNAseq.%j.err

cd $SLURM_SUBMIT_DIR

THREADS=2

##ADD a source file with path to FastqFiles
#variables imported from submission script
#accession=SRR10916163
# fastqPath="/scratch/evt82290/SRA/FastqFiles"
# fastqPath=${Run150path}
# outdir="/scratch/evt82290/MappingOutputs/Run150"
# fastqPath="/scratch/evt82290/FastqFiles/2025_Run149_ET/RNA"
# outdir="/scratch/evt82290/MappingOutputs/Run149/RNA"

# #if output directory doesn't exist, create it
if [ ! -d $outdir ]
then
    mkdir -p $outdir
fi
# ###

###################
#start
####################################################

#input file variables
  # read1=${fastqPath}/${accession}*R1_001.fastq.gz
  # read2=${fastqPath}/${accession}*R2_001.fastq.gz
  read1=${fastqPath}/${accession}*_1.fastq.gz
  read2=${fastqPath}/${accession}*_2.fastq.gz
  # unpaired=${fastqPath}/${accession}/${accession}*.fastq.gz
#For files not in SRA folder
# read1=${fastqPath}/${accession}*_1.fastq.gz
# read2=${fastqPath}/${accession}*_2.fastq.gz
# unpaired=${fastqPath}/${accession}.fastq.gz

###################################

#make output file folders
trimmed="${outdir}/TrimmedFastQs/${accession}"
mkdir $trimmed

bamdir="${outdir}/bamFiles/${accession}"
mkdir "${bamdir}"

countsdir="${outdir}/counts/${accession}"
mkdir "${countsdir}"

bwDir="${outdir}/bigWig/${accession}"
mkdir "${bwDir}"

bedgraphDir="${outdir}/bedGraph/${accession}"
mkdir "${bedgraphDir}"

mkdir ${outdir}/condensed

#pipeaccession summary: trim reads, map with STAR, get Counts

#notes
#generated a STAR genome index with the following call:
#STAR --runMode genomeGenerate --runThreadN 1 --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR --genomeFastaFiles /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO.fna --sjdbGTFfile /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_WithExtras_GFFtoGTFconversion.gtf
#STAR --runMode genomeGenerate --runThreadN 1 --genomeDir /home/evt82290/Research --genomeFastaFiles /home/evt82290/Research/GCA_000182925.2_NC12_genomic_wTetO_at_his3_CLEAN.fasta --sjdbGTFfile /home/evt82290/Research/Nc12wTetO_at_his3_CLEAN.gff
#need to rerun with normal genome assembly. The his-3 duplicated region will create multi-mappers


# #make output file folders
# trimmed= mkdir "${outdir}/TrimmedFastQs/${accession}"
#
# bamdir= mkdir "${outdir}/bamFiles/${accession}"
#
# countsdir= mkdir "${outdir}/counts/${accession}"
#
# bwDir= mkdir "${outdir}/bigWig/${accession}"

#make variables for output file names
bam="${bamdir}/${accession}_"
counts="${countsdir}/${accession}_counts.txt"
bw="${bwDir}/${accession}.bw"
bg="${bedgraphDir}/${accession}.bedGraph"
bg2="${bedgraphDir}/${accession}_log2.bedGraph"
bg2_sort="${bedgraphDir}/${accession}_log2_sorted.bedGraph"
bw2="${bwDir}/${accession}_log.bw"

############# Read Trimming ##############
#remove adaptors, trim low quality reads (default = phred 20), length > 25

##fastq files from the ebi link are in folders that either have one file with a SRR##.fastq.gz or a SRR##_1.fastq.gz ending, or have two files with a SRR##_1.fastq.gz ending or a SRR##_2.fastq.gz ending
#or have three files with a SRR##_1.fastq.gz, SRR##_2.fastq.gz and SRR##.fastq.gz ending. In this case, the third file corresponds to unpaired reads that the depositers mapped.

#if read1 file does not exist, do single-end trimming using the only file in the folder i.e. SRR##.fastq.gz filename format
##This entire section can be simplified for JGI data

if [ ! -f $read1 ]; then
#trim reads
  echo "${line} running as unpaired file only"

  module load Trim_Galore

  trim_galore --illumina --fastqc --length 25 --basename ${accession} --gzip -o $trimmed $unpaired

  wait

#map with STAR
  module load STAR

  STAR --runMode alignReads \
  --runThreadN $THREADS \
  --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
  --outFileNamePrefix ${bam} \
  --readFilesIn $trimmed/${accession}_trimmed.fq.gz \
  --readFilesCommand zcat \
  --alignIntronMax 10000 \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortingBinsN 100 \
  --outSAMunmapped Within \
  --outSAMattributes Standard \
  --limitBAMsortRAM 19990000000

  #create index
  module load SAMtools
  samtools index "${bam}Aligned.sortedByCoord.out.bam"

  ##quantify with featureCounts
  module load Subread

  featureCounts -T $THREADS \
  -t CDS \
  -g gene_name \
  -s 0 --primary \
  -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  -o $counts \
  ${bam}Aligned.sortedByCoord.out.bam


##Plot reads to visualize tracks if needed
       module load deepTools
       #Plot all reads
       bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}Aligned.sortedByCoord.out.bam" -o "${bw}"

#log-transformed BigWigs
module load ucsc

bigWigToBedGraph ${bw} ${bg}
awk '{ $4=(log($4+1)/log(2)); } 1' < ${bg} > ${bg2}
sort -k1,1 -k2,2n ${bg2} > ${bg2_sort}
bedGraphToBigWig ${bg2_sort} /home/ad45368/chrom_sizes.txt ${bw2}



#elseif read2 exists, do paired-end Trimming and PE mapping
elif [ -f $read2 ]; then
  echo "${line} running as PE file only"


  ##################
  #Trimming
  #################
  	  module load Trim_Galore

  	  trim_galore --illumina --fastqc --paired --length 25 --basename ${accession} --gzip -o $trimmed $read1 $read2
  	  wait

  ##map with STAR
  	  module load STAR
  	    STAR --runMode alignReads \
  	    --runThreadN $THREADS \
  	    --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
  	    --outFileNamePrefix ${bam} \
  	    --readFilesIn $trimmed/${accession}_val_1.fq.gz $trimmed/${accession}_val_2.fq.gz \
  	    --readFilesCommand zcat \
        --alignIntronMax 10000 \
  	    --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingBinsN 100 \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --limitBAMsortRAM 19990000000
        #create index
        module load SAMtools
        samtools index "${bam}Aligned.sortedByCoord.out.bam"

        ##quantify with featureCounts
        module load Subread

        featureCounts -T $THREADS \
        -p \
        -t CDS \
        -g gene_name \
        -s 0 --primary \
        -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
        -o $counts \
        ${bam}Aligned.sortedByCoord.out.bam


        ##Plot reads to visualize tracks if needed
             module load deepTools

             #Plot all reads
             bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}Aligned.sortedByCoord.out.bam" -o "${bw}"

#log-transformed BigWigs
module load ucsc

bigWigToBedGraph ${bw} ${bg}
awk '{ $4=(log($4+1)/log(2)); } 1' < ${bg} > ${bg2}
sort -k1,1 -k2,2n ${bg2} > ${bg2_sort}
bedGraphToBigWig ${bg2_sort} /home/ad45368/chrom_sizes.txt ${bw2}

#in rare cases there will only be a SRR##_1.fastq.gz format. Use this if nothing else exists.
else

    echo "${accesion} running as Read1 file only"

module load Trim_Galore

       # trim_galore --illumina --fastqc --length 25 --basename ${accession} --gzip -o $trimmed $read1

       #map with STAR
       module load STAR
       #   STAR --runMode alignReads \
       #   --runThreadN $THREADS \
       #   --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
       #   --outFileNamePrefix ${bam} \
       #   --readFilesIn $trimmed/${accession}_trimmed.fq.gz  \
       #   --readFilesCommand zcat \
       #   --alignIntronMax 10000 \
       #   --outSAMtype BAM SortedByCoordinate \
       #   --outSAMunmapped Within \
       #   --outSAMattributes Standard \
       #   --outBAMsortingBinsN 100 \
       #   --limitBAMsortRAM 19990000000
       #
       #   #create index
         module load SAMtools
       #   samtools index "${bam}Aligned.sortedByCoord.out.bam"

         ##quantify with featureCounts
         module load Subread

         featureCounts -T $THREADS \
         -t CDS \
         -g gene_name \
         -s 0 --primary \
         -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
         -o $counts \
         ${bam}Aligned.sortedByCoord.out.bam


         ##Plot reads to visualize tracks if needed
         	    module load deepTools
         	    #Plot all reads
         	    bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}Aligned.sortedByCoord.out.bam" -o "${bw}"

              #log-transformed BigWigs
              module load ucsc

              bigWigToBedGraph ${bw} ${bg}
              awk '{ $4=(log($4+1)/log(2)); } 1' < ${bg} > ${bg2}
              sort -k1,1 -k2,2n ${bg2} > ${bg2_sort}
              bedGraphToBigWig ${bg2_sort} /home/ad45368/chrom_sizes.txt ${bw2}

fi


fi
