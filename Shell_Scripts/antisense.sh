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
fastqPath="/scratch/evt82290/SRA/FastqFiles"
outdir="/scratch/evt82290/RNAseq/cac_aberrant_transcripts"

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
  read1=${fastqPath}/${accession}/${accession}_1.fastq.gz
  read2=${fastqPath}/${accession}/${accession}_2.fastq.gz
  unpaired=${fastqPath}/${accession}/${accession}.fastq.gz

###################################

#make output file folders
trimmed="${outdir}/TrimmedFastQs/${accession}"
mkdir $trimmed

bamdir="${outdir}/bamFiles/${accession}"
mkdir "${bamdir}"

countsdir="${outdir}/counts/${accession}"
mkdir "${countsdir}"

countsdir2="${outdir}/counts2/${accession}"
mkdir "${countsdir}"

bwDir="${outdir}/bigWig/${accession}"
mkdir "${bwDir}"

#pipeaccession summary: trim reads, map with STAR, get Counts

#notes
#generated a STAR genome index with the following call:
#STAR --runMode genomeGenerate --runThreadN 1 --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR --genomeFastaFiles /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO.fna --sjdbGTFfile /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_WithExtras_GFFtoGTFconversion.gtf
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
fowcounts="${countsdir}/${accession}_counts_fow.txt"
fowbw="${bwDir}/${accession}_fow.bw"
revcounts="${countsdir}/${accession}_counts_rev.txt"
revbw="${bwDir}/${accession}_rev.bw"

sencounts="${countsdir2}/${accession}_counts_sen.txt"
senbw="${bwDir}/${accession}_sen.bw"
antcounts="${countsdir2}/${accession}_counts_ant.txt"
antbw="${bwDir}/${accession}_ant.bw"

#modules
module load SAMtools/1.16.1-GCC-11.3.0
module load StringTie
module load BEDTools

#MAP READS TO SENSE AND ANTISENSE
# intersectBed -a genes.gff -b input.bam  -s  >  sense.bam
# intersectBed -a genes.gff -b input.bam  -S  >  antisense.bam

# first read in pair maps to reverse strand, read is mapped in proper pair, read is paired
samtools view -f 83 -b ${bam}Aligned.sortedByCoord.out.bam > ${bam}out.sorted.83.bam
samtools view -f 163 -b ${bam}Aligned.sortedByCoord.out.bam > ${bam}out.sorted.163.bam
# second read in pair maps to reverse strand, read is paired, read mapped to proper pair
samtools view -f 99 -b ${bam}Aligned.sortedByCoord.out.bam > ${bam}out.sorted.99.bam
samtools view -f 147 -b ${bam}Aligned.sortedByCoord.out.bam > ${bam}out.sorted.147.bam
# merge files
samtools merge ${bam}out.sorted.83.163.bam ${bam}out.sorted.83.bam ${bam}out.sorted.163.bam
samtools merge ${bam}out.sorted.99.147.bam ${bam}out.sorted.99.bam ${bam}out.sorted.147.bam
# stringtie to create a gtf reference - this is a non model organism
# stringtie ${bam}out.sorted.99.147.bam -o out.99.147.stringtie.gtf -p 3 --rf -m 50
# stringtie ${bam}out.sorted.83.163.bam -o out.83.163.stringtie.gtf -p 3 --rf -m 50
# stringtie merge out.83.163.stringtie.gtf out.99.147.stringtie.gtf -o stringtie_merged.gtf
# extract features that are on the + or - strand from stringtie
# cat /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf | awk '$7 == "+" { print $0 }' > forward_cds.gff
# cat /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf | awk '$7 == "-" { print $0 }' > reverse_cds.gff
# separate out reads mapping to positive and negative features from the stranded read files to get sense and antisense transcription
bedtools intersect -b forward_cds.gff -abam ${bam}out.sorted.83.163.bam > ${bam}out.sorted.83.163.sense.bam
bedtools intersect -b reverse_cds.gff -abam ${bam}out.sorted.83.163.bam > ${bam}out.sorted.83.163.antisense.bam
bedtools intersect -b forward_cds.gff -abam ${bam}out.sorted.99.147.bam > ${bam}out.sorted.99.147.antisense.bam
bedtools intersect -b reverse_cds.gff -abam ${bam}out.sorted.99.147.bam > ${bam}out.sorted.99.147.sense.bam
# merge the two sense and antisense files
samtools merge ${bam}out_flags_sense.bam ${bam}out.sorted.83.163.sense.bam ${bam}out.sorted.99.147.sense.bam
samtools merge ${bam}out_flags_antisense.bam ${bam}out.sorted.83.163.antisense.bam ${bam}out.sorted.99.147.antisense.bam
#index
samtools index ${bam}out_flags_sense.bam
samtools index ${bam}out_flags_antisense.bam
#bigwig
module load deepTools/3.5.2-foss-2022a
#Plot all reads
bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}out_flags_sense.bam" -o "${senbw}"
bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}out_flags_antisense.bam" -o "${antbw}"
