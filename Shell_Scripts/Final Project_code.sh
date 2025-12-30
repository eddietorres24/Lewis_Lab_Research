#!/bin/bash
#SBATCH --job-name=histone_chaperone_plot 	            # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=168:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/scratch/evt82290/log.%j			          # Standard output and error log - # replace cbergman with your myid
#SBATCH --error=../ChaperonePlot.%j.err                  # error output
#SBATCH --mail-user=evt82290@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/scratch/evt82290/histone_chaperone_plot_FINAL"
FASTDIR="/scratch/evt82290/downSRA/FastqFiles"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
###

#loading modules
module load SRA-Toolkit/3.0.1-centos_linux64 BWA/0.7.17-GCCcore-11.3.0 SAMtools/1.16.1-GCC-11.3.0 Subread/2.0.6-GCC-11.3.0

# #downloading reference genome
# curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Neurospora_crassa/latest_assembly_versions/GCA_000182925.2_NC12/GCA_000182925.2_NC12_genomic.fna.gz  | gunzip -c > ${OUTDIR}/NC12_genome.fna
# curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/925/GCA_000182925.2_NC12/GCA_000182925.2_NC12_genomic.gtf.gz | gunzip -c > ${OUTDIR}/NC12_annotation.gtf
#
# #downloading fastq files using accession numbers
prefetch -O ${OUTDIR} SRR8444037 SRR8444038 SRR8444043 SRR7970629 SRR7970630 SRR7970631 SRR7970598 SRR7970599 SRR7970600 SRR8269830 SRR8269647 SRR8269650 SRR10916318 SRR10916319 SRR10916320 SRR7970603 SRR7970606 SRR7970610 SRR9027727 SRR9027728 SRR9027730 SRR8269825 SRR8269775	SRR8269782 SRR8269810 SRR10916163	SRR10916164	SRR10916165 SRR10916326 SRR10916324 SRR10916325 SRR8269628	SRR8269763	SRR8269811	SRR8269812 SRR7690267 SRR7690268
#
# #prefetch -O ${OUTDIR} SRR8444037 SRR8444038 SRR8444043 SRR7970629 SRR7970630 SRR7970631 SRR7970598 SRR7970599 SRR7970600 SRR8269825 SRR8269775 SRR8269782 SRR8269810 SRR10916163 SRR10916164 SRR10916165
#
# # cac-1
# # SRR8444037 SRR8444038 SRR8444043
#
# #rtt109
# #SRR10916326 SRR10916324 SRR10916325
#
# #isw
# #SRR8269628	SRR8269763	SRR8269811	SRR8269812
#
# cac-1
fastq-dump --split-files --gzip ${OUTDIR}/SRR8444037/SRR8444037.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR8444038/SRR8444038.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR8444043/SRR8444043.sra -O ${OUTDIR}
#
# cac-2
fastq-dump --split-files --gzip ${OUTDIR}/SRR7970629/SRR7970629.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR7970630/SRR7970630.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR7970631/SRR7970631.sra -O ${OUTDIR}
#
# cac-3
fastq-dump --split-files --gzip ${OUTDIR}/SRR7970598/SRR7970598.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR7970599/SRR7970599.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR7970600/SRR7970600.sra -O ${OUTDIR}
#
# fastq-dump --split-files --gzip ${OUTDIR}/SRR8269830/SRR8269830.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR8269647/SRR8269647.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR8269650/SRR8269650.sra -O ${OUTDIR}
#
# fastq-dump --split-files --gzip ${OUTDIR}/SRR10916318/SRR10916318.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR10916319/SRR10916319.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR10916320/SRR10916320.sra -O ${OUTDIR}
# #
#  fastq-dump --split-files --gzip ${OUTDIR}/SRR10916324.sra -O ${OUTDIR}
#  fastq-dump --split-files --gzip ${OUTDIR}/SRR10916325.sra -O ${OUTDIR}
#  fastq-dump --split-files --gzip ${OUTDIR}/SRR10916326.sra -O ${OUTDIR}
# #
# fastq-dump --split-files --gzip ${OUTDIR}/SRR7970603/SRR7970603.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR7970606/SRR7970606.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR7970610/SRR7970610.sra -O ${OUTDIR}
#
# fastq-dump --split-files --gzip ${OUTDIR}/SRR9027727/SRR9027727.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR9027728/SRR9027728.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR9027730/SRR9027730.sra -O ${OUTDIR}
#
# WT
fastq-dump --split-files --gzip ${OUTDIR}/SRR8269825/SRR8269825.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR8269775/SRR8269775.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR8269782/SRR8269782.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR8269810/SRR8269810.sra -O ${OUTDIR}
#
#set-7
fastq-dump --split-files --gzip ${OUTDIR}/SRR10916163/SRR10916163.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR10916164/SRR10916164.sra -O ${OUTDIR}
fastq-dump --split-files --gzip ${OUTDIR}/SRR10916165/SRR10916165.sra -O ${OUTDIR}
#
# fastq-dump --split-files --gzip ${OUTDIR}/SRR7690267/SRR7690267.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR7690268/SRR7690268.sra -O ${OUTDIR}
#
# fastq-dump --split-files --gzip ${OUTDIR}/SRR8269628.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR8269763.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR8269811.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/SRR8269812.sra -O ${OUTDIR}
#
# #making bwa index for neurospora genome
# bwa index ${OUTDIR}/NC12_genome.fna
#
# #mapping to neurospora genome, generating sorted .bam files
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8444037/SRR8444037_1.fastq.gz ${OUTDIR}/SRR8444037/SRR8444037_2.fastq.gz > ${OUTDIR}/SRR8444037.sam
# samtools view ${OUTDIR}/SRR8444037.sam -O BAM -o ${OUTDIR}/SRR8444037.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8444037.bam -o ${OUTDIR}/SRR8444037.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8444038/SRR8444038_1.fastq.gz ${OUTDIR}/SRR8444038/SRR8444038_2.fastq.gz > ${OUTDIR}/SRR8444038.sam
# samtools view ${OUTDIR}/SRR8444038.sam -O BAM -o ${OUTDIR}/SRR8444038.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8444038.bam -o ${OUTDIR}/SRR8444038.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8444043/SRR8444043_1.fastq.gz ${OUTDIR}/SRR8444043/SRR8444043_2.fastq.gz > ${OUTDIR}/SRR8444043.sam
# samtools view ${OUTDIR}/SRR8444043.sam -O BAM -o ${OUTDIR}/SRR8444043.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8444043.bam -o ${OUTDIR}/SRR8444043.sorted.bam
#
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7970629/SRR7970629_1.fastq.gz ${OUTDIR}/SRR7970629/SRR7970629_2.fastq.gz > ${OUTDIR}/SRR7970629.sam
# samtools view ${OUTDIR}/SRR7970629.sam -O BAM -o ${OUTDIR}/SRR7970629.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7970629.bam -o ${OUTDIR}/SRR7970629.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7970630/SRR7970630_1.fastq.gz ${OUTDIR}/SRR7970630/SRR7970630_2.fastq.gz > ${OUTDIR}/SRR7970630.sam
# samtools view ${OUTDIR}/SRR7970630.sam -O BAM -o ${OUTDIR}/SRR7970630.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7970630.bam -o ${OUTDIR}/SRR7970630.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7970631/SRR7970631_1.fastq.gz ${OUTDIR}/SRR7970631/SRR7970631_2.fastq.gz > ${OUTDIR}/SRR7970631.sam
# samtools view ${OUTDIR}/SRR7970631.sam -O BAM -o ${OUTDIR}/SRR7970631.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7970631.bam -o ${OUTDIR}/SRR7970631.sorted.bam
#
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7970598/SRR7970598_1.fastq.gz ${OUTDIR}/SRR7970598/SRR7970598_2.fastq.gz > ${OUTDIR}/SRR7970598.sam
# samtools view ${OUTDIR}/SRR7970598.sam -O BAM -o ${OUTDIR}/SRR7970598.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7970598.bam -o ${OUTDIR}/SRR7970598.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7970599/SRR7970599_1.fastq.gz ${OUTDIR}/SRR7970599/SRR7970599_2.fastq.gz > ${OUTDIR}/SRR7970599.sam
# samtools view ${OUTDIR}/SRR7970599.sam -O BAM -o ${OUTDIR}/SRR7970599.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7970599.bam -o ${OUTDIR}/SRR7970599.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7970600/SRR7970600_1.fastq.gz ${OUTDIR}/SRR7970600/SRR7970600_2.fastq.gz > ${OUTDIR}/SRR7970600.sam
# samtools view ${OUTDIR}/SRR7970600.sam -O BAM -o ${OUTDIR}/SRR7970600.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7970600.bam -o ${OUTDIR}/SRR7970600.sorted.bam
#
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269830/SRR8269830_1.fastq.gz ${OUTDIR}/SRR8269830/SRR8269830_2.fastq.gz > ${OUTDIR}/SRR8269830.sam
# samtools view ${OUTDIR}/SRR8269830.sam -O BAM -o ${OUTDIR}/SRR8269830.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269830.bam -o ${OUTDIR}/SRR8269830.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269647/SRR8269647_1.fastq.gz ${OUTDIR}/SRR8269647/SRR8269647_2.fastq.gz > ${OUTDIR}/SRR8269647.sam
# samtools view ${OUTDIR}/SRR8269647.sam -O BAM -o ${OUTDIR}/SRR8269647.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269647.bam -o ${OUTDIR}/SRR8269647.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269650/SRR8269650_1.fastq.gz ${OUTDIR}/SRR8269650/SRR8269650_2.fastq.gz > ${OUTDIR}/SRR8269650.sam
# samtools view ${OUTDIR}/SRR8269650.sam -O BAM -o ${OUTDIR}/SRR8269650.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269650.bam -o ${OUTDIR}/SRR8269650.sorted.bam
#
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR10916318/SRR10916318_1.fastq.gz ${OUTDIR}/SRR10916318/SRR10916318_2.fastq.gz > ${OUTDIR}/SRR10916318.sam
# samtools view ${OUTDIR}/SRR10916318.sam -O BAM -o ${OUTDIR}/SRR10916318.bam
# samtools sort --threads 6 ${OUTDIR}/SRR10916318.bam -o ${OUTDIR}/SRR10916318.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR10916319/SRR10916319_1.fastq.gz ${OUTDIR}/SRR10916319/SRR10916319_2.fastq.gz > ${OUTDIR}/SRR10916319.sam
# samtools view ${OUTDIR}/SRR10916319.sam -O BAM -o ${OUTDIR}/SRR10916319.bam
# samtools sort --threads 6 ${OUTDIR}/SRR10916319.bam -o ${OUTDIR}/SRR10916319.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR10916320/SRR10916320_1.fastq.gz ${OUTDIR}/SRR10916320/SRR10916320_2.fastq.gz > ${OUTDIR}/SRR10916320.sam
# samtools view ${OUTDIR}/SRR10916320.sam -O BAM -o ${OUTDIR}/SRR10916320.bam
# samtools sort --threads 6 ${OUTDIR}/SRR10916320.bam -o ${OUTDIR}/SRR10916320.sorted.bam
#
#
# # bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR10916324_1.fastq.gz ${OUTDIR}/SRR10916324_2.fastq.gz > ${OUTDIR}/SRR10916324.sam
# # samtools view ${OUTDIR}/SRR10916324.sam -O BAM -o ${OUTDIR}/SRR10916324.bam
# # samtools sort --threads 6 ${OUTDIR}/SRR10916324.bam -o ${OUTDIR}/SRR10916324.sorted.bam
# #
# # bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR10916325_1.fastq.gz ${OUTDIR}/SRR10916325_2.fastq.gz > ${OUTDIR}/SRR10916325.sam
# # samtools view ${OUTDIR}/SRR10916325.sam -O BAM -o ${OUTDIR}/SRR10916325.bam
# # samtools sort --threads 6 ${OUTDIR}/SRR10916325.bam -o ${OUTDIR}/SRR10916325.sorted.bam
# #
# # bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR10916326_1.fastq.gz ${OUTDIR}/SRR10916326_2.fastq.gz > ${OUTDIR}/SRR10916326.sam
# # samtools view ${OUTDIR}/SRR10916326.sam -O BAM -o ${OUTDIR}/SRR10916326.bam
# # samtools sort --threads 6 ${OUTDIR}/SRR10916326.bam -o ${OUTDIR}/SRR10916326.sorted.bam
#
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7970603/SRR7970603_1.fastq.gz ${OUTDIR}/SRR7970603/SRR7970603_2.fastq.gz > ${OUTDIR}/SRR7970603.sam
# samtools view ${OUTDIR}/SRR7970603.sam -O BAM -o ${OUTDIR}/SRR7970603.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7970603.bam -o ${OUTDIR}/SRR7970603.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7970606/SRR7970606_1.fastq.gz ${OUTDIR}/SRR7970606/SRR7970606_2.fastq.gz > ${OUTDIR}/SRR7970606.sam
# samtools view ${OUTDIR}/SRR7970606.sam -O BAM -o ${OUTDIR}/SRR7970606.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7970606.bam -o ${OUTDIR}/SRR7970606.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7970610/SRR7970610_1.fastq.gz ${OUTDIR}/SRR7970610/SRR7970610_2.fastq.gz > ${OUTDIR}/SRR7970610.sam
# samtools view ${OUTDIR}/SRR7970610.sam -O BAM -o ${OUTDIR}/SRR7970610.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7970610.bam -o ${OUTDIR}/SRR7970610.sorted.bam
#
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR9027727/SRR9027727_1.fastq.gz ${OUTDIR}/SRR9027727/SRR9027727_2.fastq.gz > ${OUTDIR}/SRR9027727.sam
# samtools view ${OUTDIR}/SRR9027727.sam -O BAM -o ${OUTDIR}/SRR9027727.bam
# samtools sort --threads 6 ${OUTDIR}/SRR9027727.bam -o ${OUTDIR}/SRR9027727.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR9027728/SRR9027728_1.fastq.gz ${OUTDIR}/SRR9027728/SRR9027728_2.fastq.gz > ${OUTDIR}/SRR9027728.sam
# samtools view ${OUTDIR}/SRR9027728.sam -O BAM -o ${OUTDIR}/SRR9027728.bam
# samtools sort --threads 6 ${OUTDIR}/SRR9027728.bam -o ${OUTDIR}/SRR9027728.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR9027730/SRR9027730_1.fastq.gz ${OUTDIR}/SRR9027730/SRR9027730_2.fastq.gz > ${OUTDIR}/SRR9027730.sam
# samtools view ${OUTDIR}/SRR9027730.sam -O BAM -o ${OUTDIR}/SRR9027730.bam
# samtools sort --threads 6 ${OUTDIR}/SRR9027730.bam -o ${OUTDIR}/SRR9027730.sorted.bam
#
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269825/SRR8269825_1.fastq.gz ${OUTDIR}/SRR8269825/SRR8269825_2.fastq.gz > ${OUTDIR}/SRR8269825.sam
# samtools view ${OUTDIR}/SRR8269825.sam -O BAM -o ${OUTDIR}/SRR8269825.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269825.bam -o ${OUTDIR}/SRR8269825.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269775/SRR8269775_1.fastq.gz ${OUTDIR}/SRR8269775/SRR8269775_2.fastq.gz > ${OUTDIR}/SRR8269775.sam
# samtools view ${OUTDIR}/SRR8269775.sam -O BAM -o ${OUTDIR}/SRR8269775.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269775.bam -o ${OUTDIR}/SRR8269775.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269782/SRR8269782_1.fastq.gz ${OUTDIR}/SRR8269782/SRR8269782_2.fastq.gz > ${OUTDIR}/SRR8269782.sam
# samtools view ${OUTDIR}/SRR8269782.sam -O BAM -o ${OUTDIR}/SRR8269782.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269782.bam -o ${OUTDIR}/SRR8269782.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269810/SRR8269810_1.fastq.gz ${OUTDIR}/SRR8269810/SRR8269810_2.fastq.gz > ${OUTDIR}/SRR8269810.sam
# samtools view ${OUTDIR}/SRR8269810.sam -O BAM -o ${OUTDIR}/SRR8269810.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269810.bam -o ${OUTDIR}/SRR8269810.sorted.bam
#
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR10916163/SRR10916163_1.fastq.gz ${OUTDIR}/SRR10916163/SRR10916163_2.fastq.gz > ${OUTDIR}/SRR10916163.sam
# samtools view ${OUTDIR}/SRR10916163.sam -O BAM -o ${OUTDIR}/SRR10916163.bam
# samtools sort --threads 6 ${OUTDIR}/SRR10916163.bam -o ${OUTDIR}/SRR10916163.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR10916164/SRR10916164_1.fastq.gz ${OUTDIR}/SRR10916164/SRR10916164_2.fastq.gz > ${OUTDIR}/SRR10916164.sam
# samtools view ${OUTDIR}/SRR10916164.sam -O BAM -o ${OUTDIR}/SRR10916164.bam
# samtools sort --threads 6 ${OUTDIR}/SRR10916164.bam -o ${OUTDIR}/SRR10916164.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR10916165/SRR10916165_1.fastq.gz ${OUTDIR}/SRR10916165/SRR10916165_2.fastq.gz > ${OUTDIR}/SRR10916165.sam
# samtools view ${OUTDIR}/SRR10916165.sam -O BAM -o ${OUTDIR}/SRR10916165.bam
# samtools sort --threads 6 ${OUTDIR}/SRR10916165.bam -o ${OUTDIR}/SRR10916165.sorted.bam
#
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7690267/SRR7690267_1.fastq.gz  > ${OUTDIR}/SRR7690267.sam
# samtools view ${OUTDIR}/SRR7690267.sam -O BAM -o ${OUTDIR}/SRR7690267.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7690267.bam -o ${OUTDIR}/SRR7690267.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR7690268/SRR7690268_1.fastq.gz  > ${OUTDIR}/SRR7690268.sam
# samtools view ${OUTDIR}/SRR7690268.sam -O BAM -o ${OUTDIR}/SRR7690268.bam
# samtools sort --threads 6 ${OUTDIR}/SRR7690268.bam -o ${OUTDIR}/SRR7690268.sorted.bam

# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269628_1.fastq.gz ${OUTDIR}/SRR8269628_2.fastq.gz > ${OUTDIR}/SRR8269628.sam
# samtools view ${OUTDIR}/SRR8269628.sam -O BAM -o ${OUTDIR}/SRR8269628.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269628.bam -o ${OUTDIR}/SRR8269628.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269763_1.fastq.gz ${OUTDIR}/SRR8269763_2.fastq.gz > ${OUTDIR}/SRR8269763.sam
# samtools view ${OUTDIR}/SRR8269763.sam -O BAM -o ${OUTDIR}/SRR8269763.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269763.bam -o ${OUTDIR}/SRR8269763.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269811_1.fastq.gz ${OUTDIR}/SRR8269811_2.fastq.gz > ${OUTDIR}/SRR8269811.sam
# samtools view ${OUTDIR}/SRR8269811.sam -O BAM -o ${OUTDIR}/SRR8269811.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269811.bam -o ${OUTDIR}/SRR8269811.sorted.bam
#
# bwa mem -t 6 ${OUTDIR}/NC12_genome.fna ${OUTDIR}/SRR8269812_1.fastq.gz ${OUTDIR}/SRR8269812_2.fastq.gz > ${OUTDIR}/SRR8269812.sam
# samtools view ${OUTDIR}/SRR8269812.sam -O BAM -o ${OUTDIR}/SRR8269812.bam
# samtools sort --threads 6 ${OUTDIR}/SRR8269812.bam -o ${OUTDIR}/SRR8269812.sorted.bam

#making matrix using featurecounts
#featureCounts -T 6 -p -a ${OUTDIR}/NC12_annotation.gtf -o ${OUTDIR}/readcounts_FINAL.txt ${OUTDIR}/SRR8444037.sorted.bam ${OUTDIR}/SRR8444038.sorted.bam ${OUTDIR}/SRR8444043.sorted.bam ${OUTDIR}/SRR7970629.sorted.bam ${OUTDIR}/SRR7970630.sorted.bam ${OUTDIR}/SRR7970631.sorted.bam ${OUTDIR}/SRR7970598.sorted.bam ${OUTDIR}/SRR7970599.sorted.bam ${OUTDIR}/SRR7970600.sorted.bam ${OUTDIR}/SRR8269830.sorted.bam ${OUTDIR}/SRR8269647.sorted.bam ${OUTDIR}/SRR8269650.sorted.bam ${OUTDIR}/SRR10916318.sorted.bam ${OUTDIR}/SRR10916319.sorted.bam ${OUTDIR}/SRR10916320.sorted.bam ${OUTDIR}/SRR7970603.sorted.bam ${OUTDIR}/SRR7970606.sorted.bam ${OUTDIR}/SRR7970610.sorted.bam ${OUTDIR}/SRR9027727.sorted.bam ${OUTDIR}/SRR9027728.sorted.bam ${OUTDIR}/SRR9027730.sorted.bam ${OUTDIR}/SRR8269825.sorted.bam ${OUTDIR}/SRR8269775.sorted.bam ${OUTDIR}/SRR8269782.sorted.bam ${OUTDIR}/SRR8269810.sorted.bam ${OUTDIR}/SRR10916163.sorted.bam ${OUTDIR}/SRR10916164.sorted.bam ${OUTDIR}/SRR10916165.sorted.bam

#featureCounts -T 6 -a ${OUTDIR}/NC12_annotation.gtf -o ${OUTDIR}/readcounts_ash3.txt ${OUTDIR}/SRR7690267.sorted.bam ${OUTDIR}/SRR7690268.sorted.bam
#featureCounts -T 6 -a ${OUTDIR}/NC12_annotation.gtf -o ${OUTDIR}/readcounts_rtt109.txt ${OUTDIR}/SRR10916326.sorted.bam ${OUTDIR}/SRR10916324.sorted.bam ${OUTDIR}/SRR10916325.sorted.bam

#isw
#${OUTDIR}/SRR8269628.sorted.bam ${OUTDIR}/SRR8269763.sorted.bam ${OUTDIR}/SRR8269811.sorted.bam ${OUTDIR}/SRR8269812.sorted.bam

#rtt109
#${OUTDIR}/SRR10916326.sorted.bam ${OUTDIR}/SRR10916324.sorted.bam ${OUTDIR}/SRR10916325.sorted.bam
