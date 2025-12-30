#!/bin/bash
#SBATCH --job-name=homework_4		                        # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=/work/gene8940/evt82290/log.%j			    # Standard output and error log - # replace cbergman with your myid
#SBATCH --mail-user=evt82290@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/scratch/evt82290/Run124"
WORKDIR="~/NcEddieTorresFastQ"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#loading modules
module load SAMtools/1.10-GCC-8.3.0 SRA-Toolkit/2.9.6-1-centos_linux64 BWA/0.7.17-GCC-8.3.0 BCFtools/1.10.2-GCC-8.3.0

#gunzip fastq data
