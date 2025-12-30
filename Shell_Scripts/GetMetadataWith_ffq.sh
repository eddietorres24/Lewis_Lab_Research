#!/bin/bash
#SBATCH --job-name=GetMetaDataWith_ffq.sh
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=144:00:00
#SBATCH --output=./fetchfastq.%j.out
#SBATCH --error=./fetchfastq.%j.err

cd $SLURM_SUBMIT_DIR

##This script takes an input of SRR accession #s one per line. For each accessiob. Edit last line to include path to your list of accessions.
#Your accession file should include one accession per line with no other information.

accession="/home/evt82290/Research/accession_clean_ET.txt"

module load ffq/0.3.0-foss-2022a

#make a directory to store fastq metadata [.JSON] format and files.
mkdir ../FastqFiles

#For each line of the accession file, make a directory and fetch fastq metadata using ffq

while read -r line
do
#make a new directory for each accession
  mkdir ../FastqFiles/${line}
#download metadata and store in a .JSON file
  ffq -o ../FastqFiles/${line}/${line}.json --ftp "$line"

done <"${accession}"
