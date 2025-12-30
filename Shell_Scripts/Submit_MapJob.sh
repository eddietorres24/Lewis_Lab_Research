#!/bin/bash
#SBATCH --job-name=RNAseq_Submit
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --time=8:00:00
#SBATCH --output=../MappingSubmit.%j.out
#SBATCH --error=../MappingSubmit.%j.err

#check for required command line argument

if [[ -z "$1" ]]
then
		echo "
		################
		ERROR: You must include the path to your accession file in the command line call.
		eg. sh submit.sh Path/To/Your/AccessionFile.txt
		##################
		"

		exit
fi

#iterates through list of accessions and passes to mapping script
#fastq directory generate by https://github.com/UGALewisLab/downloadSRA.git
# fastqPath="/scratch/evt82290/FastqFiles/2025_Run149_ET/RNA"
# outdir="/scratch/evt82290/MappingOutputs/Run149/RNA"

fastqPath="/scratch/evt82290/SRA/FastqFiles"
# fastqPath="/scratch/evt82290/FastqFiles/misc_RNA"
# outdir="/scratch/evt82290/RNAseq/cac_aberrant_transcripts"
outdir="/scratch/evt82290/MappingOutputs/RNAseq"
# scratch/evt82290/FastqFiles/2025_Run150_ET
# /scratch/evt82290/FastqFiles/2025_Run149_ET/RNA

mkdir ${outdir}
mkdir ${outdir}/logs

#make output file folders
mkdir "${outdir}/TrimmedFastQs"
mkdir "${outdir}/bamFiles"
mkdir "${outdir}/counts"
mkdir "${outdir}/bigWig"
mkdir "${outdir}/bedGraph"
mkdir "${outdir}/beds"

while read -r line

	do
	sleep 1
	echo "$line mapping job submitted"
	sbatch --export=ALL,accession="${line}",fastqPath="${fastqPath}",outdir="${outdir}" MapRNAseq.sh & done <"$1"
