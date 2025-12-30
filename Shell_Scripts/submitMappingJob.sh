#!/bin/bash

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

fastqPath="/scratch/evt82290/FastqFiles/Run145"
outdir="/scratch/evt82290/MappingOutputs/Run145"

mkdir ${outdir}
mkdir ${outdir}/logs

#make output file folders
mkdir "${outdir}/TrimmedFastQs"
mkdir "${outdir}/bamFiles"
mkdir "${outdir}/bigWig"
mkdir "${outdir}/Peak"


while read -r line

	do
	sleep 2
	echo "$line mapping job submitted"
	sbatch --job-name="${line}" --export=ALL,accession="${line}",fastqPath="${fastqPath}",outdir="${outdir}" MapData.sh & done <"$1"
