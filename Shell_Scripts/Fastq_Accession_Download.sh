#!/bin/bash
#SBATCH --job-name=fqdump_gzip
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/evt82290/SRA/FastqFiles/logs/fqdump.%j.out
#SBATCH --error=/scratch/evt82290/SRA/FastqFiles/logs/fqdump.%j.err

set -euo pipefail

# Inputs/outputs
ACC_FILE="/scratch/evt82290/FastqAccessions/JGI_RNA.txt"            # one SRR per line
OUTDIR="/scratch/evt82290/SRA/FastqFiles"           # where FASTQs go
CACHEDIR="$OUTDIR/sra_cache"                        # per-job SRA cache

# Prep
mkdir -p "$OUTDIR" "$CACHEDIR" "$OUTDIR/logs"
module load SRA-Toolkit/

echo "Start: $(date)"
echo "Accession list: $ACC_FILE"
echo "OUTDIR: $OUTDIR"
echo "CACHE:  $CACHEDIR"

# Iterate accessions; skip blank lines and comments
while IFS= read -r acc; do
  [[ -z "$acc" || "$acc" =~ ^# ]] && continue
  echo "==> $acc"

  # 1) Prefetch to cache
  if ! prefetch -O "$CACHEDIR" "$acc"; then
    echo "[WARN] prefetch failed: $acc"
    continue
  fi

  # 2) Dump FASTQ to OUTDIR (use wildcard to match .sralite or directory)
  if fasterq-dump --split-files -e "$SLURM_CPUS_PER_TASK" -O "$OUTDIR" "$CACHEDIR/${acc}"*; then
    # 3) Compress
    if command -v pigz >/dev/null 2>&1; then
      pigz -p "$SLURM_CPUS_PER_TASK" -f "$OUTDIR/${acc}"_*.fastq
    else
      gzip -f "$OUTDIR/${acc}"_*.fastq
    fi
    echo "<== done: $acc"
  else
    echo "[WARN] fasterq-dump failed: $acc"
  fi
done < "$ACC_FILE"

echo "All done: $(date)"
