#!/usr/bin/env bash
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

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

############################
# Required env from submitter
############################
: "${accession:?Missing accession (exported by submit script)}"
: "${fastqPath:?Missing fastqPath (exported by submit script)}"
: "${outdir:?Missing outdir (exported by submit script)}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

############################
# Reference / tools
############################
STAR_INDEX="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR"
GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf"
CHROMSIZES="/home/ad45368/chrom_sizes.txt"

# The single flipped GTF you created manually in $outdir:
FLIPGTF_GLOBAL="${outdir}/$(basename "$GTF" .gtf).flipped.gtf"

# Library strandedness of your RNA-seq (NEB dUTP = reverse)
LIB_STRAND="reverse"

# Modules
module load Trim_Galore
module load STAR
module load SAMtools
module load Subread
module load BEDTools
module load deepTools
module load ucsc   # bigWigToBedGraph / bedGraphToBigWig

############################
# Output dirs
############################
TRIMDIR="${outdir}/TrimmedFastQs/${accession}"
BAMDIR="${outdir}/bamFiles/${accession}"
COUNTDIR="${outdir}/counts/${accession}"
BWDIR="${outdir}/bigWig/${accession}"
BEDGRAPHDIR="${outdir}/bedGraph/${accession}"
BEDDIR="${outdir}/beds/${accession}"
TMPDIR="${outdir}/tmp/${accession}"

mkdir -p "$TRIMDIR" "$BAMDIR" "$COUNTDIR" "$BWDIR" "$BEDGRAPHDIR" "$BEDDIR" "$TMPDIR"

# STAR tmp path (let STAR create it; do NOT pre-create)
OUTTMP="${TMPDIR}/STAR_${accession}_${SLURM_JOB_ID:-$$}"
[ -d "$OUTTMP" ] && rm -rf "$OUTTMP"
trap 'rm -rf "$OUTTMP" 2>/dev/null || true' EXIT

############################
# FASTQ discovery (supports multiple dirs)
############################
echo "[DEBUG] accession=$accession"
echo "[DEBUG] fastqPath=$fastqPath"

normalize_path_list() { echo "$1" | tr ',:' ' ' | xargs -n999 echo; }

collect_search_dirs() {
  local in=($(normalize_path_list "$fastqPath"))
  local out=() d alt
  for d in "${in[@]}"; do
    out+=("$d")
    if [[ "$d" == /scratch/* ]]; then
      alt="/lustre2${d}"
      [[ -d "$alt" ]] && out+=("$alt")
    fi
  done
  printf '%s\n' "${out[@]}"
}

basename_has_accession_prefix () {
  local f="$1" base
  base="$(basename "$f")"
  [[ "$base" == "${accession}.fastq.gz" || "$base" == ${accession}_* ]]
}

find_fastqs() {
  local paths=($(collect_search_dirs))
  local d r1 r2 ru
  shopt -s nullglob
  for d in "${paths[@]}"; do
    echo "[DEBUG] scanning dir: $d"
    ls -1 "${d}/${accession}"* 2>/dev/null | head -n 6 | sed 's/^/[DEBUG] cand: /' || true

    # 1) SRA paired
    r1="${d}/${accession}_1.fastq.gz"
    r2="${d}/${accession}_2.fastq.gz"
    if [[ -f "$r1" && -f "$r2" ]] && basename_has_accession_prefix "$r1" && basename_has_accession_prefix "$r2" ; then
      R1="$r1"; R2="$r2"; MODE="PE"
      echo "[INFO] PE (SRA)"; echo "[INFO] R1=$R1"; echo "[INFO] R2=$R2"
      shopt -u nullglob; return 0
    fi

    # 2) Illumina paired (lane-aware)
    mapfile -t candR1 < <(ls -1 "${d}/${accession}"*R1*_001.fastq.gz 2>/dev/null || true)
    mapfile -t candR2 < <(ls -1 "${d}/${accession}"*R2*_001.fastq.gz 2>/dev/null || true)
    candR1=($(for f in "${candR1[@]:-}"; do basename_has_accession_prefix "$f" && echo "$f"; done))
    candR2=($(for f in "${candR2[@]:-}"; do basename_has_accession_prefix "$f" && echo "$f"; done))
    if [[ ${#candR1[@]} -gt 0 && ${#candR2[@]} -gt 0 ]]; then
      R1="${candR1[0]}"; R2="${candR2[0]}"; MODE="PE"
      echo "[INFO] PE (Illumina)"; echo "[INFO] R1=$R1"; echo "[INFO] R2=$R2"
      shopt -u nullglob; return 0
    fi

    # 3) Single-end exact
    ru="${d}/${accession}.fastq.gz"
    if [[ -f "$ru" ]] && basename_has_accession_prefix "$ru"; then
      RU="$ru"; MODE="SE"
      echo "[INFO] SE"; echo "[INFO] RU=$RU"
      shopt -u nullglob; return 0
    fi

    # 4) SE via R1-only Illumina
    mapfile -t candR1only < <(ls -1 "${d}/${accession}"*R1*_001.fastq.gz 2>/dev/null || true)
    candR1only=($(for f in "${candR1only[@]:-}"; do basename_has_accession_prefix "$f" && echo "$f"; done))
    if [[ ${#candR1only[@]} -gt 0 ]]; then
      R1="${candR1only[0]}"; MODE="SE"
      echo "[INFO] SE (R1-only)"; echo "[INFO] R1=$R1"
      shopt -u nullglob; return 0
    fi
  done
  shopt -u nullglob
  return 1
}

# Init + discover
R1=""; R2=""; RU=""; MODE=""
if ! find_fastqs; then
  echo "[ERROR] Could not locate FASTQs for ${accession} under: ${fastqPath}" >&2
  echo "[HINT] If files are under /lustre2/scratch/..., either include that in fastqPath or rely on the auto-mirror."
  exit 1
fi
if [[ "$MODE" == "PE" && "$R1" == "$R2" ]]; then
  echo "[ERROR] R1 and R2 resolved to the same file: $R1" >&2
  exit 1
fi

############################
# Trimming
############################
if [[ "$MODE" == "PE" ]]; then
  echo "[INFO] Trimming PE (SKIPPED)"
  trim_galore --illumina --paired --length 25 --basename "${accession}" --gzip \
              -o "$TRIMDIR" "$R1" "$R2"
  R1T="${TRIMDIR}/${accession}_val_1.fq.gz"
  R2T="${TRIMDIR}/${accession}_val_2.fq.gz"
else
  echo "[INFO] Trimming SE (SKIPPED)"
  SE_IN="${RU:-$R1}"
  trim_galore --illumina --length 25 --basename "${accession}" --gzip \
              -o "$TRIMDIR" "$SE_IN"
  RUT="${TRIMDIR}/${accession}_trimmed.fq.gz"
fi

############################
# STAR alignment
############################
OUTPFX="${BAMDIR}/${accession}_"
STAR_COMMON=(
  --runThreadN "$THREADS"
  --genomeDir "$STAR_INDEX"
  --outFileNamePrefix "$OUTPFX"
  --readFilesCommand zcat
  --outSAMtype BAM SortedByCoordinate
  --twopassMode Basic
  --outFilterType BySJout
  --alignIntronMax 10000
  --outSAMunmapped Within
  --outSAMattributes Standard
  --outSAMstrandField intronMotif
  --quantMode GeneCounts
  --limitBAMsortRAM 16000000000
  --outTmpDir "$OUTTMP"
)

if [[ "$MODE" == "PE" ]]; then
  STAR "${STAR_COMMON[@]}" --readFilesIn "$R1T" "$R2T"
else
  STAR "${STAR_COMMON[@]}" --readFilesIn "$RUT"
fi

BAM="${OUTPFX}Aligned.sortedByCoord.out.bam"
samtools index "$BAM"

############################
# Stranded bigWigs (CPM + log2(CPM+1))
############################
if [[ "$LIB_STRAND" == "reverse" ]]; then
  FSTRAND="reverse"; RSTRAND="forward"
else
  FSTRAND="forward"; RSTRAND="reverse"
fi

BIN=10

# Linear CPM tracks
bamCoverage -b "$BAM" -o "${BWDIR}/${accession}.${FSTRAND}.cpm.bw" \
  --normalizeUsing CPM --binSize $BIN --filterRNAstrand "$FSTRAND" -p "$THREADS" --skipNonCoveredRegions
bamCoverage -b "$BAM" -o "${BWDIR}/${accession}.${RSTRAND}.cpm.bw" \
  --normalizeUsing CPM --binSize $BIN --filterRNAstrand "$RSTRAND" -p "$THREADS" --skipNonCoveredRegions

# Helper: log2(CPM+1) via bedGraph round-trip
log_transform_bw () {
  local inbw="$1" outbw="$2" base tmpbg tmpbg2
  base="$(basename "$inbw" .bw)"
  tmpbg="${TMPDIR}/${base}.bedGraph"
  tmpbg2="${TMPDIR}/${base}.log2p1.bedGraph"
  bigWigToBedGraph "$inbw" "$tmpbg"
  awk 'BEGIN{OFS="\t"} NF>=4{v=$4+0; if(v<0)v=0; print $1,$2,$3,(log(v+1)/log(2)); next}{print}' "$tmpbg" > "$tmpbg2"
  sort -k1,1 -k2,2n "$tmpbg2" > "${tmpbg2}.sorted"
  bedGraphToBigWig "${tmpbg2}.sorted" "$CHROMSIZES" "$outbw"
  cp "$tmpbg2" "${BEDGRAPHDIR}/${base}.log2p1.bedGraph" || true
  rm -f "$tmpbg" "$tmpbg2" "${tmpbg2}.sorted"
}

# Log tracks for forward/reverse
log_transform_bw "${BWDIR}/${accession}.${FSTRAND}.cpm.bw" "${BWDIR}/${accession}.${FSTRAND}.log2p1.bw"
log_transform_bw "${BWDIR}/${accession}.${RSTRAND}.cpm.bw" "${BWDIR}/${accession}.${RSTRAND}.log2p1.bw"

############################
# Antisense-only BAM + bigWigs  (paired-end correct)
############################

# Split exons by strand (no edits to col 9)
EXONS_PLUS="${TMPDIR}/exons.plus.gtf"
EXONS_MINUS="${TMPDIR}/exons.minus.gtf"
awk '$3=="exon" && $7=="+"' "$GTF" > "$EXONS_PLUS"
awk '$3=="exon" && $7=="-"' "$GTF" > "$EXONS_MINUS"

# Detect if BAM is paired-end
IS_PAIRED=$(samtools view -c -f 1 "$BAM" || echo 0)

if [[ "$IS_PAIRED" -gt 0 ]]; then
  # ---- Paired-end, reverse-stranded (dUTP) assumed unless you set LIB_STRAND="forward" ----
  # For reverse libraries:
  #   antisense to + genes: R1 '+'  OR R2 '-'
  #   antisense to - genes: R1 '-'  OR R2 '+'
  # For forward libraries, swap the +/- in the four filters below.

  if [[ "$LIB_STRAND" == "reverse" ]]; then
    # Intersect once per strand
    bedtools intersect -wa -abam "$BAM" -b "$EXONS_PLUS"  > "${TMPDIR}/${accession}.plus.all.bam"
    bedtools intersect -wa -abam "$BAM" -b "$EXONS_MINUS" > "${TMPDIR}/${accession}.minus.all.bam"

    # + genes: R1 forward (same as '+')  OR  R2 reverse (opposite of '+')
    samtools view -b -f 64 -F 16  "${TMPDIR}/${accession}.plus.all.bam"  > "${TMPDIR}/${accession}.plus.R1_same.bam"
    samtools view -b -f 128 -f 16 "${TMPDIR}/${accession}.plus.all.bam"  > "${TMPDIR}/${accession}.plus.R2_opp.bam"

    # - genes: R1 reverse (same as '-')  OR  R2 forward (opposite of '-')
    samtools view -b -f 64  -f 16 "${TMPDIR}/${accession}.minus.all.bam" > "${TMPDIR}/${accession}.minus.R1_same.bam"
    samtools view -b -f 128 -F 16 "${TMPDIR}/${accession}.minus.all.bam" > "${TMPDIR}/${accession}.minus.R2_opp.bam"

  else
    # ---- Forward-stranded library: invert the four selections ----
    bedtools intersect -wa -abam "$BAM" -b "$EXONS_PLUS"  > "${TMPDIR}/${accession}.plus.all.bam"
    bedtools intersect -wa -abam "$BAM" -b "$EXONS_MINUS" > "${TMPDIR}/${accession}.minus.all.bam"

    # + genes: R1 reverse (antisense) OR R2 forward (antisense)
    samtools view -b -f 64 -f 16  "${TMPDIR}/${accession}.plus.all.bam"  > "${TMPDIR}/${accession}.plus.R1_same.bam"
    samtools view -b -f 128 -F 16 "${TMPDIR}/${accession}.plus.all.bam"  > "${TMPDIR}/${accession}.plus.R2_opp.bam"

    # - genes: R1 forward (antisense) OR R2 reverse (antisense)
    samtools view -b -f 64  -F 16 "${TMPDIR}/${accession}.minus.all.bam" > "${TMPDIR}/${accession}.minus.R1_same.bam"
    samtools view -b -f 128 -f 16 "${TMPDIR}/${accession}.minus.all.bam" > "${TMPDIR}/${accession}.minus.R2_opp.bam"
  fi

  # Merge antisense pieces
  samtools merge -f "$BAMDIR/${accession}_antisense.unsorted.bam" \
    "${TMPDIR}/${accession}.plus.R1_same.bam" \
    "${TMPDIR}/${accession}.plus.R2_opp.bam" \
    "${TMPDIR}/${accession}.minus.R1_same.bam" \
    "${TMPDIR}/${accession}.minus.R2_opp.bam"

  # Cleanup chunk BAMs
  rm -f "${TMPDIR}/${accession}.plus."*.bam "${TMPDIR}/${accession}.minus."*.bam \
        "${TMPDIR}/${accession}.plus.all.bam" "${TMPDIR}/${accession}.minus.all.bam"

else
  # ---- Single-end fallback ----
  # reverse-stranded: antisense == same strand as gene (-s)
  # forward-stranded: antisense == opposite strand (-S)
  ANTISENSE_FLAG="-s"; [[ "$LIB_STRAND" == "forward" ]] && ANTISENSE_FLAG="-S"

  bedtools intersect $ANTISENSE_FLAG -abam "$BAM" -b "$GTF" \
    > "$BAMDIR/${accession}_antisense.unsorted.bam"
fi

# Sort/index final antisense BAM
samtools sort -o "$BAMDIR/${accession}_antisense.bam" "$BAMDIR/${accession}_antisense.unsorted.bam"
rm -f "$BAMDIR/${accession}_antisense.unsorted.bam"
samtools index "$BAMDIR/${accession}_antisense.bam"

# Create CPM & log2(CPM+1) bigWigs if there are mapped reads
ANTI_BAM="$BAMDIR/${accession}_antisense.bam"
ANTI_MAPPED=$(samtools view -c -F 4 "$ANTI_BAM" || echo 0)
if [[ "$ANTI_MAPPED" -eq 0 ]]; then
  echo "[WARN] $ANTI_BAM has 0 mapped reads; skipping antisense bigWig."
else
  bamCoverage -b "$ANTI_BAM" -o "${BWDIR}/${accession}.antisense.cpm.bw" \
    --normalizeUsing CPM --binSize $BIN -p "$THREADS" --skipNonCoveredRegions
  log_transform_bw "${BWDIR}/${accession}.antisense.cpm.bw" \
                   "${BWDIR}/${accession}.antisense.log2p1.bw"
fi

############################
# 5′-end CPM bigWigs (forward/reverse)
############################
bamCoverage -b "$BAM" \
  -o "${BWDIR}/${accession}.forward.5prime.cpm.bw" \
  --filterRNAstrand forward --Offset 1 --binSize 1 \
  --normalizeUsing CPM -p "$THREADS" --skipNonCoveredRegions

bamCoverage -b "$BAM" \
  -o "${BWDIR}/${accession}.reverse.5prime.cpm.bw" \
  --filterRNAstrand reverse --Offset 1 --binSize 1 \
  --normalizeUsing CPM -p "$THREADS" --skipNonCoveredRegions

############################
# featureCounts: sense vs antisense (use single flipped GTF in $outdir)
############################
# Paired-end detection
IS_PAIRED=$(samtools view -c -f 1 "$BAM" || echo 0)
if [[ "$IS_PAIRED" -gt 0 ]]; then
  FC_PE_ARGS="-p -B -C"
  echo "[INFO] featureCounts: paired-end mode"
else
  FC_PE_ARGS=""
  echo "[INFO] featureCounts: single-end mode"
fi

# Strandedness for featureCounts (reverse = -s 2)
FC_STRAND_ARG="-s 2"; [[ "$LIB_STRAND" == "forward" ]] && FC_STRAND_ARG="-s 1"

# Sense counts (original GTF)
featureCounts -T "$THREADS" $FC_STRAND_ARG $FC_PE_ARGS -t exon -g gene_name \
  -a "$GTF" -o "${COUNTDIR}/${accession}.sense.txt" "$BAM"

# Antisense counts (single flipped GTF in outdir)
if [[ ! -s "$FLIPGTF_GLOBAL" ]]; then
  echo "[FATAL] Flipped GTF not found: $FLIPGTF_GLOBAL" >&2
  exit 1
fi
featureCounts -T "$THREADS" $FC_STRAND_ARG $FC_PE_ARGS -t exon -g gene_name \
  -a "$FLIPGTF_GLOBAL" -o "${COUNTDIR}/${accession}.antisense.txt" "$BAM"

echo "[DONE] ${accession}"

############################
# Optional per-sample TMP cleanup
############################
rm -rf "${TMPDIR}" 2>/dev/null || true
