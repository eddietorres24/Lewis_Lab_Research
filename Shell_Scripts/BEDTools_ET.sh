#!/bin/bash
#SBATCH --job-name=BEDTools_ET
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../BEDTools/logs/CallPeak.%j.out
#SBATCH --error=../BEDTools/logs/CallPeak.%j.err

###THIS SCRIPT CONTAINS COMMANDS I USE TO MANIPULATE BED FILES###

#Convert Broadpeaks to bed format

for file in *.broadPeak; do cut -f1-6 "$file" > "${file%.broadPeak}_peaks.bed"; done


#bedtools
# module load BEDTools

# bedtools intersect -a -b -wa #overlap b/w a and b, keep sequence in file a
# bedtools intersect -a -b -v #no overlap b/w a and b, keep sequence in file a
# bedtools intersect -a qa-suz12_no_WT.bed -b WT_macs_24hr_rep2.bed -v > qa-suz12_no_WT_0_or_24.bed
# bedtools intersect -a internal_K27_cac3_peaks.bed -b WT_macs_0hr_rep2.bed -v > internal_K27_ectopic.bed
# bedtools intersect -a neurospora.bed -b WT_K27_narrow_rep2.bed -wa -f 0.9 > K27_mraked_genes.bed
# bedtools intersect -a CAF-1_ATAC_Peaks_merge_2.bed -b WT_ATAC_peaks.bed -wa > CAF1_ATAC_WT.bed
# bedtools intersect -a subtelomeric_K27_no_cac-3.bed -b WT_macs_0hr_rep2.bed -wa > subtelomeric_K27_normal.bed

# bedtools intersect -a cac-1-2_K27.bed -b WT_CS_H3K27me3_Rep1_peaks.bed -wa > WT_cac-1-2_K27.bed

#ATAC
# bedtools intersect -a CAF1_ATAC_WT_peaks.bed -b WT_abc_H3K27me3_Rep2_peaks.narrowPeak -wa > WT_ATAC_K27_regions.bed
# bedtools intersect -a CAF1_ATAC_NoWT.bed -b WT_abc_H3K27me3_Rep2_peaks.narrowPeak -wa > CAF-1_only_ATAC_K27_regions.bed

# bedtools intersect -a K27_narrow_genes_sorted.bed -b cac1_up.bed -wa > cac1_up_K27_genes.bed
# bedtools intersect -a K27_narrow_genes_sorted.bed -b cac2_up.bed -wa > cac2_up_K27_genes.bed
#
# bedtools intersect -a all_genes_gff.bed -b K27_narrow_genes_sorted.bed -wa > K27_genes_gff.bed
# bedtools intersect -a all_genes_gff.bed -b CAF-1_All_K4_merge_peaks.bed -wa > K4_genes_gff.bed

# bedtools intersect -a CAF1_ATAC_WT.bed -b cac-2_ATAC_peaks.bed -v > WT_ATAC_NoCAF2.bed
# bedtools intersect -a CAF1_ATAC_WT.bed -b cac-3_ATAC_peaks.bed -v > WT_ATAC_NoCAF3.bed
# bedtools intersect -a CAF1_ATAC_WT.bed -b set-7_ATAC_peaks.bed -v > WT_ATAC_NoCAF4.bed

#
# bedtools intersect -a CAF1_ATAC_WT.bed -b WT_abc_H3K27me3_Rep2_peaks.narrowPeak -wa > ATAC_WT_K27_regions.bed
# bedtools intersect -a cac1_2_set7_ATAC_NoWT.bed -b WT_abc_H3K27me3_Rep2_peaks.narrowPeak -wa > cac1_2_set7_ATAC_NoWT_K27_regions.bed
#
# bedtools intersect -a all_genes_gff_names.bed -b K27_genes_trimmed.bed -v > nonK27_genes.bed
#
# bedtools intersect -a K27_genes_stringent.bed -b cac1_H3K27me3_Rep2_peaks.bed -v > K27genes_not_in_cac1_new.bed
# bedtools intersect -a K27_genes_stringent.bed -b cac2_H3K27me3_Rep2_peaks.bed -v > K27genes_not_in_cac2_new.bed
# bedtools intersect -a K27_genes_stringent.bed -b cac3_H3K27me3_Rep2_peaks.bed -wa > K27genes_in_cac3_new.bed
# bedtools intersect -a K27_genes_stringent.bed -b K27genes_in_cac1_2_new.bed -v > K27genes_NOT_in_cac1_2_new.bed
# bedtools intersect -a K27genes_in_cac1_2_new.bed -b cac3_H3K27me3_Rep2_peaks.bed -v > K27genes_in_cac1_2_NOT_in_cac3_new.bed
# bedtools intersect -a K27genes_in_cac1_2_new.bed -b K27genes_in_cac1_2_NOT_in_cac3_new.bed -v > K27genes_in_cac1_2_3_new.bed
#
# bedtools intersect -a all_genes_gff_names.bed -b CAF-1_Ectopic_K27 -wa -f 0.8 > CAF-1_ectopic_K27_genes.bed
#
# bedtools intersect -a all_genes_gff_names.bed -b cac-1_ectopic_K27_peaks.bed -wa -f 0.8 > cac-1_ectopic_K27_genes.bed
# bedtools intersect -a all_genes_gff_names.bed -b cac-2_ectopic_K27_peaks.bed -wa -f 0.8 > cac-2_ectopic_K27_genes.bed
# bedtools intersect -a all_genes_gff_names.bed -b cac-3_ectopic_K27_peaks.bed -wa -f 0.8 > cac-3_ectopic_K27_genes.bed
#
# #K27 peaks
# bedtools intersect -a WT_H3K27me3_Rep2_peaks.bed -b cac2_H3K27me3_Rep2_peaks.bed -v > WT_K27_peaks_NO_cac2.bed
# bedtools intersect -a WT_H3K27me3_Rep2_peaks.bed -b cac2_H3K27me3_Rep2_peaks.bed -wa > WT_K27_peaks_IN_cac2.bed
# bedtools intersect -a cac2_H3K27me3_Rep2_peaks.bed -b WT_H3K27me3_Rep2_peaks.bed -v > cac-2_ectopic_K27_peaks.bed
# bedtools intersect -a WT_H3K27me3_Rep2_peaks.bed -b cac1_H3K27me3_Rep2_peaks.bed -v > WT_K27_peaks_NO_cac1.bed
# bedtools intersect -a WT_H3K27me3_Rep2_peaks.bed -b cac1_H3K27me3_Rep2_peaks.bed -wa > WT_K27_peaks_IN_cac1.bed
# bedtools intersect -a cac1_H3K27me3_Rep2_peaks.bed -b WT_H3K27me3_Rep2_peaks.bed -v > cac-1_ectopic_K27_peaks.bed
# bedtools intersect -a WT_H3K27me3_Rep2_peaks.bed -b cac3_H3K27me3_Rep2_peaks.bed -v > WT_K27_peaks_NO_cac3.bed
# bedtools intersect -a WT_H3K27me3_Rep2_peaks.bed -b cac3_H3K27me3_Rep2_peaks.bed -wa > WT_K27_peaks_IN_cac3.bed
# bedtools intersect -a cac3_H3K27me3_Rep2_peaks.bed -b WT_H3K27me3_Rep2_peaks.bed -v > cac-3_ectopic_K27_peaks.bed

#K4 marked genes
# bedtools intersect -a all_genes_gff.bed -b WT_H3K4me2_Rep1_peaks.bed -wa > WT_K4_genes.bed
# bedtools intersect -a all_genes_gff.bed -b cac-1_H3K4me2_Rep1_peaks.bed -wa > cac-1_K4_genes.bed
# bedtools intersect -a all_genes_gff.bed -b cac-2_H3K4me2_Rep1_peaks.bed -wa > cac-2_K4_genes.bed
# bedtools intersect -a all_genes_gff.bed -b cac-3_H3K4me2_Rep1_peaks.bed -wa > cac-3_K4_genes.bed
# bedtools intersect -a all_genes_gff.bed -b set-7_H3K4me2_Rep1_peaks.bed -wa > set-7_K4_genes.bed
# #retained
# bedtools intersect -a WT_K4_genes.bed -b cac-1_K4_genes.bed -wa > WT_K4_genes_IN_cac-1.bed
# bedtools intersect -a WT_K4_genes.bed -b cac-2_K4_genes.bed -wa > WT_K4_genes_IN_cac-2.bed
# bedtools intersect -a WT_K4_genes.bed -b cac-3_K4_genes.bed -wa > WT_K4_genes_IN_cac-3.bed
# bedtools intersect -a WT_K4_genes.bed -b set-7_K4_genes.bed -wa > WT_K4_genes_IN_set-7.bed
# #lost
# bedtools intersect -a WT_K4_genes.bed -b cac-1_K4_genes.bed -v > WT_K4_genes_NO_cac-1.bed
# bedtools intersect -a WT_K4_genes.bed -b cac-2_K4_genes.bed -v > WT_K4_genes_NO_cac-2.bed
# bedtools intersect -a WT_K4_genes.bed -b cac-3_K4_genes.bed -v > WT_K4_genes_NO_cac-3.bed
# bedtools intersect -a WT_K4_genes.bed -b set-7_K4_genes.bed -v > WT_K4_genes_NO_set-7.bed
# #ectopic
# bedtools intersect -a cac-1_K4_genes.bed -b WT_K4_genes.bed -v > cac-1_ectopic_K4_genes.bed
# bedtools intersect -a cac-2_K4_genes.bed -b WT_K4_genes.bed -v > cac-2_ectopic_K4_genes.bed
# bedtools intersect -a cac-3_K4_genes.bed -b WT_K4_genes.bed -v > cac-3_ectopic_K4_genes.bed
# bedtools intersect -a set-7_K4_genes.bed -b WT_K4_genes.bed -v > set-7_ectopic_K4_genes.bed
# #ectopic K4 @ K27
# bedtools intersect -a cac-1_ectopic_K4_genes.bed -b K27_genes_stringent.bed -wa > cac-1_ectopic_K4_genes_K27_regions.bed
# bedtools intersect -a cac-2_ectopic_K4_genes.bed -b K27_genes_stringent.bed -wa > cac-2_ectopic_K4_genes_K27_regions.bed
# bedtools intersect -a cac-3_ectopic_K4_genes.bed -b K27_genes_stringent.bed -wa > cac-3_ectopic_K4_genes_K27_regions.bed
# bedtools intersect -a set-7_ectopic_K4_genes.bed -b K27_genes_stringent.bed -wa > set-7_ectopic_K4_genes_K27_regions.bed

#K27 genes + & -
bedtools intersect -a plus_genes.bed -b K27_genes_stringent.bed -wa > K27_genes_plus.bed
bedtools intersect -a minus_genes.bed -b K27_genes_stringent.bed -wa > K27_genes_minus.bed

bedtools intersect -a plus_genes.bed -b minus_genes.bed -wa > plus_overlap.bed
bedtools intersect -a minus_genes.bed -b plus_genes.bed -wa > minus_overlap.bed

bedtools intersect -a plus_genes.bed -b minus_genes.bed -wa > plus_overlap.bed
bedtools intersect -a minus_genes.bed -b plus_genes.bed  -wa > minus_overlap.bed
cat plus_overlap.bed minus_overlap.bed | sort -k1,1 -k2,2n -k3,3n > overlapping_genes.bed

bedtools intersect -a nonK27_genes.bed -b H3K36me3_no_ash1_noK27.bed -v > NonK27_nonASH1_genes.bed

#Combining all overlapping peaks & merging
# CAF-1
# bedtools multiinter -header -i ${OUTDIR3}/2024_04_23_WT_peaks.bed \
#                                ${OUTDIR3}/2024_04_23_136_abcam_cac-1_peaks.bed \
#                                ${OUTDIR3}/2024_04_23_136_abcam_cac-2_peaks.bed > ${OUTDIR3}/merge_peaks.txt
#

#MERGE OVERLAPPING REGIONS WITHIN A BED FILE
bedtools sort -i WT_qa_peaks_0_24hr.bed | bedtools merge -i - > WT_qa_peaks_0_24hr_merged.bed

#test for Duplicates
# sort K27genes_in_cac1_2.bed | uniq -d > test.bed

#DEDUPLICATE FILES for whole directory
# for file in *.bed; do
#     echo "Processing $file..."
#     sort -k1,1 -k2,2n -k3,3n "$file" | uniq > "deduplicated_beds/$file"
# done

#REMOVE NON-CORE LG REGIONS AND DEDUPLICATE
for file in *.bed; do
    echo "Processing $file..."
    awk '$1 ~ /^CM/' "$file" | sort -k1,1 -k2,2n -k3,3n | uniq > "deduplicated_beds/$file"
done

### qa-suz12 ###

###rethinking ectopic peak calls
#going to keep all regions not overlapping WT samples and cut out small peaks (< 650 bp)


#Rep1
# bedtools multiinter-bams ${OUTDIR1}/WT_0hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/qa-suz12_4hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/qa-suz12_8hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/qa-suz12_12hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/qa-suz12_24hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/WT_24hr_H3K27me3_Rep1_peaks_sorted.bed > ${OUTDIR1}/qa-suz12_rep1_overlap_peaks.bed

 cat WT_macs_0hr_rep2.bed \
  WT_24hr_rep2.bed > WT_qa_peaks_0_24hr.bed
 #     qa-suz12_24hr_H3K27me3_Rep3_peaks.bed \
 #      qa-suz12_8hr_H3K27me3_Rep3_peaks.bed \
 #      qa-suz12_4hr_H3K27me3_Rep3_peaks.bed \
 #     qa-suz12_12hr_H3K27me3_Rep3_peaks.bed > qa-suz12_WT_K27.bed

###SUBETTING QA-SUZ12 BEDS
bedtools intersect -a qa-suz12_WT_K27_merged.bed -b WT_qa_peaks_0_24hr_merged.bed -v > qa-suz12_ectopic_peaks.bed

#processing csaw regions
# Normal No Recover → qa_unrecovered.bed
awk '{ if ($3 - $2 >= 301) print }' normal_no_recover_regions.bed | \
bedtools sort -i - | \
bedtools merge -i - -d 1000 | \
awk '{ if ($3 - $2 >= 1000) print }' > qa_unrecovered.bed

# Normal Recover → qa_recovered.bed
awk '{ if ($3 - $2 >= 301) print }' normal_recover_regions.bed | \
bedtools sort -i - | \
bedtools merge -i - -d 1000 | \
awk '{ if ($3 - $2 >= 1000) print }' > qa_recovered.bed

# Ectopic → qa_ectopic.bed
awk '{ if ($3 - $2 >= 301) print }' ectopic_regions.bed | \
bedtools sort -i - | \
bedtools merge -i - -d 1000 | \
awk '{ if ($3 - $2 >= 1000) print }' > qa_ectopic.bed


###splitting + and - genes
# 0) sort once (required for -sorted speedups later if you want)
bedtools sort -i all_genes_gff_names.bed > genes.sorted.bed

# 1) split by strand (BED6 assumed: col6 is strand)
awk 'BEGIN{OFS="\t"} $6=="+"' genes.sorted.bed > genes.plus.bed
awk 'BEGIN{OFS="\t"} $6=="-"' genes.sorted.bed > genes.minus.bed

# 2) drop genes whose *gene body* overlaps any opposite-strand gene body
bedtools intersect -v -a genes.minus.bed -b genes.plus.bed > minus_genes.bed
bedtools intersect -v -a genes.plus.bed  -b genes.minus.bed > plus_genes.bed

# (optional) sanity check: these should print nothing
bedtools intersect -u -a minus_genes.bed -b genes.plus.bed | head
bedtools intersect -u -a plus_genes.bed  -b genes.minus.bed | head
