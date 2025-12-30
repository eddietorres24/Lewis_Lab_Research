#!/bin/bash
#SBATCH --job-name=zl_heatmap
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=20gb
#SBATCH --time=1:00:00
#SBATCH --output=../heatmap_edge.%j.out
#SBATCH --error=../heatmap_edge.%j.err

OUTDIR2="/scratch/evt82290/Run136/BigWigs"

module load deepTools/3.5.2-foss-2022a

#Peakfile for H3K9me3 is based on CutNRun Sample 114-19

####Plotting scaled enrichment across all het domains
computeMatrix scale-regions -p 12 -b 1000 -a 1000 \
-S  $OUTDIR2/6147_136-1_ChIP_WT_H3K27me3_abcam_Rep2.bin_25.smooth_75Bulk.bw \
    $OUTDIR2/6147_136-2_ChIP_cac-1_H3K27me3_abcam_Rep2.bin_25.smooth_75Bulk.bw \
    $OUTDIR2/6147_136-3_ChIP_cac-2_H3K27me3_abcam_Rep2_S3_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
    $OUTDIR2/6147_136-4_ChIP_cac-3_H3K27me3_abcam_Rep2_S4_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
    $OUTDIR2/6147_136-83_ChIP_set-7_H3K27me3_CS_Rep2_S82_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
  -R Figure2_K27regions_Scaledcenter_FileToCheckOrderFINAL_ZL.txt -o edge_map.matrix \
	--outFileNameMatrix edge_map.txt \
	--sortRegions keep \
   --missingDataAsZero -bs 10
 # --sortUsing sum --sortUsingSamples 5 2
	#--sortRegions descend --sortUsing sum \
 	 #  --sortUsingSamples 3


#https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html
plotHeatmap -m edge_map.matrix \
	-o edge_map_K27.pdf \
		--sortRegions keep  \
			--outFileSortedRegions heatmap.bed \
			--startLabel "5'" \
				--endLabel "3'" \
			--samplesLabel "WT K27me3" "cac-1 K27me3" "cac-2 K27me3" "cac-3 K27me3" "set-7 K27me3" \
					--heatmapHeight 9 \
						--heatmapWidth 4 \
					--zMax 20 20 20 20 20 \
						--colorMap 'Blues'

#Peak file for K27me3 is based on CutNRun Sample 119-42
###Plotting scaled enrichment across all K27 domains
# computeMatrix scale-regions -p 12 -R /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/Peaks/Figure2_K27PeakFile.txt -S \
#  	/scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/119-49_CUT_RUN_WT_H3K9__S37_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
# 	/scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/133-75_CUTnRUN_WT_H3K27m2m3__S72_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
# 	/scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/119-41_CUT_RUN_WT_K27me3_.bin_25.smooth_75Bulk.bw \
# 	/scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/119-50_CUT_RUN_hda-1_H3K9__S38_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
# 	/scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/6147_136-18_CUTnRUN_HDA-1_H3K27me2m3.bin_25.smooth_75Bulk.bw \
# 	/scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/119-42_CUT_RUN_hda-1_K27me3__S30_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
# 	-b 1000 -a 1000 \
# 	-o /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/MatrixFiles/Figure2_K27regions_Scaledcenter.matrix \
# 	--outFileNameMatrix Figure2_K27regions_ScaledCenterFinal.matrix.txt \
# 	--sortRegions keep \
#    --missingDataAsZero -bs 10
#  # --sortUsing sum --sortUsingSamples 5 2
# 	#--sortRegions descend --sortUsing sum \
#  	 #  --sortUsingSamples 3


#https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html
# plotHeatmap -m /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/MatrixFiles/Figure2_K27regions_Scaledcenter.matrix \
# 	-o /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/HeatmapFiles/Figure2_K27peaks_ScaledCenterFINAL.pdf \
# 			--sortRegions keep \
# 			--outFileSortedRegions /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/Peaks/Figure2_K27regions_Scaledcenter_FileToCheckOrderFINAL.txt \
# 			--startLabel "5'" \
# 				--endLabel "3'" \
# 			--samplesLabel "WT K9" "WT K27me2/3" "WT K27me3" "hda-1 K9" "hda-1 K27me2/3"  "hda-1 K27me3" \
# 					--heatmapHeight 9 \
# 						--heatmapWidth 4 \
# 					--zMax 10 10 12 10 10 12 \
# 						--colorMap 'Blues' 'Greens' 'Greens' 'Blues' 'Greens' 'Greens'
