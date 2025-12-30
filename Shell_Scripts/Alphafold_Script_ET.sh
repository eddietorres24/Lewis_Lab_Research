#!/bin/bash
#SBATCH --job-name=PRC2
#SBATCH --partition=gpu_p
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:P100:1
#SBATCH --mem=60gb
#SBATCH --time=120:00:00
#SBATCH --output=../Alphafold/logs/%x.out
#SBATCH --error=../Alphafold/logs/%x.err

#Output directory

OUTDIR="/scratch/evt82290/alphafold"
FASTADIR="/home/evt82290/Research/FASTA_Protein_Sequences"

#Run Alphafold
# cd $SLURM_SUBMIT_DIR
#
# ml AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0
#
# alphafold --data_dir /apps/db/AlphaFold/2.3.1 --output_dir ${OUTDIR}/PRC2 --fasta_paths ${FASTADIR}/PRC2.fasta --max_template_date 3000-01-01

#exmaple
#alphafold --data_dir /apps/db/AlphaFold/2.3.1 --output_dir ./output --model_names model_1 --fasta_paths ./query.fasta --max_template_date 2021-11-17

#Code for monomer
# cd $SLURM_SUBMIT_DIR
#
# ml purge
# export SINGULARITYENV_TF_FORCE_UNIFIED_MEMORY=1
# export SINGULARITYENV_XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
#
# export ALPHAFOLD_DATA_DIR=/apps/db/AlphaFold/2.3.1
#
# singularity exec -B /apps/db/AlphaFold -B /apps/eb/CUDAcore/11.2.1 \
# --nv /apps/singularity-images/alphafold_2.3.1_cuda112.sif python /app/alphafold/run_alphafold.py \
# --use_gpu_relax \
# --data_dir=$ALPHAFOLD_DATA_DIR \
# --uniref90_database_path=$ALPHAFOLD_DATA_DIR/uniref90/uniref90.fasta \
# --mgnify_database_path=$ALPHAFOLD_DATA_DIR/mgnify/mgy_clusters.fa \
# --bfd_database_path=$ALPHAFOLD_DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
# --uniref30_database_path=$ALPHAFOLD_DATA_DIR/uniref30/UniRef30_2021_03 \
# --pdb_seqres_database_path=$ALPHAFOLD_DATA_DIR/pdb_seqres/pdb_seqres.txt \
# --template_mmcif_dir=$ALPHAFOLD_DATA_DIR/pdb_mmcif/mmcif_files \
# --obsolete_pdbs_path=$ALPHAFOLD_DATA_DIR/pdb_mmcif/obsolete.dat \
# --uniprot_database_path=$ALPHAFOLD_DATA_DIR/uniprot/uniprot.fasta \
# --model_preset=multimer \
# --max_template_date=2022-10-01 \
# --db_preset=full_dbs \
# --output_dir=./output \
# --fasta_paths=./input.fa


#Code for multimer
cd $SLURM_SUBMIT_DIR

ml purge
export SINGULARITYENV_TF_FORCE_UNIFIED_MEMORY=1
export SINGULARITYENV_XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

export ALPHAFOLD_DATA_DIR=/apps/db/AlphaFold/2.3.1

singularity exec -B /apps/db/AlphaFold -B /apps/eb/CUDAcore/11.2.1 \
--nv /apps/singularity-images/alphafold_2.3.1_cuda112.sif python /app/alphafold/run_alphafold.py \
--use_gpu_relax \
--data_dir=$ALPHAFOLD_DATA_DIR \
--uniref90_database_path=$ALPHAFOLD_DATA_DIR/uniref90/uniref90.fasta \
--mgnify_database_path=$ALPHAFOLD_DATA_DIR/mgnify/mgy_clusters.fa \
--bfd_database_path=$ALPHAFOLD_DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
--uniref30_database_path=$ALPHAFOLD_DATA_DIR/uniref30/UniRef30_2021_03 \
--pdb_seqres_database_path=$ALPHAFOLD_DATA_DIR/pdb_seqres/pdb_seqres.txt \
--template_mmcif_dir=$ALPHAFOLD_DATA_DIR/pdb_mmcif/mmcif_files \
--obsolete_pdbs_path=$ALPHAFOLD_DATA_DIR/pdb_mmcif/obsolete.dat \
--uniprot_database_path=$ALPHAFOLD_DATA_DIR/uniprot/uniprot.fasta \
--model_preset=multimer \
--max_template_date=3000-01-01 \
--db_preset=full_dbs \
--output_dir=${OUTDIR}/PRC2 \
--fasta_paths=${FASTADIR}/PRC2.fasta
