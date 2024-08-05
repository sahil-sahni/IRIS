#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=80g
#SBATCH --gres=lscratch:20
#SBATCH --partition=norm
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT_50,TIME_LIMIT_80,FAIL
#SBATCH --cpus-per-task=10

module load cytospace

# List of spatial files to process
spatial_files=("mel1_rep1_spGEP.txt" "mel1_rep2_spGEP.txt" "mel2_rep1_spGEP.txt" "mel2_rep2_spGEP.txt" "mel3_rep1_spGEP.txt" "mel3_rep2_spGEP.txt" "mel4_rep1_spGEP.txt" "mel4_rep2_spGEP.txt")

# Iterate over each spatial file
for spatial_file in "${spatial_files[@]}"; do
    # Extract the prefix (e.g., "mel1_rep1" from "mel1_rep1_spGEP.txt")
    prefix=$(basename "$spatial_file" | sed 's/_spGEP.txt//')

    cytospace \
       -sp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/singlecell reference_livnat_2018/livnat18_scGEP_naive_counts.txt" \
       -ctp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/singlecell reference_livnat_2018/livnat18_scCellLabels_naive.txt" \
       -stp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_thrane_2018/$spatial_file" \
       -cp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_thrane_2018/${prefix}_spCoord.txt" \
       -o "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_output/Thrane_2018/cell_assignment" \
       -op "$prefix" \
       -g 'square' \
       -mcn 20 \
       --downsample-off \
       -sss -nosss 200 \
       -ctfep "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_output/Thrane_2018/cell_fraction/${prefix}Seurat_weights.txt"

done
