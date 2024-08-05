#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=80g
#SBATCH --gres=lscratch:20
#SBATCH --partition=norm
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT_50,TIME_LIMIT_80,FAIL
#SBATCH --cpus-per-task=10

module load R

# Function to get the prefix from the spatial file name
get_prefix() {
  basename "$1" | sed 's/_spGEP.txt//'
}

# List of spatial files to process
spatial_files=("mel1_rep1_spGEP.txt" "mel1_rep2_spGEP.txt" "mel2_rep1_spGEP.txt" "mel2_rep2_spGEP.txt" "mel3_rep1_spGEP.txt" "mel3_rep2_spGEP.txt" "mel4_rep1_spGEP.txt" "mel4_rep2_spGEP.txt")

# Paths
sc_path="/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/singlecell reference_livnat_2018/livnat18_scGEP_naive_counts.txt"
ct_path="/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/singlecell reference_livnat_2018/livnat18_scCellLabels_naive.txt"
outdir="/data/sahnis2/final_antipdl1_project/cytospace/cytospace_output/Thrane_2018/cell_fraction"
disable_downsampling="TRUE"

# Iterate over each spatial file
for spatial_file in "${spatial_files[@]}"; do
    # Extract the prefix (e.g., "mel1_rep1" from "mel1_rep1_spGEP.txt")
    prefix=$(get_prefix "$spatial_file")

    # set paths
    st_path="/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_thrane_2018/$spatial_file"

    # Call Rscript with arguments
    Rscript /data/sahnis2/final_antipdl1_project/cytospace/code/get_cellfracs_seuratv3_edit.R --scrna-path "$sc_path" --ct-path "$ct_path" --st-path "$st_path" --outdir "$outdir" --prefix "$prefix" --disable-fraction-downsampling
done
