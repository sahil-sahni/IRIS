#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=80g
#SBATCH --gres=lscratch:20
#SBATCH --partition=norm
#SBATCH --mail-type=BEGIN,END,TIME_LIMIT_50,TIME_LIMIT_80,FAIL
#SBATCH --cpus-per-task=20

module load cytospace

# List of spatial files to process
spatial_files=(
    "ECM01.1_spGEP.txt" "ECM01.2_spGEP.txt"
    "ECM06.1_spGEP.txt"
    "ECM08.1_spGEP.txt"
    "ECM10.1_spGEP.txt"
    "MBM05.1_spGEP.txt" "MBM05.2_spGEP.txt" "MBM05.3_spGEP.txt"
    "MBM06.1_spGEP.txt"
    "MBM07.1_spGEP.txt"
    "MBM08.1_spGEP.txt"
    "MBM11.1_spGEP.txt" "MBM11.2_spGEP.txt" "MBM11.3_spGEP.txt"
    "MBM13.1_spGEP.txt"
    "MBM18.1_spGEP.txt"
)

# Iterate over each spatial file
for spatial_file in "${spatial_files[@]}"; do
    # Extract the prefix (e.g., "MBM18.1" from "MBM18.1_spGEP.txt")
    prefix="${spatial_file%_spGEP.txt}"
    
    # Print the prefix before running the cytospace command
    echo "Processing prefix: $prefix"

    # Print singlet first
    echo "first cell type"

    cytospace \
        -sp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_scGEP.txt" \
        -ctp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_scLabel.txt" \
        -stp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_spGEP_firstct.txt" \
        -cp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_spCoord_firstct.txt" \
        -o "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_output/Izar_2022/cell_assignment_v2" \
        -op "${prefix}_firstct" \
        -sc -noss 10000 -nop 18 \
        -stctp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_spCTloc_firstct.txt"

    # Print doublet next (if applicable)
    echo "second cell type"

    cytospace \
        -sp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_scGEP.txt" \
        -ctp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_scLabel.txt" \
        -stp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_spGEP_secondct.txt" \
        -cp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_spCoord_secondct.txt" \
        -o "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_output/Izar_2022/cell_assignment_v2" \
        -op "${prefix}_secondct" \
        -sc -noss 10000 -nop 18 \
        -stctp "/data/sahnis2/final_antipdl1_project/cytospace/cytospace_input/spatial reference_izar_2022_v2/${prefix}_spCTloc_secondct.txt"

    

done
