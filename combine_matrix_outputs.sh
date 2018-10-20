#!/usr/bin/bash

OUT_DIR="/home/xies/data/crispri_hamming/outputs"

# Generate numerically sorted filelist
filelist=$(ls ${OUT_DIR}/*.csv | sort --version-sort)
echo ${filelist}

cat ${filelist} >> ${OUT_DIR}/big_distance_matrix.csv
