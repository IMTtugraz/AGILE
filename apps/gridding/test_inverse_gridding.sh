#!/bin/bash

#==============================================================================
# INVERSE GRIDDING TESTS
#   you can choose between inverse gridding based on:
#      - a matrix-vector multiplication or
#      - a GPU texture interpolation
#    => choose in src/config.hpp
#==============================================================================

EXECUTABLE="./gridding_backward.exe"

# 256x256 => 300x256
MATRIX="data/csrmatrix_A_76800x65536.bin"
POSITIONS="data/vector_kspace_positions_76800(300x256).bin"
INPUT="data/vector_gspace_data_65536(256x256).bin"
OUTPUT="inverse_gridding_result_vector_300x256.bin"
echo "inverse gridding: 256x256 => 300x256"
$EXECUTABLE $MATRIX $POSITIONS $INPUT 256 256 $OUTPUT

echo "#-------------------------------------------------------------------";

# 512x512 => 512x512
MATRIX="data/csrmatrix_A_262144x262144.bin"
POSITIONS="data/vector_kspace_positions_262144(512x512).bin"
INPUT="data/vector_gspace_data_262144(512x512).bin"
OUTPUT="inverse_gridding_result_vector_512x512.bin"
echo "inverse gridding: 512x512 => 512x512"
$EXECUTABLE $MATRIX $POSITIONS $INPUT 512 512 $OUTPUT

echo "#-------------------------------------------------------------------";

# 1024x1024 => 1024x1024
MATRIX="data/csrmatrix_A_1048576x1048576.bin"
POSITIONS="data/vector_kspace_positions_1048576(1024x1024).bin"
INPUT="data/vector_gspace_data_1048576(1024x1024).bin"
OUTPUT="inverse_gridding_result_vector_1024x1024.bin"
echo "inverse gridding: 1024x1024 => 1024x1024"
$EXECUTABLE $MATRIX $POSITIONS $INPUT 1024 1024 $OUTPUT

