#!/bin/bash

#==============================================================================
# GRIDDING (lsqr) TESTS
#==============================================================================

EXECUTABLE="./gridding_forward.exe"

# 256x300 => 256x256
MATRIX="data/csrmatrix_A_76800x65536.bin"
ADJOINT="data/csrmatrix_AT_65536x76800.bin"
INPUT="data/vector_kspace_data_76800(300x256).bin"
OUTPUT="gridding_result_vector_256x256.bin"
echo "gridding: 300x256 => 256x256"
$EXECUTABLE $MATRIX $ADJOINT $INPUT $OUTPUT

echo "#-------------------------------------------------------------------";

# 512x512 => 512x512
MATRIX="data/csrmatrix_A_262144x262144.bin"
ADJOINT="data/csrmatrix_AT_262144x262144.bin"
INPUT="data/vector_kspace_data_262144(512x512).bin"
OUTPUT="gridding_result_vector_512x512.bin"
echo "gridding: 512x512 => 512x512"
$EXECUTABLE $MATRIX $ADJOINT $INPUT $OUTPUT

echo "#-------------------------------------------------------------------";

# 1024x1024 => 1024x1024
MATRIX="data/csrmatrix_A_1048576x1048576.bin"
ADJOINT="data/csrmatrix_AT_1048576x1048576.bin"
INPUT="data/vector_kspace_data_1048576(1024x1024).bin"
OUTPUT="gridding_result_vector_1024x1024.bin"
echo "gridding: 1024x1024 => 1024x1024"
$EXECUTABLE $MATRIX $ADJOINT $INPUT $OUTPUT

