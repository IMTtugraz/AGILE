#!/bin/bash

# FORWARD gridding means: kspace => gspace
# BACKWARD gridding means: gspace => kspace

EXEC_GRIDDING_FORWARD="./gridding_forward.exe"
EXEC_GRIDDING_BACKWARD="./gridding_backward.exe"

GSPACE_WIDTH="2048"
GSPACE_HEIGHT="2048"

DATA_TRANSFORMATION_MATRIX="data/csrmatrix_A_4194304x4194304.bin"
DATA_TRANSFORMATION_MATRIX_ADJOINT="data/csrmatrix_AT_4194304x4194304.bin"
DATA_KSPACE_INPUT_VECTOR="data/vector_kspace_data_4194304(2048x2048).bin"
DATA_KSPACE_POSITION_VECTOR="data/vector_kspace_positions_4194304(2048x2048).bin"

DATA_KSPACE_TEMP_RESULT_VECTOR="res_vector_kspace_data_2048.bin"
DATA_GSPACE_RESULT_VECTOR="res_vector_gspace_data_2048.bin"

GRIDDING_LOOPS=0

# remove existing *.bin files
#rm -f *.bin

# execution part from here
$EXEC_GRIDDING_FORWARD $DATA_TRANSFORMATION_MATRIX $DATA_TRANSFORMATION_MATRIX_ADJOINT $DATA_KSPACE_INPUT_VECTOR $DATA_GSPACE_RESULT_VECTOR
for (( i = 1;  i <= GRIDDING_LOOPS;  i++ ))
do
  echo "---------------------------- loop $i -------------------------------"
  $EXEC_GRIDDING_BACKWARD $DATA_TRANSFORMATION_MATRIX $DATA_KSPACE_POSITION_VECTOR $DATA_GSPACE_RESULT_VECTOR $GSPACE_WIDTH $GSPACE_HEIGHT $DATA_KSPACE_TEMP_RESULT_VECTOR
  $EXEC_GRIDDING_FORWARD $DATA_TRANSFORMATION_MATRIX $DATA_TRANSFORMATION_MATRIX_ADJOINT $DATA_KSPACE_TEMP_RESULT_VECTOR $DATA_GSPACE_RESULT_VECTOR
done

echo "----------------------------------------------------------------------"
echo "visualize with:"
echo "figure, imshow(abs(ifft2c(reshape(readbin_vector('$DATA_GSPACE_RESULT_VECTOR'), $GSPACE_WIDTH, $GSPACE_HEIGHT))),[]);"
