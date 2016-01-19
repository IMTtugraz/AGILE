#!/bin/bash

# FORWARD gridding means: kspace => gspace
# BACKWARD gridding means: gspace => kspace

EXEC_GRIDDING_FORWARD="./gridding_forward.exe"
EXEC_GRIDDING_BACKWARD="./gridding_backward.exe"

GSPACE_WIDTH="512"
GSPACE_HEIGHT="512"

DATA_TRANSFORMATION_MATRIX="data/regridding_test_from_florian/csrmatrix_A_262144x262144.bin"
DATA_TRANSFORMATION_MATRIX_ADJOINT="data/regridding_test_from_florian/csrmatrix_AT_262144x262144.bin"

DATA_KSPACE_INPUT_VECTOR="data/regridding_test_from_florian/vector_kspace_data_262144(512x512).bin"
#DATA_KSPACE_INPUT_VECTOR="data/vector_kspace_data_262144(512x512).bin"
#DATA_KSPACE_INPUT_VECTOR="res_r.bin"

DATA_KSPACE_POSITION_VECTOR="data/vector_kspace_positions_262144(512x512).bin"

DATA_KSPACE_TEMP_RESULT_VECTOR="florian_res_vector_kspace_data_262144.bin"
DATA_GSPACE_RESULT_VECTOR="florian_res_vector_gspace_data_262144.bin"

GRIDDING_LOOPS=100

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
echo "figure, imshow(abs(ifft2(reshape(readbin_vector('$DATA_GSPACE_RESULT_VECTOR'), $GSPACE_WIDTH, $GSPACE_HEIGHT))),[]);"
