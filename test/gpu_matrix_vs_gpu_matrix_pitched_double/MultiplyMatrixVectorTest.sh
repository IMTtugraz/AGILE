#!/bin/bash

#MULTIPLY MATRIX VECTOR

#ROWS = x * COLUMNS

#TEST1: MULTIPLY_MATRIX_VECTOR ||  ROWS = 2 * COLUMNS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_double
cp gpu_matrix_vs_gpu_matrix_pitched_double_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 1;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVector/R2C_Float.xls!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!num_rows = test_ctr;!num_rows = test_ctr*2;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_double/gpu_matrix_vs_gpu_matrix_pitched_double_test 
#######################################################################################################################################################
#TEST5: MULTIPLY_MATRIX_VECTOR || ROWS = 5 * COLUMNS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_double
cp gpu_matrix_vs_gpu_matrix_pitched_double_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 1;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVector/R5C_Float.xls!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!num_rows = test_ctr;!num_rows = test_ctr*5;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_double/gpu_matrix_vs_gpu_matrix_pitched_double_test 
#######################################################################################################################################################
#TEST9: MULTIPLY_MATRIX_VECTOR || ROWS = 10 * COLUMNS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_double
cp gpu_matrix_vs_gpu_matrix_pitched_double_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 1;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVector/R10C_Float.xls!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!num_rows = test_ctr;!num_rows = test_ctr*10;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_double/gpu_matrix_vs_gpu_matrix_pitched_double_test 



#COLUMNS = x * ROWS
#######################################################################################################################################################
#TEST13: MULTIPLY_MATRIX_VECTOR ||  COLUMNS = 2 * ROWS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_double
cp gpu_matrix_vs_gpu_matrix_pitched_double_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 1;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVector/C2R_Float.xls!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!num_columns = test_ctr;!num_columns = test_ctr*2;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_double/gpu_matrix_vs_gpu_matrix_pitched_double_test 
#######################################################################################################################################################
#TEST17: MULTIPLY_MATRIX_VECTOR || COLUMNS = 5 * ROWS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_double
cp gpu_matrix_vs_gpu_matrix_pitched_double_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 1;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVector/C5R_Float.xls!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!num_columns = test_ctr;!num_columns = test_ctr*5;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_double/gpu_matrix_vs_gpu_matrix_pitched_double_test 
#######################################################################################################################################################
#TEST21: MULTIPLY_MATRIX_VECTOR || COLUMNS = 10 * ROWS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_double
cp gpu_matrix_vs_gpu_matrix_pitched_double_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 1;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVector/C10R_Float.xls!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!num_columns = test_ctr;!num_columns = test_ctr*10;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_double/gpu_matrix_vs_gpu_matrix_pitched_double_test 
#######################################################################################################################################################
#TEST28: MULTIPLY_MATRIX_VECTOR || COLUMNS = ROWS
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_double
cp gpu_matrix_vs_gpu_matrix_pitched_double_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 1;!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVector/CDouble.xls!" gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp


cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_double/gpu_matrix_vs_gpu_matrix_pitched_double_test 
#######################################################################################################################################################



