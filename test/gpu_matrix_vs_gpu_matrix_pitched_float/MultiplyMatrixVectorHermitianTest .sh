#!/bin/bash

#MULTIPLY MATRIX VECTOR HERMITIAN

#ROWS = x * COLUMNS
#TEST1: MULTIPLY_MATRIX_VECTOR || ROWS = 2 * COLUMNS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_float
cp gpu_matrix_vs_gpu_matrix_pitched_float_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 2;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVectorHerm/R2C_CFloat.xls!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!num_rows = test_ctr;!num_rows = test_ctr*2;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_float/gpu_matrix_vs_gpu_matrix_pitched_float_test
#######################################################################################################################################################
#TEST2: MULTIPLY_MATRIX_VECTOR ||  ROWS = 5 * COLUMNS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_float
cp gpu_matrix_vs_gpu_matrix_pitched_float_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 2;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVectorHerm/R5C_CFloat.xls!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!num_rows = test_ctr;!num_rows = test_ctr*5;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_float/gpu_matrix_vs_gpu_matrix_pitched_float_test
#######################################################################################################################################################
#TEST3: MULTIPLY_MATRIX_VECTOR || ROWS = 10 * COLUMNS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_float
cp gpu_matrix_vs_gpu_matrix_pitched_float_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 2;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVectorHerm/R10C_CFloat.xls!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!num_rows = test_ctr;!num_rows = test_ctr*10;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_float/gpu_matrix_vs_gpu_matrix_pitched_float_test
#######################################################################################################################################################



#COLUMNS = x * ROWS


#######################################################################################################################################################
#TEST4: MULTIPLY_MATRIX_VECTOR || COLUMNS = 2 * ROWS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_float
cp gpu_matrix_vs_gpu_matrix_pitched_float_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 2;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVectorHerm/C2R_CFloat.xls!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!num_columns = test_ctr;!num_columns = test_ctr*2;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_float/gpu_matrix_vs_gpu_matrix_pitched_float_test
#######################################################################################################################################################
#TEST5: MULTIPLY_MATRIX_VECTOR || COLUMNS = 5 * ROWS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_float
cp gpu_matrix_vs_gpu_matrix_pitched_float_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 2;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVectorHerm/C5R_CFloat.xls!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!num_columns = test_ctr;!num_columns = test_ctr*5;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_float/gpu_matrix_vs_gpu_matrix_pitched_float_test
#######################################################################################################################################################
#TEST6: MULTIPLY_MATRIX_VECTOR || COLUMNS = 10 * ROWS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_float
cp gpu_matrix_vs_gpu_matrix_pitched_float_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 2;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVectorHerm/C10R_CFloat.xls!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!num_columns = test_ctr;!num_columns = test_ctr*10;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_float/gpu_matrix_vs_gpu_matrix_pitched_float_test
#######################################################################################################################################################
#TEST7: MULTIPLY_MATRIX_VECTOR || COLUMNS = ROWS 
cd ../test/gpu_matrix_vs_gpu_matrix_pitched_float
cp gpu_matrix_vs_gpu_matrix_pitched_float_edit.cpp gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp

sed -i -e "s!unsigned int test_nr = 0;!unsigned int test_nr = 2;!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp
sed -i -e "s!Timer Results!Timer Results 1404/MultiplyMatrixVectorHerm/CFloat.xls!" gpu_matrix_vs_gpu_matrix_pitched_float_test.cpp


cd ../../build_cuda32
make
test/gpu_matrix_vs_gpu_matrix_pitched_float/gpu_matrix_vs_gpu_matrix_pitched_float_test
#######################################################################################################################################################




