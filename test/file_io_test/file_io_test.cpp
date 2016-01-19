
#include "agile/gpu_environment.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_cs_matrix.hpp"
#include "agile/io/file.hpp"
#include "../test_defines.h"

#include <iostream>

// max vector, matrix dimension to display
#define DISP_DIM_MAX 20


#define TEST_VECTOR_1 "./data/x_complex_2.bin"
#define TEST_VECTOR_2 "./data/x_real_2.bin"
#define TEST_VECTOR_3 "./data/x_complex_64x64_uniform.bin"

#define TEST_MATRIX_1 "./data/A_complex_2x2.bin"
#define TEST_MATRIX_2 "./data/A_real_2x2.bin"
#define TEST_MATRIX_3 "./data/A_real_64x64_x_75x64.bin"

#define TEST_VECTOR TEST_VECTOR_3
#define TEST_MATRIX TEST_MATRIX_3


int main(int argc, char* argv[])
{
  char *matrix_file = NULL;
  char *vector_file = NULL;
  if (argc <= 2)
  {
    std::cout << "usage: " << argv[0] << " matrix.bin vector.bin" << std::endl;
    return -1;
  } else {
    matrix_file = argv[1];
    vector_file = argv[2];
  }
  
  // init the network
  agile::NetworkEnvironment environment(argc, argv);
  
  // allocate a GPU
  typedef agile::GPUCommunicator<unsigned, float, float> communicator_type;
  communicator_type com;
  com.allocateGPU();

  // GPU Information
  //agile::GPUEnvironment::printInformation(std::cout);
  //std::cout << std::endl;

  bool success;

  // read in crs matrix A
  unsigned mat_num_rows, mat_num_columns;
  std::vector<unsigned> mat_row_nnz;
  std::vector<unsigned> mat_column_index;
  std::vector<std::complex<float> > mat_data;
  success = agile::readCSMatrixFile(matrix_file,
                                      mat_num_rows, mat_num_columns,
                                      mat_row_nnz, mat_column_index,
                                      mat_data);
  if (!success)
  {
    std::cerr << "Not able to load from file: " << matrix_file << std::endl;
    exit(-1);
  }
  
  if (mat_num_rows > DISP_DIM_MAX || mat_num_columns > DISP_DIM_MAX)
    std::cout << "size of A: " << mat_num_rows << "x" << mat_num_columns
              << std::endl;
  else
  {
    PRINT_VEC("A (row_nnz)", mat_row_nnz);
    PRINT_VEC("A (column_index)", mat_column_index);
    PRINT_VEC("A (data)", mat_data);
  }

  PRINT_DELIMITER();
  
  // read in vector x
  std::vector<std::complex<float> > vec_data;
  success = agile::readVectorFile(vector_file, vec_data);
  if (!success)
  {
    std::cerr << "Not able to load from file: " << vector_file << std::endl;
    exit(-1);
  }
  
  if (vec_data.size() > DISP_DIM_MAX)
    std::cout << "size of x: " << vec_data.size() << std::endl;
  else
    PRINT_VEC("x", vec_data);
  
  
  // multiplication test from here
  PRINT_DELIMITER();
  
  // dimension check
  if (vec_data.size() != mat_num_columns) {
    std::cerr << "Error: incompatible dimensions" << std::endl;
  }
  
  // init gpu matrix and vector
  agile::GPUVector<std::complex<float> > x(vec_data.size());
  x.assignFromHost(vec_data.begin(), vec_data.end());
  
  agile::GPUCSMatrix<std::complex<float> > A(mat_row_nnz,
                                               mat_column_index,
                                               mat_data);

  // init result vector on gpu
  agile::GPUVector<std::complex<float> > y(mat_num_rows);
  agile::multiply(A, x, y);

  // transfer result from gpu to cpu
  std::vector<std::complex<float> > y_host;
  y.copyToHost(y_host);
  
  if (y_host.size() > DISP_DIM_MAX)
    std::cout << "size of y (y=Ax): " << y_host.size() << std::endl;
  else
    PRINT_VEC("y = Ax", y_host);
  
  
  PRINT_DELIMITER();
  
  agile::writeVectorFile("result.bin", y_host);
  std::cout << "matrix, vector IO test finished successfully!"
            << std::endl << std::endl;
            
  return 0;
}
