
#include "agile/gpu_environment.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_cs_matrix.hpp"
#include "agile/io/file.hpp"
#include "agile/gpu_timer.hpp"

int main(int argc, char* argv[])
{
  char *matrix_file = NULL;
  char *vector_file = NULL;
  char *result_file = NULL;

  // input checks
  if (argc <= 3)
  {
    std::cout << "usage: " << argv[0] << " <matrix file> <vector file>"
              << " <result file>" << std::endl;
    return -1;
  } else {
    matrix_file = argv[1];
    vector_file = argv[2];
    result_file = argv[3]; 
  }

  // init the network
  agile::NetworkEnvironment environment(argc, argv);

  // allocate a GPU
  typedef agile::GPUCommunicator<unsigned, float, float> communicator_type;
  communicator_type com;
  com.allocateGPU();

  bool success;

  typedef std::vector<std::complex<float> > cpu_vector_type;

  // read in crs matrix
  unsigned mat_num_rows, mat_num_columns;
  std::vector<unsigned> mat_row_nnz;
  std::vector<unsigned> mat_column_index;
  cpu_vector_type mat_data;
  success = agile::readCSMatrixFile(matrix_file,
                                      mat_num_rows, mat_num_columns,
                                      mat_row_nnz, mat_column_index,
                                      mat_data);
  if (!success)
  {
    std::cerr << "error: not able to load matrix: " << matrix_file << std::endl;
    exit(-1);
  }

  // read in vector
  cpu_vector_type x_host;
  success = agile::readVectorFile(vector_file, x_host);
  if (!success)
  {
    std::cerr << "error: not able to load vector: " << vector_file << std::endl;
    exit(-1);
  }

  // dimension check
  if (x_host.size() != mat_num_columns) {
    std::cerr << "error: incompatible dimensions" << std::endl;
  }

  // init gpu matrix and vector
  typedef agile::GPUCSMatrix<std::complex<float> > gpu_matrix_type;
  gpu_matrix_type A(mat_row_nnz, mat_column_index, mat_data);

  typedef agile::GPUVector<std::complex<float> > gpu_vector_type;
  gpu_vector_type x(mat_num_columns);
  x.assignFromHost(x_host.begin(), x_host.end());

  // init result vector on gpu
  gpu_vector_type y(mat_num_rows);

  // do the multiplication
  agile::tic();
  agile::multiply(A, x, y);
  agile::toc("multiply(A, x, y)");

  // transfer result from gpu to cpu
  cpu_vector_type y_host;
  y.copyToHost(y_host);

  // write result to file
  agile::writeVectorFile(result_file, y_host);
  return 0;
}
