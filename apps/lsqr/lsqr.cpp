
#include "agile/gpu_environment.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/io/file.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_cs_matrix.hpp"
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/lsqr.hpp"

#include "config.hpp"

int main(int argc, char* argv[])
{
  char *A_matrix_file = NULL;
  char *AT_matrix_file = NULL;
  char *vector_file = NULL;
  char *result_file = NULL;

  // input checks
  if (argc <= 4)
  {
    std::cout << "usage: " << argv[0] << " <A matrix file> <AT matrix file>"
              << " <rhs vector file> <result file>" << std::endl;
    return -1;
  } else {
    A_matrix_file = argv[1];
    AT_matrix_file = argv[2];
    vector_file = argv[3];
    result_file = argv[4]; 
  }

  // init the network
  agile::NetworkEnvironment environment(argc, argv);
  
  // allocate a GPU
  typedef agile::GPUCommunicator<unsigned, float, float> CommunicatorType;
  CommunicatorType com;
  com.allocateGPU();

  bool success;

  typedef std::vector<std::complex<float> > cpu_vector_type;

  // read in crs matrix A
  //--------------------------------------------------------------------------
  unsigned A_num_rows, A_num_columns;
  std::vector<unsigned> A_row_nnz;
  std::vector<unsigned> A_column_index;
  cpu_vector_type A_data;
  success = agile::readCSMatrixFile(A_matrix_file,
                                      A_num_rows, A_num_columns,
                                      A_row_nnz, A_column_index,
                                      A_data);
  if (!success)
  {
    std::cerr << "error: not able to load matrix A: " << A_matrix_file << std::endl;
    exit(-1);
  }

  // read in crs matrix A'
  //--------------------------------------------------------------------------
  unsigned AT_num_rows, AT_num_columns;
  std::vector<unsigned> AT_row_nnz;
  std::vector<unsigned> AT_column_index;
  cpu_vector_type AT_data;
  success = agile::readCSMatrixFile(AT_matrix_file,
                                      AT_num_rows, AT_num_columns,
                                      AT_row_nnz, AT_column_index,
                                      AT_data);
  if (!success)
  {
    std::cerr << "error: not able to load matrix A': " << AT_matrix_file << std::endl;
    exit(-1);
  }

  // read in vector
  //--------------------------------------------------------------------------
  cpu_vector_type y_host;
  success = agile::readVectorFile(vector_file, y_host);
  if (!success)
  {
    std::cerr << "error: not able to load vector: " << vector_file << std::endl;
    exit(-1);
  }
  
  // dimension check
  //--------------------------------------------------------------------------
  if (A_num_rows != AT_num_columns || A_num_columns != AT_num_rows) {
    std::cerr << "error: incompatible matrix dimensions " << std::endl
              << "       A: " << A_num_rows << "x" << A_num_columns
              << ", AT: " << AT_num_rows << "x" << AT_num_columns << std::endl;
    exit(-1);
  }
  if (y_host.size() != A_num_rows) {
    std::cerr << "error: incompatible dimensions" << std::endl;
  }

  // init gpu matrix and vector
  typedef agile::GPUCSMatrix<std::complex<float> > GPUCSMatrixType;
  GPUCSMatrixType A(A_row_nnz, A_column_index, A_data);
  typedef agile::GPUCSMatrix<std::complex<float>, true>
    GPUCSAdjointMatrixType;
  GPUCSAdjointMatrixType AT(AT_row_nnz, AT_column_index, AT_data);

  typedef agile::GPUVector<std::complex<float> > gpu_vector_type;
  gpu_vector_type y(A_num_rows);
  y.assignFromHost(y_host.begin(), y_host.end());

  // generate a forward operator
  typedef agile::ForwardMatrixWithAdjoint<CommunicatorType, GPUCSMatrixType,
    GPUCSAdjointMatrixType> ForwardType;
  ForwardType forward(com, A, AT);

  // generate a binary measure
  typedef agile::ScalarProductMeasure<CommunicatorType> MeasureType;
  MeasureType scalar_product(com);

  // init lsqr operator
  agile::LSQR<CommunicatorType, ForwardType, MeasureType> lsqr(
    com, forward, scalar_product, LSQR_ABS_TOLERANCE, LSQR_MAX_ITERATIONS);

  // init result vector on gpu
  gpu_vector_type x(A_num_columns);

#if WITH_TIMER
  struct timeval st, et;
  gettimeofday(&st, NULL);
#endif

  // do lsqr inverse computation
  lsqr(y, x);

#if WITH_TIMER
  cudaThreadSynchronize();
  gettimeofday(&et, NULL);
  float elapsed_time = ((et.tv_sec-st.tv_sec)*1000.0
    + (et.tv_usec - st.tv_usec)/1000.0);
  std::cout << "lsqr (gpu):  " << std::setprecision(5.9)
            << elapsed_time << "ms" << std::endl;
#endif

#if SHOW_LSQR_DETAILS
  std::cout << "iterations: " << lsqr.getIteration() << std::endl;
  std::cout << "final residual: " << lsqr.getRho() << std::endl;
#endif

  // transfer result from gpu to cpu
  cpu_vector_type x_host;
  x.copyToHost(x_host);

  // write result to file
  agile::writeVectorFile(result_file, x_host);
  return 0;
}
