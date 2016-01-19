#include "agile/gpu_environment.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/io/file.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_cs_matrix.hpp"
#include <iostream>

// loading config stuff and helper
#include "gridding_functions.hpp"
#include "config.hpp"


int main(int argc, char* argv[])
{
  bool success;
  char *A_matrix_file = NULL;
  char *AT_matrix_file = NULL;
  char *position_file = NULL;
  char *vector_file = NULL;
  unsigned gspace_width = 0;
  unsigned gspace_height = 0;
  char *result_file = NULL;

  // input checks
  if (argc <= 7)
  {
    std::cout << "usage: " << argv[0] << " <A matrix file> <A' matrix file> "
              << "<positon vector file> <kspace vector file> "
              << "<gspace width> <gspace height> "
              << "<result gspace vector file>" << std::endl;
    return -1;
  } else {
    A_matrix_file = argv[1];
    AT_matrix_file = argv[2];
    position_file = argv[3];
    vector_file = argv[4];
    gspace_width = strtoul(argv[5], NULL, 0);
    gspace_height = strtoul(argv[6], NULL, 0);
    result_file = argv[7];
  }

  // init the network
  agile::NetworkEnvironment environment(argc, argv);

  // allocate a GPU
  typedef agile::GPUCommunicator<unsigned, float, float> CommunicatorType;
  CommunicatorType communicator;
  communicator.allocateGPU();

  // GPU Information
  //agile::GPUEnvironment::printInformation(std::cout);
  //std::cout << std::endl;

  // Fixed types (e.g. positionas assumed to be complex, real => x, imag => y)
  //--------------------------------------------------------------------------
  typedef std::complex<float> PositionDataType;


  // read in crs matrix A
  //--------------------------------------------------------------------------
  unsigned A_num_rows, A_num_columns;
  std::vector<unsigned> A_row_nnz;
  std::vector<unsigned> A_column_index;
  std::vector<TRANSFORMATION_MATRIX_DATATYPE > A_data;
  success = agile::readCSMatrixFile(A_matrix_file,
                                      A_num_rows, A_num_columns,
                                      A_row_nnz, A_column_index,
                                      A_data);
  if (!success)
  {
    std::cerr << "error: not able to load matrix A: " << A_matrix_file
              << std::endl;
    exit(-1);
  }

  // read in crs matrix A'
  //--------------------------------------------------------------------------
  unsigned AT_num_rows, AT_num_columns;
  std::vector<unsigned> AT_row_nnz;
  std::vector<unsigned> AT_column_index;
  std::vector<TRANSFORMATION_MATRIX_DATATYPE > AT_data;
  success = agile::readCSMatrixFile(AT_matrix_file,
                                      AT_num_rows, AT_num_columns,
                                      AT_row_nnz, AT_column_index,
                                      AT_data);
  if (!success)
  {
    std::cerr << "error: not able to load matrix A': " << AT_matrix_file
              << std::endl;
    exit(-1);
  }

  // read in kspace positon vector
  //--------------------------------------------------------------------------
  std::vector<PositionDataType> pos_host;
  success = agile::readVectorFile(position_file, pos_host);
  if (!success)
  {
    std::cerr << "error: not able to load position vector: "
              << position_file << std::endl;
    exit(-1);
  }

  // read in kspace vector
  //--------------------------------------------------------------------------
  std::vector<DATA_VECTOR_DATATYPE > y_host;
  success = agile::readVectorFile(vector_file, y_host);
  if (!success)
  {
    std::cerr << "error: not able to load kspace vector: "
              << vector_file << std::endl;
    exit(-1);
  }

  // dimension checks
  //--------------------------------------------------------------------------
  if (A_num_rows != AT_num_columns || A_num_columns != AT_num_rows) {
    std::cerr << "error: incompatible matrix dimensions " << std::endl
              << "       A: " << A_num_rows << "x" << A_num_columns
              << ", AT: " << AT_num_rows << "x" << AT_num_columns << std::endl;
    exit(-1);
  }
  if (y_host.size() != A_num_rows) {
    std::cerr << "error: incompatible dimensions of matrix and kspace vector"
              << std::endl;
    exit(-1);
  }
  if (pos_host.size() != A_num_rows) {
    std::cerr << "error: incompatible dimensions of matrix and position vector"
              << std::endl;
    exit(-1);
  }
  if ((gspace_width * gspace_height) != A_num_columns) {
    std::cerr << "error: incompatible dimensions of matrix and "
              << "explicit given gspace dimensions" << std::endl;
    exit(-1);
  }

  // init gpu matrix and vector
  //--------------------------------------------------------------------------
  typedef agile::GPUCSMatrix<TRANSFORMATION_MATRIX_DATATYPE > GPUCSMatrixType;
  GPUCSMatrixType A(A_row_nnz, A_column_index, A_data);

  typedef agile::GPUCSMatrix<TRANSFORMATION_MATRIX_DATATYPE, true>
    GPUCSAdjointMatrixType;
  GPUCSAdjointMatrixType AT(AT_row_nnz, AT_column_index, AT_data);

  typedef agile::GPUVector<PositionDataType> GPUPositionVectorType;
  GPUPositionVectorType pos(pos_host.size());
  pos.assignFromHost(pos_host.begin(), pos_host.end());

  typedef agile::GPUVector<DATA_VECTOR_DATATYPE > GPUVectorType;
  GPUVectorType y(y_host.size());
  y.assignFromHost(y_host.begin(), y_host.end());

  // init result vector on gpu
  //--------------------------------------------------------------------------
  GPUVectorType x(A_num_columns);

  // gridding computation loop
  //--------------------------------------------------------------------------

  // radial data => grid data
  griddingForward(communicator, A, AT, y, x);

  for (unsigned i=1; i<=LOOP_GRIDDING_NUMBER_OF_LOOPS; ++i)
  {
    std::cout << std::endl << "------------------------- loop (cpp) ";
    std::cout.fill('0');
    std::cout.width(2);
    std::cout << i;
    std::cout << " ---------------------------" << std::endl;

    // grid data => radial data
#if LOOP_BACKWARD_GRIDDING_WITH_INTERPOLATION
    griddingBackwardInterp(communicator,x,gspace_width,gspace_height,pos,y);
#else
    griddingBackwardMult(communicator, A, x, y);
#endif

#if LOOP_GRIDDING_RESET_X_TO_ZERO
    x.assign(x.size(), 0);
#endif

    // radial data => grid data
    griddingForward(communicator, A, AT, y, x);
  }

  std::cout << std::endl;

  // transfer result from gpu to cpu and write to file
  //--------------------------------------------------------------------------
  std::vector<DATA_VECTOR_DATATYPE > x_host;
  x.copyToHost(x_host);
  agile::writeVectorFile(result_file, x_host);

  return 0;
}
