// Copyright (C) 2010-2011 Institute of Medical Engineering,
// Graz University of Technology
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses>.

// $Id: gridding_backward.cpp 476 2011-06-16 08:54:14Z freiberger $
#include "agile/gpu_environment.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_cs_matrix.hpp"
#include "agile/io/file.hpp"
#include <iostream>

// loading config stuff and helper
#include "gridding_functions.hpp"
#include "config.hpp"

int main(int argc, char* argv[])
{
  bool success;
  char *matrix_file = NULL;
  char *position_file = NULL;
  char *vector_file = NULL;
  unsigned gspace_width = 0;
  unsigned gspace_height = 0;
  char *result_file = NULL;

  // input checks
  if (argc <= 6)
  {
    std::cout << "usage: " << argv[0] << " <matrix file> <positon vector file> "
              << "<gspace vector file> <gspace width> <gspace height> "
              << " <result file>" << std::endl;
    return -1;
  } else {
    matrix_file = argv[1];
    position_file = argv[2];
    vector_file = argv[3];
    gspace_width = strtoul(argv[4], NULL, 0);
    gspace_height = strtoul(argv[5], NULL, 0);
    result_file = argv[6];
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

#if BACKWARD_GRIDDING_WITH_INTERPOLATION

  // Fixed types (e.g. positions assumed to be complex, real => x, imag => y)
  //--------------------------------------------------------------------------
  typedef std::complex<float> PositionDataType;

  std::vector<PositionDataType> pos_host;
  success = agile::readVectorFile(position_file, pos_host);
  if (!success)
  {
    std::cerr << "error: not able to load position vector: "
              << position_file << std::endl;
    exit(-1);
  }

#else

  // read in crs matrix
  //--------------------------------------------------------------------------
  unsigned A_num_rows, A_num_columns;
  std::vector<unsigned> A_row_nnz;
  std::vector<unsigned> A_column_index;
  std::vector<TRANSFORMATION_MATRIX_DATATYPE > A_data;
  success = agile::readCSMatrixFile(matrix_file,
                                      A_num_rows, A_num_columns,
                                      A_row_nnz, A_column_index,
                                      A_data);
  if (!success)
  {
    std::cerr << "error: not able to load matrix: " << matrix_file << std::endl;
    exit(-1);
  }

#endif

  // read in vector
  //--------------------------------------------------------------------------
  std::vector<DATA_VECTOR_DATATYPE > x_host;
  success = agile::readVectorFile(vector_file, x_host);
  if (!success)
  {
    std::cerr << "error: not able to load vector: " << vector_file << std::endl;
    exit(-1);
  }


  // dimension checks
  //--------------------------------------------------------------------------
#if !BACKWARD_GRIDDING_WITH_INTERPOLATION
  if (A_num_columns != x_host.size()) {
    std::cerr << "error: incompatible dimensions of matrix and "
              << "gspace vector" << std::endl;
    exit(-1);
  }
#endif

  if ((gspace_width * gspace_height) != x_host.size()) {
    std::cerr << "error: incompatible gspace vector length and "
              << "explicit given gspace dimensions" << std::endl;
    exit(-1);
  }

  // init gpu matrix and vector
  //--------------------------------------------------------------------------
#if BACKWARD_GRIDDING_WITH_INTERPOLATION
  typedef agile::GPUVector<PositionDataType> GPUPositionVectorType;
  GPUPositionVectorType pos(pos_host.size());
  pos.assignFromHost(pos_host.begin(), pos_host.end());
#else
  typedef agile::GPUCSMatrix<TRANSFORMATION_MATRIX_DATATYPE > GPUCSMatrixType;
  GPUCSMatrixType A(A_row_nnz, A_column_index, A_data);
#endif

  typedef agile::GPUVector<DATA_VECTOR_DATATYPE > GPUVectorType;
  GPUVectorType x(x_host.size());
  x.assignFromHost(x_host.begin(), x_host.end());

  // init result vector on gpu and cpu
  //--------------------------------------------------------------------------
#if BACKWARD_GRIDDING_WITH_INTERPOLATION
  GPUVectorType y(pos_host.size());
#else
  GPUVectorType y(A_num_rows);
#endif

  // backward gridding (grid data => radial data)
  //--------------------------------------------------------------------------
#if BACKWARD_GRIDDING_WITH_INTERPOLATION
  griddingBackwardInterp(communicator, x, gspace_width, gspace_height, pos, y);
#else
  griddingBackwardMult(communicator, A, x, y);
#endif

  // transfer result from gpu to cpu and write to file
  //--------------------------------------------------------------------------
  std::vector<DATA_VECTOR_DATATYPE > y_host;
  y.copyToHost(y_host);
  agile::writeVectorFile(result_file, y_host);

  return 0;
}

// End of $Id: gridding_backward.cpp 476 2011-06-16 08:54:14Z freiberger $.