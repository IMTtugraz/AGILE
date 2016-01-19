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

// $Id: gridding_forward.cpp 476 2011-06-16 08:54:14Z freiberger $

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
  char *A_matrix_file = NULL;
  char *AT_matrix_file = NULL;
  char *vector_file = NULL;
  char *result_file = NULL;

  // input checks
  if (argc <= 4)
  {
    std::cout << "usage: " << argv[0] << " <A matrix file> <A' matrix file> "
              << "<vector file> <result file>" << std::endl;
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
  CommunicatorType communicator;
  communicator.allocateGPU();

  // print out GPU information
  //agile::GPUEnvironment::printInformation(std::cout);

  bool success;

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

  // read in vector
  //--------------------------------------------------------------------------
  std::vector<DATA_VECTOR_DATATYPE > y_host;
  success = agile::readVectorFile(vector_file, y_host);
  if (!success)
  {
    std::cerr << "error: not able to load vector: " << vector_file << std::endl;
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
    std::cerr << "error: incompatible dimensions of matrix and vector"
              << std::endl;
    exit(-1);
  }

  // init gpu matrix and vector
  //--------------------------------------------------------------------------
  typedef agile::GPUCSMatrix<TRANSFORMATION_MATRIX_DATATYPE > GPUCSMatrixType;
  GPUCSMatrixType A(A_row_nnz, A_column_index, A_data);
  typedef agile::GPUCSMatrix<TRANSFORMATION_MATRIX_DATATYPE, true>
    GPUCSAdjointMatrixType;
  GPUCSAdjointMatrixType AT(AT_row_nnz, AT_column_index, AT_data);

  typedef agile::GPUVector<DATA_VECTOR_DATATYPE > GPUVectorType;
  GPUVectorType y(A_num_rows);
  y.assignFromHost(y_host.begin(), y_host.end());

  // init result vector on gpu
  //--------------------------------------------------------------------------
  GPUVectorType x(A_num_columns);

  // forward gridding (kspace => gspace)
  //--------------------------------------------------------------------------
  griddingForward(communicator, A, AT, y, x);

  // transfer result from gpu to cpu and write to file
  //--------------------------------------------------------------------------
  std::vector<DATA_VECTOR_DATATYPE > x_host;
  x.copyToHost(x_host);
  agile::writeVectorFile(result_file, x_host);

  return 0;
}


// End of $Id: gridding_forward.cpp 476 2011-06-16 08:54:14Z freiberger $.

