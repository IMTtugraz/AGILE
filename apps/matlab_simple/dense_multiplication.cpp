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

// $Id: dense_multiplication.cpp 476 2011-06-16 08:54:14Z freiberger $

#include <fstream>
#include <iostream>
#include <vector>

#include "agile/gpu_environment.hpp"
#include "agile/gpu_matrix_pitched.hpp"

struct MatrixInfo
{
  unsigned num_rows;
  unsigned num_columns;
  unsigned num_bytes_per_entry;
  unsigned is_complex;
};

struct VectorInfo
{
  unsigned size;
  unsigned num_bytes_per_entry;
  unsigned is_complex;
};

int main(int argc, char* argv[])
{
  if (argc < 4)
  {
    std::cout << "Usage: " << argv[0] << " matrix.bin x_vector.bin out.bin"
              << std::endl;
    return -1;
  }

  // acquire the device
  agile::GPUEnvironment::allocateGPU(0);
  agile::GPUEnvironment::printInformation(std::cout);

  // load the matrix from the file
  std::ifstream file(argv[1], std::ifstream::binary);
  if (!file.is_open())
  {
    std::cout << "Error: Could not open '" << argv[1] << "'." << std::endl;
    return -1;
  }
  // read the amount of rows and columns as well as the data type's size
  MatrixInfo matrix_info;
  file.read((char*)&matrix_info, sizeof(matrix_info));
  // create a temporary buffer to hold the matrix
  unsigned matrix_file_size = matrix_info.num_rows * matrix_info.num_columns
                              * matrix_info.num_bytes_per_entry
                              * (matrix_info.is_complex ? 2 : 1);
  std::vector<char> matrix_buffer(matrix_file_size);
  file.read(&matrix_buffer[0], matrix_file_size);
  file.close();

  // load the vector
  file.open(argv[2], std::ifstream::binary);
  if (!file.is_open())
  {
    std::cout << "Error: Could not open '" << argv[2] << "'." << std::endl;
    return -1;
  }
  // read the vector info
  VectorInfo x_info;
  file.read((char*)&x_info, sizeof(x_info));
  // again a temporary buffer to hold the data
  unsigned x_file_size = x_info.size * x_info.num_bytes_per_entry
                         * (x_info.is_complex ? 2 : 1);
  std::vector<char> x_buffer(x_file_size);
  file.read(&x_buffer[0], x_file_size);
  file.close();

  if (matrix_info.num_bytes_per_entry != x_info.num_bytes_per_entry)
  {
    std::cout << "Error: Matrix and vector must be of the same type."
              << std::endl;
    return -1;
  }

  if ((matrix_info.num_bytes_per_entry != 4)
      || (x_info.num_bytes_per_entry != 4))
  {
    std::cout << "Error: Only single precision is implemented." << std::endl;
    return -1;
  }

  // prepare the output file
  VectorInfo y_info;
  y_info.size = matrix_info.num_rows;
  y_info.num_bytes_per_entry = matrix_info.num_bytes_per_entry;
  y_info.is_complex = matrix_info.is_complex || x_info.is_complex;
  std::ofstream result_file(argv[3], std::ofstream::binary);
  if (!result_file.is_open())
  {
    std::cout << "Error: Could not open '" << argv[3] << "' for output."
              << std::endl;
    return -1;
  }
  result_file.write((char*)&y_info, sizeof(y_info));

  clock_t start_time, end_time;
  clock_t mult_start_time, mult_end_time;

  // do the multiplication
  if (!matrix_info.is_complex && !x_info.is_complex)
  {
    start_time = clock();
    agile::GPUMatrixPitched<float> A(matrix_info.num_rows, matrix_info.num_columns,
                                (float*)&matrix_buffer[0]);
    agile::GPUVector<float> x;
    x.assignFromHost((float*)&x_buffer[0], (float*)&x_buffer[0] + x_info.size);
    agile::GPUVector<float> y(matrix_info.num_rows);
    mult_start_time = clock();
    agile::multiply(A, x, y);
    mult_end_time = clock();
    // copy the result to the host
    std::vector<float> y_host;
    y.copyToHost(y_host);
    end_time = clock();
    result_file.write((char*)&y_host[0],
                      y_info.size * y_info.num_bytes_per_entry);
  }
  else
  {
    // either one is complex and so will be the result
    start_time = clock();
    agile::GPUVector<std::complex<float> > y(matrix_info.num_rows);
    if (matrix_info.is_complex)
    {
      // complex matrix
      agile::GPUMatrixPitched<std::complex<float> > A(
        matrix_info.num_rows, matrix_info.num_columns,
        (std::complex<float>*)&matrix_buffer[0]);
      if (x_info.is_complex)
      {
        agile::GPUVector<std::complex<float> > x;
        x.assignFromHost((std::complex<float>*)&x_buffer[0],
                         (std::complex<float>*)&x_buffer[0] + x_info.size);
        mult_start_time = clock();
        agile::multiply(A, x, y);
        mult_end_time = clock();
      }
      else
      {
        agile::GPUVector<float> x;
        x.assignFromHost((float*)&x_buffer[0], (float*)&x_buffer[0] + x_info.size);
        mult_start_time = clock();
        agile::multiply(A, x, y);
        mult_end_time = clock();
      }
    }
    else
    {
      // matrix is real so the x-vector is complex
      agile::GPUMatrixPitched<float> A(matrix_info.num_rows, matrix_info.num_columns,
                                  (float*)&matrix_buffer[0]);
      agile::GPUVector<std::complex<float> > x;
        x.assignFromHost((std::complex<float>*)&x_buffer[0],
                         (std::complex<float>*)&x_buffer[0] + x_info.size);
      mult_start_time = clock();
      agile::multiply(A, x, y);
      mult_end_time = clock();
    }
    // copy the result to the host
    std::vector<std::complex<float> > y_host;
    y.copyToHost(y_host);
    end_time = clock();
    result_file.write((char*)&y_host[0],
                      y_info.size * y_info.num_bytes_per_entry * 2);
  }
  result_file.close();

  std::cout << "Time needed incl data transfer (ms): "
            << double(end_time - start_time) * 1000.0 / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "Time needed for multiplication (ms): "
            << double(mult_end_time - mult_start_time) * 1000.0 / CLOCKS_PER_SEC
            << std::endl;
  std::cout << std::endl << std::endl;

  return 0;
}

// End of $Id: dense_multiplication.cpp 476 2011-06-16 08:54:14Z freiberger $.
