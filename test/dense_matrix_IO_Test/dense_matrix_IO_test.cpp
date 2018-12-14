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

// $Id: dense_matrix_IO_test.cpp 476 2011-06-16 08:54:14Z freiberger $

// Read in a binary file that was written from Matlab, construct an
// identical copy of the matrix and write the result out again.
// Tests the matlab functions writeMatlab2bin.m and readBin2Matlab.m.
// The whole experiment is controlled with the scrip: dense_matrix_test.m

#include <fstream>
#include <iostream>
#include <vector>

#include "agile/gpu_environment.hpp"
#include "agile/gpu_matrix_pitched.hpp"
#include "agile/operator/fft2d.hpp"

// This is needed for image and kspace matrices
struct MatrixInfo
{
  unsigned num_rows;
  unsigned num_columns;
  unsigned num_bytes_per_entry;
  unsigned is_complex;
};

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " matrix.bin matrix_out.bin"
              << std::endl;
    return -1;
  }

  // acquire the device
  agile::GPUEnvironment::allocateGPU(0);
  agile::GPUEnvironment::printInformation(std::cout);

  // load the image matrix from the file
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

  if ((matrix_info.num_bytes_per_entry != 4))
  {
    std::cout << "Error: Only single precision is implemented." << std::endl;
    return -1;
  }

  // prepare the output file
  MatrixInfo matrix_out_info;
  matrix_out_info.num_rows = matrix_info.num_rows;
  matrix_out_info.num_columns = matrix_info.num_columns;
  matrix_out_info.num_bytes_per_entry = matrix_info.num_bytes_per_entry;
  matrix_out_info.is_complex = matrix_info.is_complex;
  
  unsigned matrix_out_file_size = matrix_out_info.num_rows 
                              * matrix_out_info.num_columns
                              * matrix_out_info.num_bytes_per_entry
                              * (matrix_out_info.is_complex ? 2 : 1);
  std::vector<char> matrix_out_buffer(matrix_out_file_size);
  std::ofstream result_file(argv[2], std::ofstream::binary);
  if (!result_file.is_open())
  {
    std::cout << "Error: Could not open '" << argv[2] << "' for output."
              << std::endl;
    return -1;
  }
  result_file.write((char*)&matrix_out_info, sizeof(matrix_out_info));

  // Do test for complex case
  if (matrix_out_info.is_complex == 1)
  {
    std::cout << "matrix_out_info.is_complex = " << matrix_out_info.is_complex
              << std::endl;
    std::cout << "You are running the test for a complex matrix." << std::endl;
  
    // Initialize GPUMatrixPitched Elements for original matrix and output matrix
    agile::GPUMatrixPitched<std::complex<float> > matrix(matrix_info.num_rows,
      matrix_info.num_columns, (std::complex<float>*)&matrix_buffer[0]);
    agile::GPUMatrixPitched<std::complex<float> > matrix_out(
      matrix_out_info.num_rows, matrix_out_info.num_columns,
      (std::complex<float>*)&matrix_out_buffer[0]);
    
    // Just make idetical copy: if you want to expand it, just do anything on 
     // the GPU you want here
    matrix_out = matrix;
    std::vector<std::complex<float> > matrix_out_host;
    matrix_out.copyToHost(matrix_out_host);
    
    const std::complex<float>* ptr = &matrix_out_host[0];
    for (unsigned row = 0; row < matrix_out.getNumRows(); ++row)
    {
      // do not write the zero-padded elements
      result_file.write((char*)ptr, matrix_out.getNumColumns() *
        sizeof(std::complex<float>));
      ptr += matrix_out.getPitchElements();
    }
  }
  // Do test for real case
  else
  {
    std::cout << "matrix_out_info.is_complex = " << matrix_out_info.is_complex
              << std::endl;
    std::cout << "You are running the test for a real matrix." << std::endl;
  
    // Initialize GPUMatrixPitched Elements for original matrix and output matrix
    agile::GPUMatrixPitched<float> matrix(matrix_info.num_rows,
      matrix_info.num_columns, (float*)&matrix_buffer[0]);
    agile::GPUMatrixPitched<float> matrix_out(matrix_out_info.num_rows,
      matrix_out_info.num_columns, (float*)&matrix_out_buffer[0]);
  
    // Just make idetical copy: if you want to expand it,
    // just do anything on the GPU you want here
    matrix_out = matrix;
    std::vector<float> matrix_out_host;
    matrix_out.copyToHost(matrix_out_host);
  
    const float* ptr = &matrix_out_host[0];
    for (unsigned row = 0; row < matrix_out.getNumRows(); ++row)
    {
      // do not write the zero-padded elements
      result_file.write((char*)ptr,
        matrix_out.getNumColumns() * sizeof(float));
      ptr += matrix_out.getPitchElements();
    }
  }
   
  result_file.close();
  return 0;
}

// End of $Id: dense_matrix_IO_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
