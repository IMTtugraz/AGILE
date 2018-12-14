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

// $Id: 3D_matrix_I_O.cpp 476 2011-06-16 08:54:14Z freiberger $

#include <fstream>
#include <iostream>
#include <vector>

#include "agile/gpu_environment.hpp"
#include "agile/gpu_matrix_pitched.hpp"

// This is needed for image matrices
struct MatrixInfo
{
  unsigned num_rows;
  unsigned num_columns;
  unsigned num_coils;
  unsigned num_bytes_per_entry;
  unsigned is_complex;
};

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " img.bin img_output.bin"
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

  // read the amount of rows, columns and coils as well as the data type's size
  MatrixInfo img_info;
  file.read((char*)&img_info, sizeof(img_info));

  // create a temporary buffer to hold the matrix
  unsigned img_file_size = img_info.num_rows * img_info.num_columns
			      * img_info.num_coils
                              * img_info.num_bytes_per_entry
                              * (img_info.is_complex ? 2 : 1);
  std::vector<char> img_buffer(img_file_size);
  file.read(&img_buffer[0], img_file_size);
  file.close();
  
  //std::cout << "img_info.is_complex: " << img_info.is_complex << std::endl;
  
  if ((img_info.num_bytes_per_entry != 4))
  {
    std::cout << "Error: Only single precision is implemented." << std::endl;
    return -1;
  }

  // prepare the output image file
  MatrixInfo img_output_info;
  img_output_info.num_rows = img_info.num_rows;
  img_output_info.num_columns = img_info.num_columns;
  img_output_info.num_coils = img_info.num_coils;
  img_output_info.num_bytes_per_entry = img_info.num_bytes_per_entry;
  img_output_info.is_complex = img_info.is_complex;

  unsigned img_output_file_size = img_output_info.num_rows * img_output_info.num_columns
                              * img_output_info.num_coils
			      * img_output_info.num_bytes_per_entry
                              * (img_output_info.is_complex ? 2 : 1);

  std::vector<char> img_output_buffer(img_output_file_size);
  
  std::cout << "img_output_info.num_rows: " << img_output_info.num_rows << std::endl;

  std::ofstream result_file(argv[2], std::ofstream::binary);
  if (!result_file.is_open())
  {
    std::cout << "Error: Could not open '" << argv[2] << "' for output."
              << std::endl;
    return -1;
  }
  result_file.write((char*)&img_output_info, sizeof(img_output_info));

  clock_t start_time, end_time;
  clock_t assignement_start_time, assignement_end_time;

  // GPU work starts here
  start_time = clock();
  // Initialize GPUMatrixPitched Elements for the original image and the output image
  agile::GPUMatrixPitched<std::complex<float> > img(img_info.num_rows*img_info.num_coils, img_info.num_columns,
			      (std::complex<float>*)&img_buffer[0]);
  agile::GPUMatrixPitched<std::complex<float> > img_output(img_output_info.num_rows*img_info.num_coils, img_output_info.num_columns,
			      (std::complex<float>*)&img_output_buffer[0]);

  // Calculations happen here: This is a proxy, nothing is done
  assignement_start_time = clock();
  img_output = img;
  assignement_end_time = clock();

  // copy the result to the host
  std::vector<std::complex<float> > img_output_host;
  img_output.copyToHost(img_output_host);
  end_time = clock();

  std::cout << "img_output.getNumRows: " << img_output.getNumRows() << std::endl;
  
  const std::complex<float>* ptr = &img_output_host[0];
  for (unsigned row = 0; row < img_output.getNumRows(); ++row)
  {
    // do not write the zero-padded elements
    result_file.write((char*)ptr, img_output.getNumColumns() * sizeof(std::complex<float>));
    ptr += img_output.getNumColumns();
  }

  result_file.close();

  std::cout << "Time needed incl data transfer (ms): "
            << double(end_time - start_time) * 1000.0 / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "Time needed for multiplication (ms): "
            << double(assignement_end_time - assignement_start_time) * 1000.0 / CLOCKS_PER_SEC
            << std::endl;
  std::cout << std::endl << std::endl;

  return 0;
}

// End of $Id: 3D_matrix_I_O.cpp 476 2011-06-16 08:54:14Z freiberger $.
