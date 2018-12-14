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

// $Id: fft_test_kspace.cpp 476 2011-06-16 08:54:14Z freiberger $

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
    std::cout << "Usage: " << argv[0] << " img.bin img_recon.bin"
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
  MatrixInfo img_info;
  file.read((char*)&img_info, sizeof(img_info));

  // create a temporary buffer to hold the matrix
  unsigned img_file_size = img_info.num_rows * img_info.num_columns
                              * img_info.num_bytes_per_entry
                              * (img_info.is_complex ? 2 : 1);
  std::vector<char> img_buffer(img_file_size);
  file.read(&img_buffer[0], img_file_size);
  file.close();

  if ((img_info.num_bytes_per_entry != 4))
  {
    std::cout << "Error: Only single precision is implemented." << std::endl;
    return -1;
  }

  // Prepare the temporary kspace (not sure if this is really needed)
  MatrixInfo kspace_info;
  kspace_info.num_rows = img_info.num_rows;
  kspace_info.num_columns = img_info.num_columns;
  kspace_info.num_bytes_per_entry = img_info.num_bytes_per_entry;
  kspace_info.is_complex = 1;

  unsigned kspace_file_size = kspace_info.num_rows * kspace_info.num_columns
                              * kspace_info.num_bytes_per_entry
                              * (kspace_info.is_complex ? 2 : 1);

  std::vector<char> kspace_buffer(kspace_file_size);

  // prepare the output kspace file
  MatrixInfo img_recon_info;
  img_recon_info.num_rows = img_info.num_rows;
  img_recon_info.num_columns = img_info.num_columns;
  img_recon_info.num_bytes_per_entry = img_info.num_bytes_per_entry;
  img_recon_info.is_complex = 1;

  unsigned img_recon_file_size = img_recon_info.num_rows * img_recon_info.num_columns
                              * img_recon_info.num_bytes_per_entry
                              * (img_recon_info.is_complex ? 2 : 1);

  std::vector<char> img_recon_buffer(img_recon_file_size);

  std::ofstream result_file(argv[2], std::ofstream::binary);
  if (!result_file.is_open())
  {
    std::cout << "Error: Could not open '" << argv[2] << "' for output."
              << std::endl;
    return -1;
  }
  result_file.write((char*)&img_recon_info, sizeof(img_recon_info));

  clock_t start_time, end_time;
  clock_t fft_start_time, fft_end_time;

  // do the fft
  start_time = clock();
  // Initialize GPUMatrixPitched Elements for the original image, kspace, and the retransformed image
  agile::GPUMatrixPitched<std::complex<float> > img(img_info.num_rows, img_info.num_columns,
			      (std::complex<float>*)&img_buffer[0]);
  agile::GPUMatrixPitched<std::complex<float> > kspace(kspace_info.num_rows, kspace_info.num_columns,
			      (std::complex<float>*)&kspace_buffer[0]);
  agile::GPUMatrixPitched<std::complex<float> > img_recon(img_recon_info.num_rows, img_recon_info.num_columns,
			      (std::complex<float>*)&img_recon_buffer[0]);

  // generate the fft operator
  agile::FFT2D<std::complex<float>, std::complex<float> > cudafft2;

  // Transforms happen here
  fft_start_time = clock();
  fftshift(img);
  cudafft2.forward(img,kspace);
  fftshift(kspace);
  fft_end_time = clock();

  // copy the result to the host
  std::vector<std::complex<float> > kspace_host;
  kspace.copyToHost(kspace_host);
  end_time = clock();

  const std::complex<float>* ptr = &kspace_host[0];
  for (unsigned row = 0; row < kspace.getNumRows(); ++row)
  {
    // do not write the zero-padded elements
    result_file.write((char*)ptr, kspace.getNumColumns() * sizeof(std::complex<float>));
    ptr += kspace.getPitchElements();
  }

  result_file.close();

  std::cout << "Time needed incl data transfer (ms): "
            << double(end_time - start_time) * 1000.0 / CLOCKS_PER_SEC
            << std::endl;
  std::cout << "Time needed for multiplication (ms): "
            << double(fft_end_time - fft_start_time) * 1000.0 / CLOCKS_PER_SEC
            << std::endl;
  std::cout << std::endl << std::endl;

  return 0;
}

// End of $Id: fft_test_kspace.cpp 476 2011-06-16 08:54:14Z freiberger $.
