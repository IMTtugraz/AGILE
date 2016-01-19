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

// $Id: gpu_matrix_pitched_test.cpp 476 2011-06-16 08:54:14Z freiberger $

#include <iostream>

#include "agile/gpu_environment.hpp"
#include "agile/gpu_matrix_pitched.hpp"

int main()
{
  // acquire the device and print some specs
  agile::GPUEnvironment::allocateGPU(0);
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;

  // create the matrix
  float matrixA[5][5] = {{1, 0, 3, 4, 7},
                         {8, 2, 0, 0, 0},
                         {0, 0, 0, 7, 2},
                         {0, 0, 0, 1, 0},
                         {1, 9, 4, 0, 0}};
  agile::GPUMatrixPitched<float> A(5, 5, (float*) matrixA);

  // print testmatrix (DEBUG)
  //-------------------------------------------------------------------
  std::cout << "full test matrix A:" << std::endl;
  for (unsigned row = 0; row < 5; ++row)
  {
    for (unsigned column = 0; column < 5; ++column)
      std::cout << " " << matrixA[row][column];
    std::cout << ";" << std::endl;
  }
  //-------------------------------------------------------------------

  // create vector for multiplication
  std::vector<float> x_host(5);
  for (unsigned counter = 0; counter < 5; ++counter)
    x_host[counter] = (counter + 1) * (counter + 1);
  agile::GPUVector<float> x;
  x.assignFromHost(x_host.begin(), x_host.end());

  // perform multiplication
  agile::GPUVector<float> y(5);
  agile::multiply(A, x, y);
  // copy the result to the host
  std::vector<float> y_host;
  y.copyToHost(y_host);
  for (unsigned counter = 0; counter < y_host.size(); ++counter)
    std::cout << y_host[counter] << " ";
  std::cout << std::endl;

  // create a complex matrix
  std::complex<float> matrixB[5][5];
  for (unsigned row = 0; row < 5; ++row)
    for (unsigned column = 0; column < 5; ++column)
      matrixB[row][column]
        = std::complex<float>(matrixA[row][column],
                              (float(row) + 1) * (float(column) - 3));
  agile::GPUMatrixPitched<std::complex<float> > B(5, 5,
                                             (std::complex<float>*) matrixB);

  // print testmatrix (DEBUG)
  //-------------------------------------------------------------------
  std::cout << "full test matrix B:" << std::endl;
  for (unsigned row = 0; row < 5; ++row)
  {
    for (unsigned column = 0; column < 5; ++column)
      std::cout << " " << real(matrixB[row][column])
                 << "+i*" << imag(matrixB[row][column]);
    std::cout << ";" << std::endl;
  }
  //-------------------------------------------------------------------

  // complex multiplication
  agile::GPUVector<std::complex<float> > z(5);
  agile::multiply(B, x, z);
  // copy the result to the host
  std::vector<std::complex<float> > z_host;
  z.copyToHost(z_host);
  for (unsigned counter = 0; counter < z_host.size(); ++counter)
    std::cout << z_host[counter] << " ";
  std::cout << std::endl;

  // complex multiplication the other way round
  agile::GPUVector<std::complex<float> > u(5);
  std::vector<std::complex<float> > u_host(5);
  for (unsigned counter = 0; counter < u_host.size(); ++counter)
    u_host[counter] = std::complex<float>(
                        float(counter) * 3 - 4, -float(counter) * counter);
  u.assignFromHost(u_host.begin(), u_host.end());
  agile::multiply(A, u, z);
  // copy the result to the host
  z.copyToHost(z_host);
  for (unsigned counter = 0; counter < z_host.size(); ++counter)
    std::cout << z_host[counter] << " ";
  std::cout << std::endl;

  // double complex multiplication
  agile::multiply(B, u, z);
  // copy the result to the host
  z.copyToHost(z_host);
  for (unsigned counter = 0; counter < z_host.size(); ++counter)
    std::cout << z_host[counter] << " ";
  std::cout << std::endl;

  return 0;
}

// End of $Id: gpu_matrix_pitched_test.cpp 476 2011-06-16 08:54:14Z freiberger $.

