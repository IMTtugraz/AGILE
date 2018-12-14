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

// $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix_pitched.hpp"

#include <iostream>
#include <iomanip>

// Include the header for vector operations on the CPU.
#include "agile/cpu_vector.hpp"
// There is also a header implementing a dense matrix on the CPU.
#include "agile/cpu_matrix.hpp"

// Again our vector output function.
void output(const char* string, const agile::CPUVector<float>& x)
{
  std::cout << string;
  for (unsigned counter = 0; counter < x.size(); ++counter)
    std::cout << x[counter] << " ";
  std::cout << std::endl;
}

// The function to print a matrix.
void output(const char* string, unsigned num_rows, unsigned num_columns,
            agile::CPUVector<float> data)
{
  agile::CPUVector<float>::iterator iter = data.begin();
  std::cout << string << std::endl;
  for (unsigned row = 0; row < num_rows; ++row)
  {
    std::cout << "  ";
    for (unsigned column = 0; column < num_columns; ++column)
      std::cout << std::setw(4) << *iter++;
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

int main()
{
  // The initialization is the same as always.
  agile::GPUEnvironment::allocateGPU(0);
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;

  // We create two vectors on the host an output them.
  agile::CPUVector<float> x_host, y_host;
  for (unsigned counter = 0; counter < 10; ++counter)
  {
    x_host.push_back(counter * 2 + 1);
    y_host.push_back(counter * 2 + 2);
  }
  output("x: ", x_host);
  output("y: ", y_host);

  // Transfer the vectors to the GPU.
  agile::GPUVector<float> x, y;
  x.assignFromHost(x_host.begin(), x_host.end());
  y.assignFromHost(y_host.begin(), y_host.end());

  // As an example we perform the summation \f$ z \leftarrow x + y \f$.
  agile::GPUVector<float> z(x.size());
  agile::addVector(x, y, z);

  // Print this result.
  agile::CPUVector<float> z_host;
  z.copyToHost(z_host);
  output("GPU - x + y: ", z_host);

  // How do we know that the GPU implementation is correct? Well, simply do
  // the same computation on the CPU.
  agile::CPUVector<float> z_reference(x_host.size());
  agile::addVector(x_host, y_host, z_reference);
  output("CPU - x + y: ", z_reference);

  // Et voila... The results are identical. The header \p cpu_vector.hpp
  // provides all the functionality of \p gpu_vector.hpp for the CPU.
  // Have a look at the scalar product.
  std::cout << "GPU - (x, y) = " << agile::getScalarProduct(x, y)
            << std::endl;
  std::cout << "CPU - (x, y) = " << agile::getScalarProduct(x_host, y_host)
            << std::endl;

  // Let's also try the l2-norm of the x vectors.
  std::cout << "GPU - norm2(x) = " << norm2(x) << std::endl;
  std::cout << "CPU - norm2(x) = " << norm2(x_host) << std::endl << std::endl;

  // Now, we perform some test with a GPU matrix.
  // First, create a matrix and print it:
  agile::CPUVector<float> matrix_data;
  for (unsigned row = 0; row < x.size(); ++row)
    for (unsigned column = 0; column < x.size(); ++column)
      matrix_data.push_back(float(row + 1) * float(2 * column + 1));
  output("A: ", x.size(), x.size(), matrix_data);

  // Transfer the matrix to the GPU and perform a multiplication with \p x.
  agile::GPUMatrixPitched<float> A(x.size(), x.size(), &matrix_data[0]);
  agile::multiply(A, x, z);
  z.copyToHost(z_host);
  output("GPU - A * x: ", z_host);

  // Do a reference calculation on the CPU.
  agile::CPUMatrix<float> A_host(x_host.size(), x_host.size(), &matrix_data[0]);
  agile::multiply(A_host, x_host, z_reference);
  output("CPU - A * x: ", z_reference);

  // The transposed product should also bring the same result.
  agile::multiply(x, A, z);
  z.copyToHost(z_host);
  output("GPU - A^H * x: ", z_host);
  agile::multiply(x_host, A_host, z_reference);
  output("CPU - A^H * x: ", z_reference);

  // This tutorial demonstrated how to make reference solutions on the CPU.
  return 0;
}

// End of $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $.
