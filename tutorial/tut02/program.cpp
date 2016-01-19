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

// We have to include the headers for the environment, for the GPU vector and
// for the GPU matrix.
#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"

#include "agile/gpu_matrix_pitched.hpp"


#include <iostream>
#include <iomanip>

// Define a small function that prints a vector to \p std::cout.
void output(const char* string, std::vector<float> x)
{
  std::cout << string;
  for (unsigned counter = 0; counter < x.size(); ++counter)
    std::cout << x[counter] << " ";
  std::cout << std::endl;
}

// Another function to print a matrix.
void output(const char* string, unsigned num_rows, unsigned num_columns,
            std::vector<float> data)
{
  std::vector<float>::iterator iter = data.begin();
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

// Here starts the main program.
int main()
{
  // Initialize the first GPU and print information as done in the first step.
  agile::GPUEnvironment::allocateGPU(0);
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;
  
  // Create a vector first on the CPU and transfer it to the GPU.
  std::vector<float> x_host;
  for (unsigned counter = 0; counter < 10; ++counter)
    x_host.push_back(counter * 2 + 1);
  output("x: ", x_host);
  agile::GPUVector<float> x;
  x.assignFromHost(x_host.begin(), x_host.end());

  // Now we create a dense matrix. The elements are stored in row-major order.
  std::vector<float> matrix_data;
  for (unsigned row = 0; row < x.size(); ++row)
    for (unsigned column = 0; column < x.size(); ++column)
      matrix_data.push_back(float(row + 1) * float(2 * column + 1));


  // Print the matrix to \p std::cout.
  output("A: ", x.size(), x.size(), matrix_data);


  // We transfer the matrix to the GPU. This can be done using the constructor
  // of \p GPUMatrixPitched, which takes the number of rows, the number of columns
  // and a pointer to an array of size (rows * columns) holding the matrix
  // elements.
  agile::GPUMatrixPitched<float> A(x.size(), x.size(), &matrix_data[0]);


  // We need another vector to store the result of our matrix vector
  // multiplications. NOTE: This matrix has to have the correct dimensions
  // because the library is too lazy to check this!
  agile::GPUVector<float> y(x.size());


  // The hard stuff is done. Now we can use our matrix. Perform the product
  // \f$ y \leftarrow Ax \f$.
  multiply(A, x, y);

  // Transfer the result back to the host and print it.
  std::vector<float> y_host;
  y.copyToHost(y_host);
  output("A * x: ", y_host);

  // Also the multiplication with the hermitian matrix \f$ A^H = \bar A^T \f$
  // is implemented. It can be evaluated by changing the order of arguments
  // to the \p multiply function (i.e. vector-in, matrix-in, vector-out).
  multiply(x, A, y);

  // Output this result, too.
  y.copyToHost(y_host);
  output("A^H * x: ", y_host);



  // Hopefully this tutorial made clear how to use GPU vectors together with
  // GPU matrices.
  return 0;
}

// End of $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $.
