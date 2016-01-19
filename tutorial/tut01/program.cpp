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

// \subsection headers Headers

// In the most basic version we need only two include files: one for the
// environment and another one for the GPU vector.
#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"

// \p iostream is no bad idea to have some output possibility.
#include <iostream>

// \subsection printing Printing function

// Define a small function that prints a vector to \p std::cout.
void output(const char* string, const std::vector<float>& x)
{
  std::cout << string;
  for (unsigned counter = 0; counter < x.size(); ++counter)
    std::cout << x[counter] << " ";
  std::cout << std::endl;
}

// \subsection main Main program.
int main()
{
  // \subsubsection init GPU initialization

  // The GPU is initialized using the singleton \p GPUEnvironment. This object
  // holds basic GPU information like the maximum number of threads, for
  // example. It is vital to initialize the environment before calling other
  // functions of the AGILE library. The allocation/initialization of the
  // GPU is achieved by calling \p allocateGPU(). We allocate GPU no. 0. If
  // you have more GPU's in your PC, try another number here.
  agile::GPUEnvironment::allocateGPU(0);
  
  // Now, we can print some information about our GPU to std::out.
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;

  // \subsubsection creation Creation of GPU vectors

  // The GPU is initialized now and we can start using it. First we create two
  // vectors on the host which will be transfered to the GPU later on.
  std::vector<float> x_host, y_host;
  for (unsigned counter = 0; counter < 10; ++counter)
  {
    x_host.push_back(counter * 2 + 1);
    y_host.push_back(counter * 2 + 2);
  }

  // We print the vectors to \p std::cout:
  output("x: ", x_host);
  output("y: ", y_host);
  
  // Now transfer these vectors to the GPU. You can use the \p assignFromHost
  // method to do so. This method takes two iterators \p begin and \p end
  // to a host vector and copies the elements in the range [\p begin, \p end).
  agile::GPUVector<float> x, y;
  x.assignFromHost(x_host.begin(), x_host.end());
  y.assignFromHost(y_host.begin(), y_host.end());

  // \subsubsection computations Basic vector operations

  // Now we can use the GPU vectors for calculations. Let's start by adding
  // them. The result shall be stored in another GPU vector \p z. \b NOTE: It
  // is really important to make sure that all the vectors have the correct
  // size because the library does not perform sanity checks!
  agile::GPUVector<float> z(x.size());
  addVector(x, y, z);
  
  // We want to print the result of \p z. As it is not possible to access
  // the GPU memory directly, we copy the whole vector back to the host.
  // You can use the \p copyToHost method for this task.
  std::vector<float> z_host;
  z.copyToHost(z_host);
  
  // Output the result.
  output("x + y: ", z_host);

  // We can also subtract two vectors.
  subVector(x, y, z);
  z.copyToHost(z_host);
  output("x - y: ", z_host);
  
  // Or we can add half of y to x:
  addScaledVector(x, float(0.5), y, z);
  z.copyToHost(z_host);
  output("x + 0.5 * y: ", z_host);

  // Also the scalar product is implemented. Its definition reads
  // \f$ \left ( x, y \right ) := \bar x^T y \f$.
  float scalar_product = getScalarProduct(x, y);
  std::cout << "(x, y) = " << scalar_product << std::endl;

  // Well, that concludes the first example.
  return 0;
}

// End of $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $.
