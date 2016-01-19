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

// $Id: pcg_test.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix_pitched.hpp"
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/jacobi.hpp"
#include "agile/operator/pcg.hpp"
#include "agile/network/gpu_communicator.hpp"

#include <iostream>

#define JACOBI_PRECONDITIONER 1

int main(int argc, char* argv[])
{
  // init the network
  agile::NetworkEnvironment environment(argc, argv);

  // allocate a GPU
  typedef agile::GPUCommunicator<unsigned, float, float> communicator_type;
  communicator_type com;
  com.allocateGPU();
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;

  // generate a matrix; NOTE: this matrix has to be symmetric positive definite
  // or there is no guarantee that the CG converges
  const unsigned SIZE = 20;
  float A_host[SIZE][SIZE];
  for (unsigned row = 0; row < SIZE; ++row)
    for (unsigned column = 0; column <= row; ++column)
    {
      A_host[row][column] = (float(SIZE) - float(row) + float(SIZE) / 2.0)
                            * (float(column) + 1.0);
      A_host[column][row] = A_host[row][column];
      if (row == column)
        A_host[row][column] = 2.0 * float(SIZE) + float(row) + float(column);
    }
  agile::GPUMatrixPitched<float> A(SIZE, SIZE, (float*)A_host);

  // generate a reference vector
  std::vector<float> x_reference_host(SIZE);
  for (unsigned counter = 0; counter < SIZE; ++counter)
    x_reference_host[counter] = float(SIZE) - float(counter) + float(SIZE/3);
  agile::GPUVector<float> x_reference;
  x_reference.assignFromHost(x_reference_host.begin(), x_reference_host.end());

  // generate a forward operator
  typedef agile::ForwardMatrix<communicator_type, agile::GPUMatrixPitched<float> >
    forward_type;
  forward_type forward(com, A);
#if JACOBI_PRECONDITIONER
  typedef agile::JacobiPreconditioner<communicator_type, float>
    preconditioner_type;
  std::vector<float> diagonal(SIZE);
  for (unsigned row = 0; row < SIZE; ++row)
    diagonal[row] = A_host[row][row];
  preconditioner_type preconditioner(com, diagonal);
#else
  typedef agile::InverseIdentity<communicator_type> preconditioner_type;
  preconditioner_type preconditioner(com);
#endif
  // generate a binary measure
  typedef agile::ScalarProductMeasure<communicator_type> measure_type;
  measure_type scalar_product(com);

  // generate the PCG solver
  // compute the inverse using PCG
  const double REL_TOLERANCE = 1e-12;
  const double ABS_TOLERANCE = 1e-6;
  const unsigned MAX_ITERATIONS = 100;
  agile::PreconditionedConjugateGradient<communicator_type, forward_type,
                                           preconditioner_type, measure_type>
   pcg(com, forward, preconditioner, scalar_product,
       REL_TOLERANCE, ABS_TOLERANCE, MAX_ITERATIONS);
  // generate the rhs
  agile::GPUVector<float> y(SIZE);
  forward(x_reference, y);
  // solve Ax = y using CG
  agile::GPUVector<float> x(SIZE);
  pcg(y, x);

  // some CG statistics
  if (pcg.convergence())
    std::cout << "CG converged in ";
  else
    std::cout << "Error: CG did not converge in ";
  std::cout << pcg.getIteration() + 1 << " iterations." << std::endl;
  std::cout << "Initial residual    = " << pcg.getRho0() << std::endl;
  std::cout << "Final residual      = " << pcg.getRho() << std::endl;
  std::cout << "Ratio rho_k / rho_0 = " << pcg.getRho() / pcg.getRho0()
            << std::endl;

  // output the reference
  std::cout << "Reference: " << std::endl << "  ";
  for (unsigned counter = 0; counter < x_reference_host.size(); ++counter)
    std::cout << x_reference_host[counter] << " ";
  std::cout << std::endl;
  // and the solution
  std::vector<float> x_host;
  x.copyToHost(x_host);
  std::cout << "CG solution: " << std::endl << "  ";
  for (unsigned counter = 0; counter < x_host.size(); ++counter)
    std::cout << x_host[counter] << " ";
  std::cout << std::endl;

  // calculate and output the difference
  agile::GPUVector<float> difference(SIZE);
  subVector(x_reference, x, difference);
  agile::GPUVector<float> difference_dist(difference);
  com.distribute(difference_dist);
  std::cout << "L2 of difference: "
            << std::sqrt(std::abs(scalar_product(difference, difference_dist)))
            << std::endl;

  return 0;
}

// End of $Id: pcg_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
