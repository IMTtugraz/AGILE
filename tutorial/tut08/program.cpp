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

// We start with including the important stuff.
#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix_pitched.hpp"
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/jacobi.hpp"
#include "agile/operator/pcg.hpp"
#include "agile/network/gpu_communicator.hpp"

#include <iostream>

// The next macro can be used to switch between a Jacobi (diagonal)
// or an identity preconditioner. More on that later...
#define JACOBI_PRECONDITIONER 1

int main(int argc, char* argv[])
{
  // Initialization of the network, the communicator and the allocation of
  // the GPU is done as in previous tutorials.
  agile::NetworkEnvironment environment(argc, argv);
  typedef agile::GPUCommunicator<unsigned, float, float> communicator_type;
  communicator_type com;
  com.allocateGPU();
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;

  // We are interested in solving the linear problem \f$ Ax = y \f$, with a
  // given matrix \f$ A \f$ and a right-hand side vector \f$ y \f$. The unknown
  // is the vector \f$ x \f$.
  // Now, we can generate a matrix that shall be inverted (actually we do not
  // invert the matrix but use the CG algorithm). Note that CG requires a
  // symmetric positive definite (SPD) matrix and it is not too trivial to
  // write down a SPD matrix. If you fail to provide a SPD matrix to the CG
  // algorithm there is no guarantee that it will converge. You might be lucky,
  // you might be not...
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

  // The matrix is still in the host's memory and has to be transfered to the
  // GPU. This is done automatically by the constructor of \p GPUMatrixPitched.
  agile::GPUMatrixPitched<float> A(SIZE, SIZE, (float*)A_host);

  // Next we need a reference solution. We can create any vector we like at
  // this place.
  std::vector<float> x_reference_host(SIZE);
  for (unsigned counter = 0; counter < SIZE; ++counter)
    x_reference_host[counter] = float(SIZE) - float(counter) + float(SIZE/3);

  // This vector has to be transfered to the GPU memory too. For vectors, this
  // can be achieved by the member function \p assignFromHost.
  agile::GPUVector<float> x_reference;
  x_reference.assignFromHost(x_reference_host.begin(), x_reference_host.end());

  // We wrap the GPU matrix from above into a forward operator called
  // \p ForwardMatrix. Forward operators are simply objects that implement
  // the parenthesis-operator \p operator() which takes an
  // \p accumulated vector and returns a \p distributed one. In all other
  // respects the operator is a black box for us.
  // The \p ForwardMatrix operator requires a reference to the communicator
  // when constructing the object so that it has access to the network.
  typedef agile::ForwardMatrix<communicator_type, agile::GPUMatrixPitched<float> >
    forward_type;
  forward_type forward(com, A);

  // What we also want to use a preconditioner, which means that we change from
  // the original problem \f$ Ax = y \f$ to the equivalent one
  // \f$ PAx = Py \f$, where \f$ P \f$ is a preconditioner. The rationale is
  // that most often the matrix \f$ A \f$ is ill-conditioned and the CG algorithm
  // does not converge properly at all or it needs many iterations. The use of
  // a preconditioner makes the whole system better conditioned. The simplest
  // choice is to use the identity \f$ P = I \f$ (which means no preconditioning
  // at all). The best choice would be \f$ P = A^{-1} \f$ as we would have the
  // solution for \f$ x \f$ in the first step already (but then we need again
  // to find the inverse of \f$ A \f$ which we wanted to avoid). An
  // 'intermediate' possibility is to take \f$ P = diag(A)^{-1} \f$ which is
  // easy and fast to invert and gives better results than the identity.
  // A preconditioner belongs to the inverse operator. All inverse operators
  // implement a parenthesis-operator which takes a \p distributed vector
  // as input and returns an \p accumulated one (opposite to the forward
  // operators, thus).
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

  // The last operator needed is a measure. A measure operator has again
  // a parenthesis-operator. This timeis takes an \p accumulated vector as first
  // input and a \p distributed one as second input and returns a scalar
  // measuring somehow the size of the vectors. An example is the scalar
  // product operator.
  typedef agile::ScalarProductMeasure<communicator_type> measure_type;
  measure_type scalar_product(com);

  // Finally, generate the PCG solver. It needs the absolute and relative
  // tolerances as input so that it knows when the solution is good enough for
  // our purposes. Furthermore it requires the maximum amount of iterations
  // after which it simply capitulates without having found a solution.
  const double REL_TOLERANCE = 1e-12;
  const double ABS_TOLERANCE = 1e-6;
  const unsigned MAX_ITERATIONS = 100;
  agile::PreconditionedConjugateGradient<communicator_type, forward_type,
                                           preconditioner_type, measure_type>
    pcg(com, forward, preconditioner, scalar_product,
        REL_TOLERANCE, ABS_TOLERANCE, MAX_ITERATIONS);

  // What we have not generated, yet, is the right hand side \f$ y \f$. This is
  // simply one call to our forward operator.
  agile::GPUVector<float> y(SIZE);
  forward(x_reference, y);

  // We need one more vector to hold the result of the CG algorithm. Note that
  // we also supply the initial guess for the solution via this vector.
  agile::GPUVector<float> x(SIZE);

  // Finally, we have constructed, initialized, wrapped... everything. The only
  // thing left to do is to call the CG operator.
  pcg(y, x);

  // Print some statistics (and hope that the operator actually converged).
  if (pcg.convergence())
    std::cout << "CG converged in ";
  else
    std::cout << "Error: CG did not converge in ";
  std::cout << pcg.getIteration() + 1 << " iterations." << std::endl;
  std::cout << "Initial residual    = " << pcg.getRho0() << std::endl;
  std::cout << "Final residual      = " << pcg.getRho() << std::endl;
  std::cout << "Ratio rho_k / rho_0 = " << pcg.getRho() / pcg.getRho0()
            << std::endl;

  // As the vectors in this example were quite small we can even print them to
  // standard output.
  std::cout << "Reference: " << std::endl << "  ";
  for (unsigned counter = 0; counter < x_reference_host.size(); ++counter)
    std::cout << x_reference_host[counter] << " ";
  std::cout << std::endl;

  // The solution is still on the GPU and has to be transfered to the CPU memory.
  // This is accomplished using \p copyToHost.
  std::vector<float> x_host;
  x.copyToHost(x_host);

  // Output the solution, too.
  std::cout << "CG solution: " << std::endl << "  ";
  for (unsigned counter = 0; counter < x_host.size(); ++counter)
    std::cout << x_host[counter] << " ";
  std::cout << std::endl;

  // Finally, we also compute the difference between the reference solution and
  // the true solution (of course, we do this on the GPU).
  agile::GPUVector<float> difference(SIZE);
  subVector(x_reference, x, difference);

  // To measure the distance, we use the scalar product measure we have
  // introduced above. Note, that this operator wants the first vector in
  // accumulated format and the second one in distributed format. The solution
  // we got from the CG algorithm is accumulated (because CG is an inverse
  // operator). This means, we have to distribute the solution to have mixed
  // formats.
  agile::GPUVector<float> difference_dist(difference);
  com.distribute(difference_dist);
  std::cout << "L2 of difference: "
            << std::sqrt(std::abs(scalar_product(difference, difference_dist)))
            << std::endl;

  // So, that's it.
  return 0;
}

// End of $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $.
