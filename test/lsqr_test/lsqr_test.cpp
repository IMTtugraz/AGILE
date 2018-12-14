
#include "agile/gpu_environment.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix_pitched.hpp"
#include "agile/gpu_cs_matrix.hpp"
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/jacobi.hpp"
#include "agile/operator/lsqr.hpp"
#include "../test_defines.h"
#include <iostream>

int main(int argc, char* argv[])
{
  // init the network
  agile::NetworkEnvironment environment(argc, argv);

  // allocate a GPU
  typedef agile::GPUCommunicator<unsigned, float, float> communicator_type;
  communicator_type com;
  com.allocateGPU();

  // GPU Information
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;

  // generate a matrix; NOTE: this matrix has to be symmetric positive definite
  // or there is no guarantee that the LSQR converges
  const unsigned SIZE = 20;
  float A_host[SIZE][SIZE];
  for (unsigned row = 0; row < SIZE; ++row)
  {
    for (unsigned column = 0; column <= row; ++column)
    {
      A_host[row][column] = (float(SIZE) - float(row) + float(SIZE) / 2.0)
                            * (float(column) + 1.0);
      A_host[column][row] = A_host[row][column];
      if (row == column)
        A_host[row][column] = 2.0 * float(SIZE) + float(row) + float(column);
    }
  }

  // print out testmatrix (DEBUG)
  //-------------------------------------------------------------------
  for (unsigned row = 0; row < SIZE; ++row)
  {
    for (unsigned column = 0; column < SIZE; ++column)
      std::cout << A_host[row][column] << " ";      
    std::cout << std::endl;
  }
  std::cout << std::endl;
  //-------------------------------------------------------------------

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

  // generate a binary measure
  typedef agile::ScalarProductMeasure<communicator_type> measure_type;
  measure_type scalar_product(com);

  // generate the LSQR solver
  // compute the inverse using LSQR
  const double ABS_TOLERANCE = 1e-6;
  const unsigned MAX_ITERATIONS = 100;

  // init Operator
  agile::LSQR<communicator_type, forward_type, measure_type> lsqr(
    com, forward, scalar_product, ABS_TOLERANCE, MAX_ITERATIONS);

  // generate the rhs
  agile::GPUVector<float> y(SIZE);
  forward(x_reference, y);
  agile::GPUVector<float> x(SIZE);
  
  // solve Ax = y for x using LSQR
  lsqr(y, x);


  std::cout << std::endl;
  
  // some LSQR statistics
  if (lsqr.convergence())
    std::cout << "LSQR converged in ";
  else
    std::cout << "Error: LSQR did not converge in ";
  std::cout << lsqr.getIteration() + 1 << " iterations." << std::endl;
  std::cout << "Initial residual    = " << lsqr.getRho0() << std::endl;
  std::cout << "Final residual      = " << lsqr.getRho() << std::endl;
  std::cout << "Ratio rho_k / rho_0 = " << lsqr.getRho() / lsqr.getRho0()
            << std::endl;

  // rhs y
  std::vector<float> y_host;
  y.copyToHost(y_host);
  PRINT_VEC("y", y_host);

  // output the reference
  PRINT_VEC("Reference", x_reference_host);

  // and the solution
  std::vector<float> x_host;
  x.copyToHost(x_host);
  PRINT_VEC("LSQR solution for x", x_host);

  // calculate and output the difference
  agile::GPUVector<float> difference(SIZE);
  subVector(x_reference, x, difference);
  agile::GPUVector<float> difference_dist(difference);
  com.distribute(difference_dist);
  std::cout << std::endl
            << "L2 of difference: "
            << std::sqrt(std::abs(getScalarProduct(difference,
                                                   difference_dist)))
            << std::endl;
  return 0;
}
