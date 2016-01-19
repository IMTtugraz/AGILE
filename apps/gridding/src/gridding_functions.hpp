// $Id: gridding_functions.hpp 476 2011-06-16 08:54:14Z freiberger $
#ifndef _GRIDDING_FUNCTIONS_HPP_
#define _GRIDDING_FUNCTIONS_HPP_

#include "agile/operator/forward_operator.hpp"
#include "agile/operator/binary_measure_operator.hpp"
#include "agile/operator/lsqr.hpp"

// loading config stuff
#include "config.hpp"


#if WITH_TIMER
  #include <iomanip>
  #include <sys/time.h>
#endif


template <typename CommunicatorType, typename GPUCSMatrixType,
  typename GPUCSAdjointMatrixType, typename GPUVectorType>
void griddingForward(CommunicatorType &communicator, GPUCSMatrixType &A,
  GPUCSAdjointMatrixType& AT, GPUVectorType &y, GPUVectorType &x)
{
#if WITH_TIMER
  struct timeval st, et;
#endif

  // generate a forward operator
  //--------------------------------------------------------------------------
  //typedef agile::ForwardMatrix<communicator_type, gpu_matrix_type>
  //  forward_type;
  //forward_type forward(com, A);

  // optimized version with predefined adjoint matrix
  typedef agile::ForwardMatrixWithAdjoint<CommunicatorType, GPUCSMatrixType,
    GPUCSAdjointMatrixType> ForwardType;
  ForwardType forward(communicator, A, AT);

  // generate a binary measure
  typedef agile::ScalarProductMeasure<CommunicatorType> MeasureType;
  MeasureType scalar_product(communicator);

  // generate the LSQR solver
  // init lsqr operator
  agile::LSQR<CommunicatorType, ForwardType, MeasureType>
    lsqr(communicator, forward, scalar_product, LSQR_ABS_TOLERANCE,
      LSQR_MAX_ITERATIONS);

#if WITH_TIMER
  gettimeofday(&st, NULL);
#endif

  // do lsqr inverse computation
  lsqr(y, x);

#if WITH_TIMER
  cudaThreadSynchronize();
  gettimeofday(&et, NULL);
  float elapsed_time = ((et.tv_sec-st.tv_sec)*1000.0
    + (et.tv_usec - st.tv_usec)/1000.0);
  std::cout << "forward-time (lsqr):";
  std::cout.fill(' ');
  std::cout.width(45);
  std::cout << std::right << elapsed_time << "ms" << std::endl;
#endif

#if WITH_LSQR_DETAILS
  std::cout << "iterations: " << lsqr.getIteration()
            << " (max. iter:" << LSQR_MAX_ITERATIONS << ", "
            << "tol. " << LSQR_ABS_TOLERANCE << ")" << std::endl
            << "initial residual: " << std::setprecision(9) << lsqr.getRho0()
            << std::endl
            << "final residual: " << std::setprecision(9) << lsqr.getRho()
            << std::endl;
#endif
}



template <typename CommunicatorType, typename GPUCSMatrixType,
  typename GPUVectorType>
void griddingBackwardMult(CommunicatorType &communicator, GPUCSMatrixType &A,
  GPUVectorType &x, GPUVectorType &y)
{
#if WITH_TIMER
  struct timeval st, et;
  gettimeofday(&st, NULL);
#endif

  agile::multiply(A, x, y);

#if WITH_TIMER
  cudaThreadSynchronize();
  gettimeofday(&et, NULL);
  float elapsed_time = ((et.tv_sec-st.tv_sec)*1000.0
    + (et.tv_usec - st.tv_usec)/1000.0);
  std::cout << "backward-time (multiplication):";
  std::cout.fill(' ');
  std::cout.width(34);
  std::cout << std::right << elapsed_time << "ms" << std::endl;
#endif
}


template <typename CommunicatorType, typename GPUVectorType,
  typename GPUPositionVectorType>
void griddingBackwardInterp(CommunicatorType &communicator, GPUVectorType &x,
  unsigned gspace_width, unsigned gspace_height, GPUPositionVectorType &pos,
  GPUVectorType &y)
{
#if WITH_TIMER
  struct timeval st, et;
  gettimeofday(&st, NULL);
#endif

  // x is a matrix serialized into the vector x by matlab with "matrix(:)"
  // and matlab is using column major ordering.
  agile::interpolate2d(x, gspace_width, gspace_height, false, pos, y);

#if WITH_TIMER
  cudaThreadSynchronize();
  gettimeofday(&et, NULL);
  float elapsed_time = ((et.tv_sec-st.tv_sec)*1000.0
    + (et.tv_usec - st.tv_usec)/1000.0);
  std::cout << "backward-time (interpolation):";
  std::cout.fill(' ');
  std::cout.width(35);
  std::cout << std::right << elapsed_time << "ms" << std::endl;
#endif
}

#endif

// End of $Id: gridding_functions.hpp 476 2011-06-16 08:54:14Z freiberger $.
