// $Id: gpu_matrix.cpp.in 450 2011-05-31 11:22:00Z freiberger $

/* This file was generated automatically by CMake. You have to modify '/home2/GIT/AGILE/src/gpu_matrix.cpp.in' if you want to make changes. */

#define TType1IsComplex 0
#define TType2IsComplex 

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>



#include "gpu_matrix.ipp"

namespace agile
{


  template void multiply<double >(
    const GPUMatrix<double >& A, const GPUVector<double >& x,
    GPUVector<double >& y);

  template void multiply<double >(
    const GPUVector<double >& x, const GPUMatrix<double >& A,
    GPUVector<double >& y);

  template void multiplyElementwise<double >(
    const GPUMatrix<double >& A, const GPUMatrix<double >& B,
    GPUMatrix<double >& Z);

  template to_real_type<double >::type norm2(
    const GPUMatrix<double >& A);

  template void copy<double >(const GPUMatrix<double >& A,
                                      GPUMatrix<double >& Z);
  //template void copy<double >(const GPUMatrix<double >& A,
  //                                    GPUMatrix<double >* Z);

  template void multiplyTransMatrixMatrix<double >(const GPUMatrix<double >& X, const GPUMatrix<double >& Y,
                GPUMatrix<double >& Z);

  template void dotProduct<double >(const GPUMatrix<double >& X, const GPUMatrix<double >& Y,
                                      double & val);

//#if !TType1IsComplex
//#if !TType2IsComplex

  template void scale<double >(
     const to_real_type<double >::type& alpha, const GPUMatrix<double >& A,
    GPUMatrix<double >& B);


//#endif  // !TType2IsComplex
//#endif  // !TType1IsComplex


} // namespace agile

// End of $Id: gpu_matrix.cpp.in 450 2011-05-31 11:22:00Z freiberger $.
