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


  template void multiply<float >(
    const GPUMatrix<float >& A, const GPUVector<float >& x,
    GPUVector<float >& y);

  template void multiply<float >(
    const GPUVector<float >& x, const GPUMatrix<float >& A,
    GPUVector<float >& y);

  template void multiplyElementwise<float >(
    const GPUMatrix<float >& A, const GPUMatrix<float >& B,
    GPUMatrix<float >& Z);

  template to_real_type<float >::type norm2(
    const GPUMatrix<float >& A);

  template void copy<float >(const GPUMatrix<float >& A,
                                      GPUMatrix<float >& Z);
  //template void copy<float >(const GPUMatrix<float >& A,
  //                                    GPUMatrix<float >* Z);

  template void multiplyTransMatrixMatrix<float >(const GPUMatrix<float >& X, const GPUMatrix<float >& Y,
                GPUMatrix<float >& Z);

  template void dotProduct<float >(const GPUMatrix<float >& X, const GPUMatrix<float >& Y,
                                      float & val);

//#if !TType1IsComplex
//#if !TType2IsComplex

  template void scale<float >(
     const to_real_type<float >::type& alpha, const GPUMatrix<float >& A,
    GPUMatrix<float >& B);


//#endif  // !TType2IsComplex
//#endif  // !TType1IsComplex


} // namespace agile

// End of $Id: gpu_matrix.cpp.in 450 2011-05-31 11:22:00Z freiberger $.
