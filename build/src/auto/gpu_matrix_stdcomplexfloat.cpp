// $Id: gpu_matrix.cpp.in 450 2011-05-31 11:22:00Z freiberger $

/* This file was generated automatically by CMake. You have to modify '/home2/GIT/AGILE/src/gpu_matrix.cpp.in' if you want to make changes. */

#define TType1IsComplex 1
#define TType2IsComplex 

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>



#include "gpu_matrix.ipp"

namespace agile
{


  template void multiply<std::complex<float> >(
    const GPUMatrix<std::complex<float> >& A, const GPUVector<std::complex<float> >& x,
    GPUVector<std::complex<float> >& y);

  template void multiply<std::complex<float> >(
    const GPUVector<std::complex<float> >& x, const GPUMatrix<std::complex<float> >& A,
    GPUVector<std::complex<float> >& y);

  template void multiplyElementwise<std::complex<float> >(
    const GPUMatrix<std::complex<float> >& A, const GPUMatrix<std::complex<float> >& B,
    GPUMatrix<std::complex<float> >& Z);

  template to_real_type<std::complex<float> >::type norm2(
    const GPUMatrix<std::complex<float> >& A);

  template void copy<std::complex<float> >(const GPUMatrix<std::complex<float> >& A,
                                      GPUMatrix<std::complex<float> >& Z);
  //template void copy<std::complex<float> >(const GPUMatrix<std::complex<float> >& A,
  //                                    GPUMatrix<std::complex<float> >* Z);

  template void multiplyTransMatrixMatrix<std::complex<float> >(const GPUMatrix<std::complex<float> >& X, const GPUMatrix<std::complex<float> >& Y,
                GPUMatrix<std::complex<float> >& Z);

  template void dotProduct<std::complex<float> >(const GPUMatrix<std::complex<float> >& X, const GPUMatrix<std::complex<float> >& Y,
                                      std::complex<float> & val);

//#if !TType1IsComplex
//#if !TType2IsComplex

  template void scale<std::complex<float> >(
     const to_real_type<std::complex<float> >::type& alpha, const GPUMatrix<std::complex<float> >& A,
    GPUMatrix<std::complex<float> >& B);


//#endif  // !TType2IsComplex
//#endif  // !TType1IsComplex


} // namespace agile

// End of $Id: gpu_matrix.cpp.in 450 2011-05-31 11:22:00Z freiberger $.
