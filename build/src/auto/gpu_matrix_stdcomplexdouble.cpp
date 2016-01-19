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


  template void multiply<std::complex<double> >(
    const GPUMatrix<std::complex<double> >& A, const GPUVector<std::complex<double> >& x,
    GPUVector<std::complex<double> >& y);

  template void multiply<std::complex<double> >(
    const GPUVector<std::complex<double> >& x, const GPUMatrix<std::complex<double> >& A,
    GPUVector<std::complex<double> >& y);

  template void multiplyElementwise<std::complex<double> >(
    const GPUMatrix<std::complex<double> >& A, const GPUMatrix<std::complex<double> >& B,
    GPUMatrix<std::complex<double> >& Z);

  template to_real_type<std::complex<double> >::type norm2(
    const GPUMatrix<std::complex<double> >& A);

  template void copy<std::complex<double> >(const GPUMatrix<std::complex<double> >& A,
                                      GPUMatrix<std::complex<double> >& Z);
  //template void copy<std::complex<double> >(const GPUMatrix<std::complex<double> >& A,
  //                                    GPUMatrix<std::complex<double> >* Z);

  template void multiplyTransMatrixMatrix<std::complex<double> >(const GPUMatrix<std::complex<double> >& X, const GPUMatrix<std::complex<double> >& Y,
                GPUMatrix<std::complex<double> >& Z);

  template void dotProduct<std::complex<double> >(const GPUMatrix<std::complex<double> >& X, const GPUMatrix<std::complex<double> >& Y,
                                      std::complex<double> & val);

//#if !TType1IsComplex
//#if !TType2IsComplex

  template void scale<std::complex<double> >(
     const to_real_type<std::complex<double> >::type& alpha, const GPUMatrix<std::complex<double> >& A,
    GPUMatrix<std::complex<double> >& B);


//#endif  // !TType2IsComplex
//#endif  // !TType1IsComplex


} // namespace agile

// End of $Id: gpu_matrix.cpp.in 450 2011-05-31 11:22:00Z freiberger $.
