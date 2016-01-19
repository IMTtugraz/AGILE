
/* This file was generated automatically by CMake. You have to modify '/home2/GIT/AGILE/src/calc/l2solve.cpp.in' if you want to make changes. */

#define TType1IsComplex 1
#define TType2IsComplex 0

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

#include "calc/l2solve.cpp"

namespace agile
{
    template void L2Solve<std::complex<double>, double >::Init();


    //Pure Virtual Solve-Methode
    template void  L2Solve<std::complex<double> ,double >::Solve(const GPUMatrix<std::complex<double> >* u, const GPUMatrix<std::complex<double> >* c,
      const GPUMatrix<std::complex<double> >* rhs, const GPUMatrix<std::complex<double> >* u0,
      unsigned maxits, double alpha, double beta,
      GPUMatrix<std::complex<double> >* du, GPUMatrix<std::complex<double> >* dc);

} // namespace agile

