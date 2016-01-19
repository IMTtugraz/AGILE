
/* This file was generated automatically by CMake. You have to modify '/home2/GIT/AGILE/src/calc/tvsolve.cpp.in' if you want to make changes. */

#define TType1IsComplex 1
#define TType2IsComplex 0

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

#include "calc/tvsolve.cpp"

namespace agile
{
    template void TVSolve<std::complex<float>, float >::Init();


    //Pure Virtual Solve-Methode
    template void  TVSolve<std::complex<float> ,float >::Solve(const GPUMatrix<std::complex<float> >* u, const GPUMatrix<std::complex<float> >* c,
      const GPUMatrix<std::complex<float> >* rhs, const GPUMatrix<std::complex<float> >* u0,
      unsigned maxits, float alpha, float beta,
      GPUMatrix<std::complex<float> >* du, GPUMatrix<std::complex<float> >* dc);

} // namespace agile

