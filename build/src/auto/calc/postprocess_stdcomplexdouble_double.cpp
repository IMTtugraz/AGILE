
/* This file was generated automatically by CMake. You have to modify '/home/dieheart/workspace/AGILE/src/calc/postprocess.cpp.in' if you want to make changes. */

#define TType1IsComplex 1
#define TType2IsComplex 0

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

#include "calc/postprocess.cpp"

namespace agile
{
    template void PostProcess<std::complex<double> ,double >::calc_abs(const GPUMatrix<std::complex<double> >* in_mat,
            GPUMatrix<double >& out_mat);

    template void PostProcess<std::complex<double> ,double >::calc_phase(const GPUMatrix<std::complex<double> >* in_mat,
            GPUMatrix<double >& out_mat);

    template void PostProcess<std::complex<double> ,double >::calc_imag(const GPUMatrix<std::complex<double> >* in_mat,
            GPUMatrix<double >& out_mat);

    template void PostProcess<std::complex<double> ,double >::calc_real(const GPUMatrix<std::complex<double> >* in_mat,
            GPUMatrix<double >& out_mat);

    template void PostProcess<std::complex<double> ,double >::init();

} // namespace agile

