
/* This file was generated automatically by CMake. You have to modify '/home/dieheart/workspace/AGILE/src/calc/genkspacefov.cpp.in' if you want to make changes. */

#define TType1IsComplex 1

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

#include "calc/genkspacefov.cpp"

namespace agile
{
    template int KSpaceFOV<std::complex<float> >::genkspace_fov(
    const GPUMatrix<std::complex<float> >& in_mat, GPUMatrix<std::complex<float> >& out_mat);

} // namespace agile

